const mathLib = require("./mathLib.js");

const { randwalk_perturbation } = require("./mif2cRW.js");
const { rprocessInternal } = require("./rprocessInternal.js");
const { dmeasureInternal } = require("./dmeasureInternal.js");
const { pfilter_computations } = require("./pfilterComputations.js");
const { cooling, partrans, coef } = require("./mif2Helpers.js");

exports.mif2Internal = function (arg) {
  let pomp = arg.object;
  let Nmif = arg.Nmif;
  let start = arg.start;
  let transform = arg.transform? arg.transform: FALSE;
  let ivps =arg.ivps;
  let pars= arg.pars;
  let rw_sd = arg.rw_sd;
  let Np = arg.Np;
  let factorvar = arg.factorvar;
  let coolingType = arg.coolingType ? arg.coolingType: ["hyperbolic", "geometric"];
  let coolingFraction = arg.coolingFraction;
  let tol = arg.tol? arg.tol: 1e-17;
  let maxFail = arg.maxFail? arg.maxFail: Infinity;
  let verbose = arg.verbose? argv: false;
  let _paramMatrix = arg._paramMatrix? arg._paramMatrix: null;
  let _ndone = arg._ndone? arg._ndone: 0;
  let _indices = arg._indices? arg._indices: 0;
  if (Array.isArray(Nmif) || !isFinite(Nmif) || Nmif < 1)
    throw new Error("Nmif must be a positive integer.");
    
  Nmif = parseInt(Nmif);
  if (_paramMatrix === null) {
    if (!start) start = coef(pomp);
    // if (Object.getPrototypeOf(start) == Object.prototype) start = Object.values(start);
  } else { 
    throw new Error("Not translated, if paramMatrix is supplied");  
  }
  
  if (Object.keys(start).length === 0 || !Object.values(start).every(element => {return typeof element === 'number'}))
    throw new Error("parameters must be specified as a named numeric vector.")

  let ntimes = pomp.times.length;
  if (typeof Np === "function" || Np === undefined || Np <= 0) {
    throw new Error(`Number of particles should be a positive number. ${Np} is not translated`);
  }

  if (!rw_sd) throw new Error("rw_sd function must be specified!");
  // pkern.sd: it produces the matrix of params rw in ntime 
  let rw_sd_matrix = [];
  for (let i = 0; i < ntimes; i++) {
    rw_sd_matrix.push(rw_sd(pomp.times[i]))
  }

  if (Array.isArray(coolingFraction) || !isFinite(coolingFraction) || coolingFraction <= 0 || coolingFraction > 1)
    throw new Error(`coolingFraction must be in (0,1]. coolingFraction = ${coolingFraction}`);
  let coolingFn = cooling(coolingType, coolingFraction, ntimes);
  let paramMatrix;
  if (_paramMatrix === null) {
    paramMatrix = new Array(Np).fill(null).map(a => start);
  } else {
    paramMatrix = _paramMatrix;
  }

  convRec = new Array(Nmif + 1).fill(null).map(a => new Object());
  convRec[0] = Object.assign({loglik: null, nfail: null}, start)
  if (transform)
    paramMatrix = partrans(pomp, paramMatrix, dir="toEstimationScale");
  
  // Iterate the filtering main loop 
  let pfp;
  for (let n = 0; n < Nmif; n++) {
    try {
      pfp = mif2Pfilter(
        object=pomp,
        params=paramMatrix,
        Np=Np,
        mifiter=_ndone + n + 1,
        coolingFn,
        rw_sd= rw_sd_matrix,
        tol=tol,
        maxFail,
        transform,
        _indices=_indices,
      )
      
    } catch (error) {
      console.error(`Iterate the filtering stoped: ${error}`)
    }
    paramMatrix = pfp.paramMatrix;
    convRec[n].loglik = pfp.loglik;
    convRec[n].nfail = pfp.nfail;

    convRec[n + 1].loglik = null;
    convRec[n + 1].nfail = null;
    Object.assign(convRec[n + 1], coef(pfp));
    
    _indices = pfp.indices;
  }
  
  if (transform)
    pfp.paramMatrix = partrans(pomp,paramMatrix,dir="fromEstimationScale");
  
  
  return{
    ...pfp,
    Nmif: Nmif,
    rw_sd: rw_sd,
    coolingType: coolingType,
    coolingFraction: coolingFraction,
    transform: transform,
    convRec: convRec
  }
}

const mif2Pfilter = function (object, params, Np, mifiter, coolingFn, rw_sd,
  tol = 1e-17, maxFail = Inf, transform, _indices = 0)
{
  if ((Array.isArray(tol) && tol.length !== 1) || tol === Infinity || tol < 0)
    throw new Error(`${tol} should be a small positive number in mif2Pfilter.`);

  let do_ta = !!_indices.length ;
  if (do_ta)
    throw new Error(` ${_indices} has improper length in mif2Pfilter.`);

  let times = [object.t0, ...object.times];
  let ntimes = times.length - 1;
  let loglik = new Array(ntimes);
  let effSampleSize = Number(ntimes);
  let nfail = 0;
  let alpha, pmag;
  let x = [];

  for (let nt = 0; nt < ntimes; nt++) {//ntimes
    alpha = coolingFn(nt + 1,mifiter).alpha;
    pmg = Object.assign({}, rw_sd[nt]);
    Object.keys(pmg).map(key => pmg[key] *= alpha);
    params = randwalk_perturbation(params, pmg);
    
    if (transform)
      tparams = partrans(object, params, dir="fromEstimationScale");
    
    if (nt === 0) {
      //get initial states
      let initparams = transform ? tparams : params;
      for (let i = 0; i < params.length; i++) {
        x.push(object.initializer(initparams[i], object.interpolator(object.t0)));
      }
    } 
    
    let X = [];
    try {
      X = rprocessInternal(
        object,
        xstart = x,
        [times[nt],times[nt + 1]],
        transform ? tparams : params,
        offset=1
      );
    } catch (error) {
      console.error(`In mif2.js: process simulation error: ${error}`);
    }
    
    let weights = [];
    try {
      weights = dmeasureInternal(
        object,
        y = object.data[nt],
        X,
        times[nt + 1],
        transform ? tparams : params,
        log = false
      ); 
    } catch (error) {
      console.error(`In mif2.js: error in calculation of weights: ${error}`);
    }
    
    let allFinite = weights.map(w => isFinite(w)).reduce((a, b) => a & b, 1);
    if (!allFinite) {
      throw new Error("In dmeasure: weights returns non-finite value");
    }

    // compute weighted mean at last timestep
    if (nt === ntimes - 1) {
      if (weights.map(w => w>0).reduce((a, b) => a || b, 0)) {
        // replace and fill object.params instead of coef(object). This is the same thing.
        object.params = mathLib.mean(params, w = weights);
        object.params = coef(object, transform = transform);
      } else {
        console.warn("filtering failure at last filter iteration, using unweighted mean for 'coef' ");
        object.params = mathLib.mean(params);
        object.params = coef(object, transform = transform);
      }
    }

    // compute effective sample size, log-likelihood
    // also do resampling if filtering has not failed
    let xx;
    try {
      xx  = pfilter_computations(
        X,
        params,
        Np = Np,
        0, //rw_sd
        predmean = false,
        predvar = false,
        filtmean = false,
        trackancestry = do_ta,
        onepar = false,
        weights = weights,
        tol = tol
      );
    } catch (error) {
      console.error(`particle-filter error: ${error}`) 
    }

    let allFail = xx.fail;
    loglik[nt] = xx.loglik;
    effSampleSize[nt] = xx.ess;
    if (do_ta) {
      _indices = _indices[xx.ancestry];
    }

    x = xx.states;
    params = xx.params;// should be in toScale.

    if (allFail) { // all particles are lost
      nfail = nfail + 1;
      if (nfail > maxFail)
        throw new Error("In mif2Pfilter: too many filtering failures")
    }    
  } // end of nt loop   
    

  if (nfail > 0) {
    console.log("warning! filtering failure occurred.");
  }

  return {
    ...object,
    paramMatrix: params,
    effSamplesize: effSampleSize,
    condLoglik: loglik,
    indices: _indices,
    Np: Np,
    tol: tol,
    nfail: Number(nfail),
    loglik: loglik.reduce((a,b) => a+b, 0)
  }

}


