const mathLib = require("./mathLib.js");

const { randwalk_perturbation } = require("./mif2cRW.js");
const { rprocessInternal } = require("./rprocessInternal.js");
const { dmeasureInternal } = require("./dmeasureInternal.js");
const { pfilter_computations } = require("./pfilterComputations.js");
const { cooling, partrans } = require("./mif2Helpers.js");

exports.mif2Internal = function (object)
  {
  let pomp = object.pomp;
  let Nmif = object.Nmif;
  let start = object.start;
  let transform = object.transform? object.transform: FALSE;
  let ivps =object.ivps;
  let pars= object.pars;
  let rw_sd = object.rw_sd;
  let param_rwIndex = object.param_rwIndex;
  let Np = object.Np;
  let factorvar = object.factorvar;
  let coolingType = object.coolingType ? object.coolingType: ["hyperbolic", "geometric"];
  let coolingFraction = object.coolingFraction;
  let tol = object.tol? object.tol: 1e-17;
  let maxFail = object.maxFail? object.maxFail: Infinity;
  let verbose = object.verbose? objectv: false;
  let _paramMatrix = object._paramMatrix? object._paramMatrix: null;
  let _ndone = object._ndone? object._ndone: 0;
  let _indices = object._indices? object._indices: 0;
  if (Array.isArray(Nmif) || !isFinite(Nmif) || Nmif < 1)
    throw new Error("Nmif must be a positive integer.");
    
  Nmif = parseInt(Nmif);
  if (_paramMatrix === null) {
    if (!start) start = pomp.coef;
    if (Object.getPrototypeOf(start) == Object.prototype) start = Object.values(start);
  } else { 
    throw new Error("Not translated, if paramMatrix is supplied");  
  }
  if (start.length ===0 || !start.every(element => {return typeof element === 'number'}))
    throw new Error("parameters must be specified as a named numeric vector.")

  let ntimes = pomp.times.length;
  if (typeof Np === "function" || Np === undefined || Np <= 0) {
    throw new Error(`Number of particles should be a positive number. ${Np} is not translated`);
  }

  if (!rw_sd) throw new Error("rw_sd function must be specified!");
  //TODO: pkern.sd: it produce the matrix of params rw in ntime 
  let rw_sd_matrix = [];
  for (let i = 0; i < ntimes; i++) {
    rw_sd_matrix.push(rw_sd(pomp.times[i]))
  }

  if (Array.isArray(coolingFraction) || !isFinite(coolingFraction) || coolingFraction <= 0 || coolingFraction > 1)
    throw new Error(`coolingFraction must be in (0,1]. coolingFraction = ${coolingFraction}`);
  let coolingFn = cooling(coolingType, coolingFraction, ntimes);
  let paramMatrix;
  if (_paramMatrix === null) {
    paramMatrix = new Array(Np).fill(start);
  } else {
    paramMatrix = _paramMatrix;
  }

  convRec = new Array(Nmif + 1).fill(Array(start.length + 2));
  convRec[1] = [null, null, ...start];//[loglik,nfail,...start]
  if (transform)
    paramMatrix = partrans(pomp,paramMatrix,dir="toEstimationScale");
  
  // Iterate the filtering main loop 
  let pfp;
  for (let n = 0; n < Nmif; n++) {
    // try {
      pfp = mif2Pfilter(
        object=pomp,
        params=paramMatrix,
        Np=Np,
        mifiter=_ndone + n + 1,
        coolingFn,
        rw_sd= rw_sd_matrix,
        param_rwIndex,
        tol=tol,
        maxFail,
        transform,
        _indices=_indices,
      )
      
    // } catch (error) {
    //   throw new Error(`Iterate the filtering stoped: ${error}`)
    // }
    paramMatrix = pfp.paramMatrix;
    // convRec[n+1,-c(1,2)] = pfp.coef;
    // convRec[n,c(1,2)] = [pfp.loglik, pfp.nfail];
    _indices = pfp.indices;
  }
  
  if (transform)
    pfp.paramMatrix = partrans(pomp,paramMatrix,dir="fromEstimationScale");
  console.log(pfp)
  
  return{
    // "mif2d.pomp",
    pfp,
    Nmif: Nmif,
    rw_sd: rw_sd,
    coolingType: coolingType,
    coolingFraction: coolingFraction,
    transform: transform,
    convRec: convRec
  }
}

mif2Pfilter = function (object, params, Np, mifiter, coolingFn, rw_sd, param_rwIndex,
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
    pmag =  rw_sd[nt].map(val => alpha * val);
    
    params = randwalk_perturbation(params, pmag, param_rwIndex);
    
    if (transform)
      tparams = partrans(object, params, dir="fromEstimationScale");
    
    if (nt === 0) {
      //get initial states
      let initparams = transform ? tparams : params;
      for (let i = 0; i < params.length; i++) {
        x.push(object.initializer(object, initparams[i]))
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
      console.error(`In mif2.js: error in calculation of weights: ${error}`)
    }
    
    let allFinite = weights.map(w => isFinite(w)).reduce((a, b) => a & b, 1)
    if (!allFinite) {
      throw new Error("In dmeasure: weights returns non-finite value")
    }

    // compute weighted mean at last timestep
    if (nt === ntimes - 1) {
      if (weights.map(w => w>0).reduce((a, b) => a || b, 0)) {
        object.coef = partrans(object, mathLib.mean(params, w = weights), dir="fromEstimationScale");
        
      } else {
        console.warn("filtering failure at last filter iteration, using unweighted mean for 'coef' ");
        object.coef = partrans(object, mathLib.mean(params), dir="fromEstimationScale");
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
        tol = tol,
        param_rwIndex
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
    params = xx.params;

    if (allFail) { // all particles are lost
      nfail = nfail + 1;
      if (nfail > maxFail)
        throw new Error("In mif2Pfilter: too many filtering failures")
    }    
  } // end of nt loop   
    
// console.log(loglik)
  if (nfail > 0) {
    console.log("warning! filtering failure occurred.");
  }

  return {
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


