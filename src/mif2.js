const snippet = require("./modelSnippet");
const { randwalk_perturbation } = require("./mif2cRW.js");
const { rprocessInternal } = require("./rprocessInternal.js")

const cooling = function(type, fraction, ntimes) {
  switch(type){
  case "geometric":
    let factor = Math.pow(fraction, (1/50));
    return function (nt, m) {
      let alpha = Math.pow(factor, (nt / ntimes + m - 1));
      return {alpha: alpha, gamma: alpha ** 2};
    }
    
  case "hyperbolic":
    if (fraction < 1) {
      let scal = (50 * ntimes * fraction - 1) / (1 - fraction);
      return function (nt, m) {
        let alpha = (1 + scal) / (scal + nt + ntimes * (m - 1));
        return {alpha: alpha, gamma: alpha ** 2};
      }
    } else {
      return function (nt, m) {
        return {alpha: 1, gamma: 1};
      }
    }
  }
}

/**
 * Rescale the parameters.
 * @param {object} pomp 
 * @param {array} params 
 * @param {string} dir 
 */
const partrans = function (pomp, params, dir = ["fromEstimationScale","toEstimationScale"]) {
  if (!Array.isArray(params[0])) params = new Array(params);
  
  switch(dir){
  case "fromEstimationScale":
    for (let i = 0 ; i < params.length; i++) {
      params[i] = pomp.fromEstimationScale(params[i]);
    }
    break;
    
  case "toEstimationScale":
    for (let i = 0 ; i < params.length; i++) {
      params[i] = pomp.toEstimationScale(params[i]);
    }
    break;
  } 
  
  if (params.length === 1) return [...params];
  return params;
}

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

  convRec = new Array(Nmif+1).fill(Array(start.length + 2));
  convRec[1] = [null, null, ...start];//[loglik,nfail,...start]
  if (transform)
    paramMatrix = partrans(pomp,paramMatrix,dir="toEstimationScale");
  
  // Iterate the filtering main loop 
  let pfp;
  for (let n =0; n <= Nmif; n++) {
    try {
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
        verbose=verbose,
        transform,
        _indices=_indices,
      )
    } catch (error) {
      throw new Error(`Iterate the filtering stoped: ${error}`)
    }
    // paramMatrix = pfp@paramMatrix
    // conv.rec[n+1,-c(1,2)] <- coef(pfp)
    // conv.rec[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)
    // .indices <- pfp@indices

    // if (verbose) cat("mif2 it
  }
  let aa
  if (transform)
    aa = partrans(pomp,paramMatrix,dir="fromEstimationScale");
  // console.log(aa)
}

mif2Pfilter = function (object, params, Np, mifiter, coolingFn, rw_sd, param_rwIndex,
  tol = 1e-17, maxFail = Inf, verbose, transform, _indices = 0)
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
  for (let nt = 0; nt < 1; nt++) {
    alpha = coolingFn(nt + 1,mifiter).alpha;
    pmag =  rw_sd[nt].map(val => alpha * val);
    params = randwalk_perturbation(params, pmag, param_rwIndex);
    
    if (transform)
      tparams = partrans(object, params, dir="fromEstimationScale");
    
    if (nt == 0) {
      //get initial states
      params = transform ? tparams : params;
      let x = [];
      for (let i = 0; i < params.length; i++) {
        x.push(object.initializer(object, params[i]))
      }
      let X = [], xarray;
      for (let i = 0; i < params.length; i++) {
        params[i] = transform ? tparams[i] : params[i];
        xarray = rprocessInternal(
          object,
          xstart=x[i],
          times=[times[nt],times[nt+1]],
          params[i],
          offset=1
        )
        X.push(xarray)
      }
      
      
      // console.log(x)
    } 
  }

}


