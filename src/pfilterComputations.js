// examines weights for filtering failure.
// computes conditional log likelihood and effective sample size.
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// weights are used in filtering mean computation.
// if length(weights) == 1, an unweighted average is computed.
// tracks ancestry of particles if desired.
// returns all of the above in a named list.

// predmean: calculate prediction means?
// predvar: calculate prediction variances?
// filtmean: calculate filtering means?
// trackancestry: track ancestry?
// onepar: are all cols of 'params' the same?
const mathLib = require("./mathLib.js")
exports.pfilter_computations = function (x, params, Np, rw_sd, predmean, predvar,
  filtmean, trackancestry, onepar, weights, toler, param_rwIndex)
{
  let nvars = x[0].length;
  let nreps = x.length;
  let npars = params[0].length;
  let nparreps = params.length;

  if (nreps % params.length != 0)
    throw new Error("In pfilter_computations: ncol('states') should be a multiple of ncol('params')"); // # nocov
  let nrw = Array.isArray(rw_sd) ? rw_sd.length : 0;	     // number of parameters that are variable
  let do_rw = nrw > 0;                                     // do random walk in parameters?
  let rw_names = null;	                                   
  if (do_rw) {
    // names of parameters undergoing random walk
    rw_names = rw_sd_name;
  }
  let do_par_resamp = !onepar || do_rw; // should we do parameter resampling?

  if (do_par_resamp) {
    if (nparreps !== nreps)
      throw new Error("In pfilter_computations: ncol('states') should be equal to ncol('params')"); // # nocov
  }
  let ess ;     // effective sample size
  let loglik ;  // log likelihood
  let fail ;  	// particle failure?

  // check the weights and compute sum and sum of squares
  let nlost = 0;
  let w = 0;
  let ws = 0;
  let all_fail = 0;
  for (let k = 0; k < nreps; k++) {
    if (weights[k] > toler) {
      w += weights[k];
      ws += weights[k] * weights[k];
    } else {			// this particle is lost
      weights[k] = 0;
      nlost++;
    }
  }
  if (nlost >= nreps) all_fail = 1; // all particles are lost
  if (all_fail) {
    loglik = Math.log(toler); // minimum log-likelihood
    ess = 0;		              // zero effective sample size
  } else {
    loglik = Math.log(w / nreps); // mean of weights is likelihood
    ess = w * w / ws;	            // effective sample size
  }
  
  fail = all_fail;
  let pidx;                       // indices of parameters undergoing random walk
  let xp, lv;
  
  if (do_rw) {
    pidx = param_rwIndex;
    xp = params;
    lv = nvars + nrw;
  } else {
    pidx = null;
    lv = nvars;
  }

  let xpm;
  if (predmean || predvar) {
    xpm = new Array(lv);
  }

  let xpv;
  if (predvar) {
    xpv = new Array(lv);
  }

  let xfm;
  if (filtmean) {
    if (do_rw) {
      xfm = new Array(nvars + npars);
    } else {
      xfm = new Array(nvars);
    }
  }

  if (trackancestry) {
    PROTECT(anc = NEW_INTEGER(np)); nprotect++;
    xanc = INTEGER(anc);
  }
  let sum = 0;
  for (let j = 0; j < nvars; j++) {	// state variables

    // compute prediction mean
    if (predmean || predvar) {
      sum = 0;
      for (let k = 0; k < nreps; k++) {
        sum += x[k][j];
      }
      sum /= nreps;
      xpm[j] = sum;
    }

    // compute prediction variance
    if (predvar) {
      let sumsq = 0;
      let vsq;
      for (let k = 0; k < nreps; k++) {
        vsq = x[k][j] - sum;
        sumsq += vsq * vsq;
      }
      xpv[j] = sumsq / (nreps - 1);
    }

    //  compute filter mean
    if (filtmean) {
      if (all_fail) {		// unweighted average
        let ws = 0;
        for (let k = 0; k < nreps; k++) {
          ws += x[k][j];
        }
        xfm[j] = ws/ nreps;
      } else { 			// weighted average
        let ws = 0;
        for (let k = 0; k < nreps; k++) {
          ws += x[k][j] * weights[k];
        }
        xfm[j] = ws / w;
      }
    }
  }

  // compute means and variances for parameters (if needed)
  let offset;
  
  if (do_rw) {
    for (let j = 0; j < nrw; j++) {
      offset = pidx[j];		// position of the parameter

      if (predmean || predvar) {
        let sum = 0;
        for (let k = 0; k < nreps; k++) {
          sum += xp[k][offset];
        }
        sum /= nreps;
        xpm[j] = sum;
      }

      if (predvar) {
        let sumsq = 0;
        let vsq;
        for (let k = 0; k < nreps; k++) {
          vsq = xp[k][offset] - sum;
          sumsq += vsq * vsq;
        }
        xpv[j] = sumsq / (nreps - 1);
      }
    }

    if (filtmean) {
      let ws;
      for (let j = 0; j < npars; j++) {
        if (all_fail) {		// unweighted average
          ws = 0;
          for (let k = 0; k < nreps; k++) {
            ws += xp[k][j];
          }
          xfm[j] = ws/ nreps;
        } else {		// weighted average
          ws = 0;
          for (let k = 0; k < nreps; k++) {
            ws += xp[k][j] * xw[k];
          }
          xfm[j] = ws / w;
        }
      }
    }
  }
  let newstates, newparams;
  if (!all_fail) { // resample the particles unless we have filtering failure
    let  xdim = new Array(2);
    let  sample =  new Array(Np);
    
    // create storage for new states
    xdim[0] = nvars; xdim[1] = Np;
    newstates = new Array(Np).fill(null).map(()=>new Array(nvars));

    // create storage for new parameters
    if (do_par_resamp) {
      newparams = new Array(Np).fill(null).map(()=>new Array(npars));
    }

    // resample
    mathLib.nosortResamp(nreps, weights, Np, sample, 0);
    
    for (let k = 0; k < Np; k++) { // copy the particles
      for (let j = 0; j < nvars; j++) {
        newstates[k][j] = x[sample[k]][j];
      }
      
      if (do_par_resamp) {
        for (let j = 0; j < npars; j++){
          newparams[k][j] = params[sample[k]][j];
        } 
      }
      if (trackancestry) xanc[k] = sample[k] + 1;
    }
    
  } else { // don't resample: just drop 3rd dimension in x prior to return
    if (trackancestry)
      for (let k = 0; k < np; k++) {
        xanc[k] = k + 1;
      }
  }
  
  let retval = {"fail": fail, "loglik": loglik, "ess": ess,
   "states": null, "params": null, "pm": null, "pv": null, "fm": null, "ancestry": null};

  if (all_fail) {
    retval.states = x;
  } else {
    retval.states = newstates;
  }
  
  if (all_fail || !do_par_resamp) {
    retval.params = params;
  } else {
    retval.params = newparams;
  }

  if (predmean) {
    retval.pm = xpm;
  }
  if (predvar) {
    retval.pv = xpv;
  }
  if (filtmean) {
    retval.fm = xfm;
  }
  if (trackancestry) {
    retval.ancestry = xanc;
  }

  return retval;
}