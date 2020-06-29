const mathLib = require("./mathLib");
const { Z_BINARY } = require("zlib");


exports.euler_model_simulator  = function(func, xstart, times, params, deltat,
   method, object)
{
  let zeronames = object.zeronames;
  let statenames = object.statenames; 
  if (deltat <= 0)
    throw new Error("In euler.js: 'delta.t' should be a positive number");
  let nvars = xstart[0].length;
  let nreps = xstart.length;
  let npars = params.length;
  let ntimes = times.length;

  let zidx = mathLib.index(statenames, zeronames);

  let t = times[0];
  let xt = xstart;

  for (let step = 1; step < ntimes; step++) {
    if (t > times[step]) {
      throw new Error("In euler.js: 'times' is not an increasing sequence");
    }
    
    // set accumulator variables to zero
    for (j = 0; j < nreps; j++) {
      for (i = 0; i < nvars; i++) {
        if(zidx[i] !== null)
         xt[j][i] = 0;
      } 
    }
    switch (method) {
      case 1:			// one step
        dt = times[step] - t;
        nstep = (dt > 0) ? 1 : 0;
        break;
      case 2:			// fixed step
        dt = deltat;
        nstep = numMapSteps(t, times[step], dt);
        break;
      case 3:			// Euler method
        dt = deltat;
        nstep = numEulerSteps(t, times[step], dt);
        break;
      default:
        throw new Error("In euler.js: unrecognized 'method'"); // # nocov
    }

    for (let k = 0; k < nstep; k++) { // loop over Euler steps
      //interpolate
      for (let j = 0 ; j < nreps; j++) { // loop over replicates
         xt[j] = func(object, xt[j], params[j], t, dt);
      }
      t += dt;

      if ((method == 3) && (k == nstep-2)) { // penultimate step
        dt = times[step] - t;
        t = times[step] - dt;
      }  
    }
  }
  
  return xt;
}

// helpers
const numMapSteps = function (t1, t2, dt) {
  let DOUBLE_EPS = 10e-8
  let tol = Math.sqrt(DOUBLE_EPS)
  let nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}

const numEulerSteps = function (t1, t2, dt) {
  let DOUBLE_EPS = 10e-8
  let tol = Math.sqrt(DOUBLE_EPS)
  let nstep;
  // nstep will be the number of Euler steps to take in going from t1 to t2.
  // note also that the stepsize changes.
  // this choice is meant to be conservative
  // (i.e., so that the actual dt does not exceed the specified dt
  // by more than the relative tolerance 'tol')
  // and to counteract roundoff error.
  // It seems to work well, but is not guaranteed:
  // suggestions would be appreciated.

  if (t1 >= t2) {
    dt = 0;
    nstep = 0;
  } else if (t1 + dt >= t2) {
    dt = t2 - t1;
    nstep = 1;
  } else {
    nstep = Math.ceil((t2 - t1) / dt / (1 + tol));
    dt = (t2 - t1)/ nstep;
  }
  return nstep;
}