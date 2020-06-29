const { euler_model_simulator } = require("./euler.js")
exports.rprocessInternal  = function (object, xstart, times, params, offset = 0, args) {
  let rv = do_rprocess(object, xstart, times, params, offset, args);
  return rv;
}

const do_rprocess = function (object, xstart, times, params, offset, args) {
  let ntimes = times.length;
  if (ntimes < 2) {
    throw new Error("in 'rprocess': length(times) < 2: with no transitions, there is no work to do.");
  }

  let off = Number(offset) ;
  if ((off < 0)||(off >= ntimes))
    throw new Error(`illegal 'offset' value ${off}`);
  let nvars = xstart[0].length;
  let nrepsx = xstart.length;  
  let npars = params[0].length;
  let nreps = params.length;

  // if (nrepsx > nreps) {		// more ICs than parameters
  //   if (nrepsx % nreps !== 0) {
  //     throw new Error("in 'rprocess': the larger number of replicates is not a multiple of smaller.");
  //   } else {
  //     double *src, *tgt;
  //     let dims = new Array(3);
  //     int j, k;
  //     dims[0] = npars;
  //     dims[1] = nrepsx;
  //     PROTECT(copy = duplicate(params)); nprotect++;
  //     PROTECT(params = makearray(2,dims)); nprotect++;
  //     setrownames(params,GET_ROWNAMES(GET_DIMNAMES(copy)),2);
  //     src = REAL(copy);
  //     tgt = REAL(params);
  //     for (j = 0; j < nrepsx; j++) {
  //       for (k = 0; k < npars; k++, tgt++) {
  //         *tgt = src[k+npars*(j%nreps)];
  //       }
  //     }
  //   }
  //   nreps = nrepsx;
  // } else if (nrepsx < nreps) {	// more parameters than ICs
  //   if (nreps % nrepsx != 0) {
  //     errorcall(R_NilValue,"in 'rprocess': the larger number of replicates is not a multiple of smaller.");
  //   } else {
  //     double *src, *tgt;
  //     int dims[2];
  //     int j, k;
  //     dims[0] = nvars; dims[1] = nreps;
  //     PROTECT(copy = duplicate(xstart)); nprotect++;
  //     PROTECT(xstart = makearray(2,dims)); nprotect++;
  //     setrownames(xstart,GET_ROWNAMES(GET_DIMNAMES(copy)),2);
  //     src = REAL(copy);
  //     tgt = REAL(xstart);
  //     for (j = 0; j < nreps; j++) {
  //       for (k = 0; k < nvars; k++, tgt++) {
  //         *tgt = src[k+nvars*(j%nrepsx)];
  //       }
  //     }
  //   }
  // }
  // extract the process function. NOTE: only discrete-time translated
  let type = object.rprocess.type === "euler_sim" ? 2: 0;
  let X;
  
  switch (type) {
    case 1: // one-step simulator
      fn = object.rprocess.stepFunction;
      deltat = 1.0;
      X = euler_model_simulator(fn, xstart, times, params, deltat, type, object)
      break;

    case 2: case 3: // discrete-time and Euler
      fn = object.rprocess.stepFunction;
      deltat = Number(object.rprocess.deltaT);        
      X = euler_model_simulator(fn, xstart, times, params, deltat, type, object)
      break;  

    case 4: // Gillespie's method
      throw new Error("in 'rprocess': Gillespie's method is not translated")
      
    case 0: default:
      throw new Error("'rprocess' is undefined. Note: only 'euler_sim' (discrete-time Euler) method is translated");
  }
      
  return X; 
}  







