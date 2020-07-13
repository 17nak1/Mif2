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

  if (nrepsx > nreps) {		// more ICs than parameters
    if (params.length === 1) params = new Array(nrepsx).fill(null).map(a => params[0]);
  } else if (nrepsx < nreps) {	// more parameters than ICs
    throw new Error("More parameters than ICs is not translated!")
  }
  // extract the process function. NOTE: only discrete-time translated
  let type = object.rprocess.type === "euler_sim" ? 3: 0;//TODO: type = *(INTEGER(GET_SLOT(rproc,install("type"))));
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
