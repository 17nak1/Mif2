exports.dmeasureInternal  = function (object, y, x, times, params, log) {
  let rv = do_dmeasure(object, y, x, times, params, log);
  return rv;
}

const do_dmeasure = function (object, y, x, times, params, give_log) {
  ntimes = Array.isArray(times) ? length(times) : 1;
  if (!times || ntimes < 1)
    throw new Error("In 'dmeasureInternal': times is not defined");
  
  let nvars = Object.keys(x[0]).length;
  let nrepsx = x.length; 
  let npars = Object.keys(params[0]).length;
  let nrepsp = params.length;
  
  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp != 0) || (nreps % nrepsx != 0))
    throw new Error("in 'dmeasureInternal': larger number of replicates is not a multiple of smaller"); 

  let ff = object.dmeasure;
   let F = new Array(nreps);

   for (k = 0; k < ntimes; k++) { // loop over times.Note:Only ntimes = 1 in ths translation
    for (let j = 0; j < nreps; j++) { // loop over replicates
      if (nrepsp === 1) params[j] = params[0];
      F[j] = ff(y, x[j], params[j],give_log);
    }
  }
  return F;
}