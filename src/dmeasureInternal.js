exports.dmeasureInternal  = function (object, y, x, times, params, log) {
  let rv = do_dmeasure(object, y, x, times, params, log);
  return rv;
}

const do_dmeasure = function (object, y, x, times, params, give_log) {
  
  ntimes = Array.isArray(times) ? length(times) : 1;
  if (!times || ntimes < 1)
    throw new Error("In 'dmeasureInternal': times is not defined");
  
  let nvars = x[0].length;
  let nrepsx = x.length; 
  let npars = params[0].length;
  let nrepsp = params.length;

  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp != 0) || (nreps % nrepsx != 0))
    throw new Error("in 'dmeasureInternal': larger number of replicates is not a multiple of smaller"); 

  //  // construct state, parameter, covariate, observable indices
  //  oidx = INTEGER(PROTECT(name_index(Onames,pompfun,"obsnames","observables"))); nprotect++;
  //  sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames","state variables"))); nprotect++;
  //  pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;
  //  cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames","covariates"))); nprotect++;
   let ff = object.dmeasure;
   let F = new Array(nreps);

   for (k = 0; k < ntimes; k++) { // loop over times.Note:Only ntimes = 1 in ths translation
    for (let j = 0; j < nreps; j++) { // loop over replicates
      F[j] = ff(object, y, x[j], params[j],give_log);
    }
  }
  return F;
}