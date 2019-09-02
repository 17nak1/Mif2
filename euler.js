





euler_model_simulator = function (func, xstart, t0, times, params, deltat, rprocmode method,  accumIndexArr,  covar,  args) {
  var nvars, npars, nreps, ntimes, nzeros, ncovars
  var cvec, X, fn

  if (deltat <= 0) throw "'delta.t' should be a positive number."
  var nvars = xstart[0].length, nreps = xstart.length
  var npars = dim[0]
  var ntimes = times.length

  
  // t0 = tstart
  if (t0 > times[0]) throw "'t0' must be no later than 'times[1]'."


  // indices of accumulator variables
  var nzeros = accumIndexArr.length
  // var *zidx = INTEGER(PROTECT(matchnames(Snames,accumvars,"state variables"))); nprotect++;

  // extract user function
  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames)); nprotect++;

  // array to hold results
  PROTECT(X = ret_array(nvars,nreps,ntimes,Snames)); nprotect++;


  var *pidx = 0, *sidx = 0, *cidx = 0;
  pomp_onestep_sim *ff = NULL;

  switch (mode) {

  case Rfun: {

    // construct list of all arguments
    PROTECT(args = add_args(args,Snames,Pnames,Cnames)); nprotect++;

  }

    break;

  case native: case regNative: {

    // construct state, parameter, covariate indices
    sidx = INTEGER(GET_SLOT(func,install("stateindex")));
    pidx = INTEGER(GET_SLOT(func,install("paramindex")));
    cidx = INTEGER(GET_SLOT(func,install("covarindex")));

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    set_pomp_userdata(args);
    GetRNGstate();

  }

    break;

  default:

    errorcall(R_NilValue,"unrecognized 'mode' %d",mode); // # nocov

  }

  // main computation loop
  
   
  var xt = X
  var time = times
  var t = t0
  var dt, nstep
  for (step = 0 ; step < ntimes; step++, xt += nvars*nreps) {
    var posn = NULL;
    var nstep = 0;
    
    if (t > time[step]) throw "'times' must be an increasing sequence."

    // set accumulator variables to zero
    for (j = 0; j < nreps; j++) {
      for (i = 0; i < nzeros; i++) {
        xt[accumIndexArr[i]][j] = 0
      }
    }

    // determine size and number of time-steps
    // switch (method) {
    // case onestep: default:  // one step
      dt = time[step]-t;
      nstep = (dt > 0) ? 1 : 0;
    //   break;
    // case discrete:      // fixed step
    //   dt = deltat;
    //   nstep = num_map_steps(t,time[step],dt);
    //   break;
    // case euler:     // Euler method
    //   dt = deltat;
    //   nstep = num_euler_steps(t,time[step],&dt);
    //   break;
    // }

    // loop over individual steps
    for (k = 0; k < nstep; k++) {

      // interpolate the covar functions for the covariates
      table_lookup(&covariate_table,t,cov);

      // loop over replicates
       // *ap, *pm, *xm, *ps = REAL(params);

      for (j = 0, pm = ps, xm = xt; j < nreps; j++, pm += npars, xm += nvars) {

        // switch (mode) {

        // case Rfun: {

           ans, nm;

          if (j == 0 && k == 0) {

            ans = eval_call(fn,args,&t,&dt,xm,nvars,pm,npars,cov,ncovars)

            // PROTECT(nm = GET_NAMES(ans)); nprotect++;
            // if (invalid_names(nm))
            //   errorcall(R_NilValue,"'rprocess' must return a named numeric vector.");
            // posn = INTEGER(PROTECT(matchnames(Snames,nm,"state variables"))); nprotect++;

            // ap = REAL(AS_NUMERIC(ans));

            for (i = 0; i < nvars; i++) xm[posn[i]] = ans[i];

          } else {

            PROTECT(ans = eval_call(fn,args,&t,&dt,xm,nvars,pm,npars,cov,ncovars));
            ap = REAL(AS_NUMERIC(ans));
            for (i = 0; i < nvars; i++) xm[posn[i]] = ap[i];
            UNPROTECT(1);

          }

        // }

        //   break;

        // case native: case regNative: {

        //   (*ff)(xm,pm,sidx,pidx,cidx,cov,t,dt);

        // }

        //   break;

        // default:

        //   errorcall(R_NilValue,"unrecognized 'mode' %d",mode); // # nocov

        // }

      }

      t += dt;

      if ((method == euler) && (k == nstep-2)) { // penultimate step
        dt = time[step]-t;
        t = time[step]-dt;
      }

    }

    if (step < ntimes-1)
      memcpy(xt+nvars*nreps,xt,nreps*nvars*sizeof());

  }

  // clean up
  switch (mode) {

  case native: case regNative: {
    PutRNGstate();
    unset_pomp_userdata();

  }

    break;

  case Rfun: default:

    break;

  }

  UNPROTECT(nprotect);
  return X;

}

var num_euler_steps ( t1,  t2,  *dt) {
   tol = sqrt(_EPS);
  var nstep;
  // nstep will be the number of Euler steps to take in going from t1 to t2.
  // note also that the stepsize changes.
  // this choice is meant to be conservative
  // (i.e., so that the actual dt does not exceed the specified dt
  // by more than the relative tolerance 'tol')
  // and to counteract roundoff error.
  // It seems to work well, but is not guaranteed:
  // suggestions would be appreciated.

  if (t1 >= t2) {
    *dt = 0.0;
    nstep = 0;
  } else if (t1+*dt >= t2) {
    *dt = t2-t1;
    nstep = 1;
  } else {
    nstep = (int) ceil((t2-t1)/(*dt)/(1+tol));
    *dt = (t2-t1)/(() nstep);
  }
  return nstep;
}

var num_map_steps ( t1,  t2,  dt) {
   tol = sqrt(_EPS);
  var nstep;
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = (int) floor((t2-t1)/dt/(1-tol));
  return (nstep > 0) ? nstep : 0;
}
© 2019 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
