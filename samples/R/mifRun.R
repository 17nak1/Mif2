setwd(getwd())
# library(pomp)
start_time <- Sys.time()


######################################################  Model Snippet
rproc <- Csnippet("
                  double seas, beta, foi;
                  double births, va, tt;
                  double rate[6], trans[6];


                  va = 0;


                  // term-time seasonality
                  tt = (t-floor(t))*365.25;
                  if ((tt>=7&&tt<=100) || (tt>=115&&tt<=199) || (tt>=252&&tt<=300) || (tt>=308&&tt<=356))
                  seas = 1.0+amplitude*0.2411/0.7589;
                  else
                  seas = 1.0-amplitude;

                  // transmission rate
                  beta = R0*(gamma+mu)*(sigma+mu)*seas/sigma;  //seasonal transmission rate
                  // expected force of infection
                  foi = beta*I/pop;

                  rate[0] = foi;  //         force of infection
                  rate[1] = mu;             // natural S death
                  rate[2] = sigma;        // rate of ending of latent stage
                  rate[3] = mu;             // natural E death
                  rate[4] = gamma;        // recovery
                  rate[5] = mu;             // natural I death

                  //if( t>= 1934.9 && t< 1944.1)
                  //printf(\"%f  %f\\n\", t, pop);
                  // Poisson births
                  births = rpois(birthrate*(1-va)*dt);
                  //printf(\" %f, %f \\n\", t, births);
                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  S += births - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  R = pop - S - E - I;
                  H += trans[4];           // true incidence
                  "
                  )

# The above uses true incidence for reporting

initz <- Csnippet("
                  double m = pop/(S_0+E_0+I_0+R_0);
                  S = nearbyint(m*S_0);
                  E = nearbyint(m*E_0);
                  I = nearbyint(m*I_0);
                  R = nearbyint(m*R_0);
                  H = 0;
                  ")
# Sampling from the normal approximation of the binomial distribution
dmeas <- Csnippet("
                  double m = rho*H;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (R_FINITE(cases)) {
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  if (give_log) lik = log(lik);
                  } else {
                  lik = (give_log) ? 0 : 1;
                  }
                  ")
rmeas <- Csnippet("

                  double m = rho*H;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")




toEst <- Csnippet("
                  Tmu = log(mu);
                  Tpsi = log(psi);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  TR0 = log(R0);
                  Trho = logit(rho);
                  Tamplitude = logit(amplitude);
                  to_log_barycentric (&TS_0, &S_0, 4);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tpsi = exp(psi);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    TR0 = exp(R0);
                    Trho = expit(rho);
                    Tamplitude = expit(amplitude);
                    from_log_barycentric (&TS_0, &S_0, 4);
                    ")

params_mod <- c("R0","amplitude","gamma","mu","sigma","rho","psi")

params_ic <- c("S_0","E_0","R_0","I_0")

statenames <- c("S","E","I","R","H")

zeronames <- c("H")
#########################################            data
cureentfolder <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
London_BiData <- read.csv(file.path(cureentfolder, "London_BiDataMainsh.csv"))
London_covar <- read.csv(file.path(cureentfolder, "London_covar.csv"))
rw_size <- 0.05
rw <- rw.sd (R0=ifelse(time<1944,0,rw_size),
             amplitude=ifelse(time<1944,0,rw_size),
             mu=ifelse(time<1944,0,rw_size),
             rho=ifelse(time<1944,0,rw_size),
             psi=ifelse(time<1944,0,rw_size),
             S_0=ifelse(time<1944,0,rw_size),
             E_0=ifelse(time<1944,0,rw_size),
             I_0=ifelse(time<1944,0,rw_size),
             R_0=ifelse(time<1944,0,rw_size))
#########################################       make pomp
pomp(
  data = London_BiData,
  times="time",
  t0=1940,#with(get(paste0(name,"_BiData")),2*time[1]-time[2]),
  rprocess = euler.sim(rproc,delta.t=1/365.25),
  rmeasure=rmeas,
  covar=London_covar,
  tcovar="time",
  dmeasure=dmeas,
  zeronames=zeronames,
  initializer=initz,
  toEstimationScale=toEst,
  fromEstimationScale=fromEst,
  statenames=statenames,
  paramnames=c(params_mod,params_ic)
) -> m1


current_params= c(R0=3.132490e+01 , amplitude=3.883620e-01 , gamma=7.305000e+01 , mu=6.469830e-04 , sigma=4.566000e+01 ,rho= 4.598709e-01 ,psi= 1.462546e-01 ,S_0= 3.399189e-02 ,E_0=2.336327e-04 ,R_0=9.657741e-01,I_0=4.221789e-07)
params_ic_fit <- character(0)
params_mod_fit <- c("R0", "amplitude", "mu" ,"rho" ,"psi")
no_mif <- 1
no_particles <- 100
cool_fraction <- 0.05
tstart = 1
tend = 40
# datasetj <- as.data.frame(read.csv(file.path(cureentfolder, "predmean.csv")))
start_time <- Sys.time()


  coef(m1) <- c(current_params)
  mif2(m1,
       Nmif=no_mif,
       start=current_params,
       transform=T,
       ivps=params_ic_fit,
       pars=params_mod_fit,
       rw.sd=rw,
       Np=no_particles,
       var.factor=2,
       cooling.type="hyperbolic",
       cooling.fraction=cool_fraction
  ) -> sets_mf


  current_params <- coef(sets_mf)
  Sys.time()- start_time

#
# datapredict <- as.data.frame(ss@pred.mean)
# pm=c();for(i in tstart:tend){pm[i]=datapredict[1,i]}
# pfilter(m1,params=current_params,Np=np,filter.mean = T,pred.mean=T, max.fail=3000) -> ss; ss@loglik
# datapredict2 <- as.data.frame(ss@pred.mean)
# pm2=c();for(i in tstart:tend){pm2[i]=datapredict2[1,i]}
# plot(datasetj$S[c(tstart :tend)], type ="l",main= np,col = "red", ylab = "JS")
# points(pm[c(tstart :tend)], type ="l")
# plot(pm2[c(tstart :tend)], type ="l",main= np,col="blue", ylab = "R")
# points(pm[c(tstart :tend)], type ="l")

