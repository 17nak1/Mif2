
snippet = {}
let mathLib = require('./mathLib')
//* Set the seed for rnorm-In R:RNGkind("L'Ecuyer-CMRG", normal.kind="Box-Muller");set.seed(1234) 
// const libR = require('lib-r-math.js')
// const {
//   Poisson,
//   rng: { MersenneTwister },
//   rng: { normal: { Inversion } }
// } = libR
// const mt = new MersenneTwister(0)// 
// const { rpois } = Poisson(new Inversion(mt))
// mt.init(1234)
snippet.rprocess = function (params, t, del_t, [S,E,I,R,H], pop, births) {
  var seas, beta, beta0, foi, R0, tt, va
  var trans = new Array(6).fill(0)
  var rate = new Array(6) 
  var deltaT = 14 / 365.25
  var dt = 1 / 365.25 
  
  R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4] 
  beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma
  
  va = snippet.rprocessVaccine(t)
  tt = (t - Math.floor(t)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - amplitude
  }                 
  beta = R0 * (gamma + mu) * (sigma + mu) * seas / sigma //seasonal transmission rate
  foi = beta * I / pop
  rate[0] = foi//force of infection
  rate[1] = mu// natural S death
  rate[2] = sigma// rate of ending of latent stage
  rate[3] = mu// natural E death
  rate[4] = gamma// recovery
  rate[5] = mu// natural I death 
   
  // births = mathLib.rpois(birthrate * (1 - va) * del_t )// Poisson births
  mathLib.reulermultinom(2, Math.round(S), 0, del_t, 0, rate, trans)
  mathLib.reulermultinom(2, Math.round(E), 2, del_t, 2, rate, trans)
  mathLib.reulermultinom(2, Math.round(I), 4, del_t, 4, rate, trans)
  S += (births - trans[0] - trans[1])
  E += (trans[0] - trans[2] - trans[3]) 
  I += (trans[2] - trans[4] - trans[5]) 
  R = pop - S - E - I
  H += trans[4] 
  return [S, E, I, R, H]
}
snippet.rprocessVaccine = function(t) {
  var vaccineRate
  if (t < 1968)
    vaccineRate = 0
  else if (t >= 1968 && t <= 1969)
    vaccineRate = 0.33
  else if (t >= 1969 && t <= 1970)
    vaccineRate = 0.46
  else if (t >= 1970 && t <= 1971)
    vaccineRate = 0.51
  else if (t >= 1971 && t <= 1972)
    vaccineRate = 0.53
  else if (t >= 1972 && t <= 1973)
    vaccineRate = 0.52
  else if (t >= 1973 && t <= 1974)
    vaccineRate = 0.46
  else if (t >= 1974 && t <= 1975)
    vaccineRate = 0.46
  else if (t >= 1975 && t <= 1976)
    vaccineRate = 0.48
  else if (t >= 1976 && t <= 1977)
    vaccineRate = 0.48
  else if (t >= 1977 && t <= 1978)
    vaccineRate = 0.51
  else if (t >= 1978 && t <= 1979)
    vaccineRate = 0.53;
  else if (t >= 1979 && t <= 1980)
    vaccineRate = 0.55;
  else if (t >= 1980 && t <= 1981)
    vaccineRate = 0.58;
  else if (t >= 1981 && t <= 1982)
    vaccineRate = 0.60
  else if (t >= 1982 && t <= 1983)
    vaccineRate = 0.63
  else if (t >= 1983 && t <= 1984)
    vaccineRate = 0.68
  else if (t >= 1984 && t <= 1985)
    vaccineRate = 0.71
  else if (t >= 1985 && t <= 1988)
    vaccineRate = 0.76
  else if (t >= 1988 && t <= 1989)
    vaccineRate = 0.814
  else if (t >= 1989 && t <= 1990)
    vaccineRate = 0.9488
  else if (t >= 1990 && t <= 1991)
    vaccineRate = 0.9818
  else if (t >= 1991 && t <= 1992)
    vaccineRate = 0.90
  else if (t >= 1992 && t <= 1993)
    vaccineRate = 0.92
  else if (t >= 1993 && t <= 1994)
    vaccineRate = 0.91
  else if (t >= 1994 && t <= 1995)
    vaccineRate = 0.91
  else if (t >= 1995 && t <= 1996)
    vaccineRate = 0.92
  else if (t >= 1996 && t <= 1997)
    vaccineRate = 0.92
  else if (t >= 1997 && t <= 1998)
    vaccineRate = 0.91
  else if (t >= 1998 && t <= 1999)
    vaccineRate = 0.88
  else if (t >= 1999 && t <= 2000)
    vaccineRate = 0.88
  else if (t >= 2000 && t <= 2001)
    vaccineRate = 0.87
  else if (t >= 2001 && t <= 2002)
    vaccineRate = 0.84
  else if (t >= 2002 && t <= 2003)
    vaccineRate = 0.82
  else if (t >= 2003 && t <= 2004)
    vaccineRate = 0.80
  else if (t >= 2004 && t <= 2005)
    vaccineRate = 0.81
  else if (t >= 2005 && t <= 2006)
    vaccineRate = 0.84
  else if (t >= 2006 && t <= 2007)
    vaccineRate = 0.85
  else if (t >= 2007 && t <= 2008)
    vaccineRate = 0.85
  else if (t >= 2008 && t <= 2009)
    vaccineRate = 0.85
  else if (t >= 2009 && t <= 2010)
    vaccineRate = 0.88
  else
    vaccineRate = 0.89
  return vaccineRate
}

snippet.initz = function(pop, S_0, E_0, I_0, R_0, H) {
  var m = pop / (S_0 + E_0 + R_0 + I_0),
    S = Math.round(m * S_0),
    E = Math.round(m * E_0),
    I = Math.round(m * I_0),
    R = Math.round(m * R_0),
    H = 0
  return [S, E, I, R, H]
}

snippet.dmeasure = function (rho, psi, H, dCases, giveLog = 1) {
  var lik
  var mn = rho * H
  var v = mn * (1.0 - rho + psi * psi * mn)
  var tol = 1.0e-18
  var modelCases = Number(dCases)
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
    }
    if (giveLog) lik = Math.log(lik)
  } else {
    lik = (giveLog) ? 0 : 1;
  }
  return lik
}

snippet.rmeasure = function (H, rho, psi) {
  var mn = rho * H
  var v = mn * (1.0 - rho + psi * psi * mn)
  var tol = 1.0e-18
  var cases = mathLib.rnorm(mn, Math.sqrt(v) + tol)
  if (cases > 0) {
    cases = Math.round(cases)
  } else {
    cases = 0
  }
  return cases
}

module.exports = snippet
