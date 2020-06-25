/**
 *  @file        modeSnippet.js
 *               Makes the model and its dependencies.
 *                 
 *  @autor       Nazila Akhavan, nazila@kingsds.network
 *  @date        Feb 2019
 */

snippet = {}
let mathLib = require('./mathLib')
let rpois = require('./rpois')

snippet.rprocess = function (pomp, params, t, del_t, [S,E,I,R,H]) {
  let pop =pomp.population(t);
  let birthrate = pomp.birthrate(t)
  let seas, beta, beta0, foi, R0, tt, va
  let trans = new Array(6).fill(0)
  let rate = new Array(6) 
  let deltaT = 14 / 365.25
  let dt = 1 / 365.25 
  
  R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4] 
  beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma
  
  va = 0;
  tt = (t - Math.floor(t)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - amplitude
  }                 
  beta = R0 * (gamma + mu) * (sigma + mu) * seas / sigma  //seasonal transmission rate
  foi = beta * I / pop
  rate[0] = foi            //force of infection
  rate[1] = mu             // natural S death
  rate[2] = sigma          // rate of ending of latent stage
  rate[3] = mu             // natural E death
  rate[4] = gamma          // recovery
  rate[5] = mu             // natural I death 
   
  let births = rpois.rpoisOne(birthrate * (1 - va) * del_t )// Poisson births
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

snippet.initz = function(pomp, params) {
  let ind = pomp.pIndex;
  let m = pomp.population(pomp.t0) / (params[ind["S"]] + params[ind["E"]] + params[ind["R"]] + params[ind["I"]]);
  let S_0 = Math.round(m * params[ind["S"]]);
  let E_0 = Math.round(m * params[ind["E"]]);
  let I_0 = Math.round(m * params[ind["I"]]);
  let R_0 = Math.round(m * params[ind["R"]]);
  let H_0 = 0;
  return [S_0, E_0, I_0, R_0, H_0];
}

snippet.dmeasure = function (rho, psi, H, dCases, giveLog = 1) {
  let lik
  let mn = rho * H
  let v = mn * (1.0 - rho + psi * psi * mn)
  let tol = 1.0e-18
  let modelCases = Number(dCases)
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
  let mn = rho * H
  let v = mn * (1.0 - rho + psi * psi * mn)
  let tol = 1.0e-18
  let cases = mathLib.rnorm(mn, Math.sqrt(v) + tol)
  if (cases > 0) {
    cases = Math.round(cases)
  } else {
    cases = 0
  }
  return cases
}
/**
 * paramnames is an array of all parameters in the model and all code is based on the order
 * of params in this array. The first index is zero i.e R0 = params[0]
 */
snippet.paramnames = ["R0","amplitude","gamma","mu","sigma","rho","psi", "S_0", "E_0", "I_0", "R_0"];
snippet.zeronames = ["H"];
snippet.statenames = ["S","E","I","R","H"];

snippet.toEst = function(params) {
  let mu = Math.log(params[3]);
  let psi = Math.log(params[6]);
  let sigma = Math.log(params[4]);
  let gamma = Math.log(params[2]);
  let R0 = Math.log(params[0]);
  let rho = mathLib.logit(params[5]);
  let amplitude = mathLib.logit(params[1]);
  let states = mathLib.toLogBarycentric([params[7], params[8], params[9], params[10]]);
  //Parameters order should be the same as paramnames.
  return [R0, amplitude, gamma, mu, sigma, rho, psi, ...states];
}

snippet.fromEst = function(params) {
  let mu = Math.exp(params[3]);
  let psi = Math.exp(params[6]);
  let sigma = Math.exp(params[4]);
  let gamma = Math.exp(params[2]);
  let R0 = Math.exp(params[0]);
  let rho = mathLib.expit(params[5]);
  let amplitude = mathLib.expit(params[1]);
  let states = mathLib.fromLogBarycentric([params[7], params[8], params[9], params[10]]);
  //Parameters order should be the same as paramnames.
  return [R0, amplitude, gamma, mu, sigma, rho, psi, ...states];
}

module.exports = snippet
