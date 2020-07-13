/**
 *  @file        modeSnippet.js
 *               Makes the model and its dependencies.
 *                 
 *  @autor       Nazila Akhavan, nazila@kingsds.network
 *  @date        Feb 2019
 */

snippet = {}
let mathLib = require('./mathLib');
let rpois = require('./rpois');

snippet.rprocess = function (states, params, t, dt, covar) {
  
  let S = states.S;
  let E = states.E;
  let I = states.I;
  let R = states.R;
  let H = states.H;

  let R0 = params.R0;
  let amplitude = params.amplitude;
  let gamma = params.gamma;
  let mu = params.mu;
  let sigma = params.sigma ;

  let pop = covar.pop;
  let birthrate = covar.birthrate;
  let seas, beta, beta0, foi, tt, va;
  let length = 3;
  let trans = new Array(length * 2).fill(0);
  let rate = new Array(length * 2); 

  
  
  beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma;
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
   
  let births = rpois.rpoisOne(birthrate * (1 - va) * dt )// Poisson births
  mathLib.reulermultinom(2, Math.round(S), 0, dt, 0, rate, trans)
  mathLib.reulermultinom(2, Math.round(E), 2, dt, 2, rate, trans)
  mathLib.reulermultinom(2, Math.round(I), 4, dt, 4, rate, trans)
  S += (births - trans[0] - trans[1])
  E += (trans[0] - trans[2] - trans[3]) 
  I += (trans[2] - trans[4] - trans[5]) 
  R = pop - S - E - I
  H += trans[4] 
  return {S, E, I, R, H}
}

snippet.initz = function(args, covar) {
  
  let m = covar.pop / (args.S_0 + args.E_0 + args.R_0 + args.I_0);
  let S = Math.round(m * args.S_0);
  let E = Math.round(m * args.E_0);
  let I = Math.round(m * args.I_0);
  let R = Math.round(m * args.R_0);
  let H = 0;
  return {S: S, E: E, I: I, R: R, H: H};
}

snippet.dmeasure = function (data ,hiddenState, params, giveLog = 1) {
  
  let lik
  let rho = params.rho;
  let psi = params.psi;
  let H = hiddenState.H;
  let cases = data.cases;

  let tol = 1.0e-18
  let mn = rho * H;
  let v = mn * (1.0 - rho + psi * psi * mn);
  
  let modelCases = Number(cases);
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol;
    }
    if (giveLog) lik = Math.log(lik);
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

snippet.paramsMod = ["R0","amplitude","gamma","mu","sigma","rho","psi"];
snippet.paramsIc = ["S_0", "E_0", "I_0", "R_0"];
snippet.zeronames = ["H"];
snippet.statenames = ["S","E","I","R","H"];

snippet.toEst = function(params) {
  let mu = Math.log(params.mu);
  let psi = Math.log(params.psi);
  let sigma = Math.log(params.sigma);
  let gamma = Math.log(params.gamma);
  let R0 = Math.log(params.R0);
  let rho = mathLib.logit(params.rho);
  let amplitude = mathLib.logit(params.amplitude);
  let states = mathLib.toLogBarycentric([params.S_0, params.E_0, params.I_0, params.R_0]);
  //Parameters order should be the same as paramnames.
  return {R0: R0, amplitude: amplitude, gamma: gamma, mu: mu, sigma: sigma,
     rho: rho, psi: psi, S_0: states[0], E_0: states[1], I_0: states[2], R_0: states[3]};
}

snippet.fromEst = function(params) {
  let mu = Math.exp(params.mu);
  let psi = Math.exp(params.psi);
  let sigma = Math.exp(params.sigma);
  let gamma = Math.exp(params.gamma);
  let R0 = Math.exp(params.R0);
  let rho = mathLib.expit(params.rho);
  let amplitude = mathLib.expit(params.amplitude);
  let states = mathLib.fromLogBarycentric([params.S_0, params.E_0, params.I_0, params.R_0]);
  //Parameters order should be the same as paramnames.
  return {R0: R0, amplitude: amplitude, gamma: gamma, mu: mu, sigma: sigma,
    rho: rho, psi: psi, S_0: states[0], E_0: states[1], I_0: states[2], R_0: states[3]};
}

module.exports = snippet
