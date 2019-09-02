
fs = require('fs')
let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
const mathjs = require('mathjs') 
var linearInterpolator = require('linear-interpolator/node_main')

//* Set the seed for rnorm-In R:RNGkind("L'Ecuyer-CMRG", normal.kind="Box-Muller");set.seed(1234) 
const libR = require('lib-r-math.js')
const {
  Poisson,
  rng: { MersenneTwister },
  rng: { normal: { Inversion } }
} = libR
const mt = new MersenneTwister(0)// 
const { rpois } = Poisson(new Inversion(mt))
mt.init(1234)

var {
  Normal,
  rng: {
    LecuyerCMRG,
    normal: { BoxMuller }
  }
} = libR
const ad = new LecuyerCMRG(1234)
const { rnorm } = Normal(new BoxMuller(ad))

//////////////////////////////////////////data///////////////////////////////////////
var LondonBidata = [], LondonCovar = []
var rate = new Array(6) 
var params = [2.879243e+01 , 3.669313e-01 , 7.305000e+01 , 7.020601e-03 , 4.566000e+01 , 4.920423e-01 , 1.709873e-01 , 2.560858e-02 , 1.496166e-05 , 5.523655e-08,9.743764e-01]
// var params = [2.878407e+01 , 4.034646e-01,  7.305000e+01  ,1.030760e-03  ,4.566000e+01  ,4.429944e-01  ,1.647752e-01,  4.094919e-02  ,1.912634e-05 ,7.823758e-08, 9.590316e-01]
var times =[1940, 1944]
var Np = 10
var nvars = 4
var toler = 1e-17
var nlost = 0
//begin mif2
var Nmif = 1
var rw_size = .05, rw = new Array(Nmif).fill(null).map(() => Array(Np).fill(0.05))//change it to matrix
var coolFrac = 0.5 
var rwIndex = new Array(11).fill(0)
rwIndex[7] = 1; rwIndex[8] = 1; rwIndex[9] = 1; rwIndex[10] = 1 // only states


//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonCovar.push(lines[i].split(','))
}
var dataCovar = [LondonCovar][0]
//* 2nd data set
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonBidata.push(lines[i].split(','))
}
var dataCases = [LondonBidata][0]

//* main function

var [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, R_0, I_0] = params
var [t0, tdata] = times

var deltaT = 0.03832991102// = 2/52 
var dt = 1/ 365.25
var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
var interpolPop = linearInterpolator(d1)
var interpolBirth = linearInterpolator(d2)
var timeLen = dataCases.length 

  var va = 0, seas
  var particles = []
  var sampleNum = new Array(Np)
  var predictionMean = Array(timeLen).fill(null).map(() => Array(nvars)), doPredictionVariance = [], filterMean = []
  var timeCounter = 0, ws = 0, vsq, sumsq, ess
  var doPredictionVariance = 0, doPredictionMean = 1, doFilterMean = 0, allFail = 0
  var maxFail = Infinity
  //***********************************************TIME LOOP************************************
  for (k = t0; k < Number(dataCases[dataCases.length - 2][0]) + deltaT / 3; k += deltaT){//Number(dataCases[dataCases.length - 2][0]) + deltaT / 3
    if (k <= tdata && k > tdata - deltaT) {
      k = tdata
    }
   if ( k === t0) {
     var Nlog = mathLib.toLogBarycentric([params[7], params[8], params[9], params[10]],4)
     var N = mathLib.fromLogBarycentric(Nlog, 4)
     var m = interpolPop(k) / (N[0] + N[1] + N[2] + N[3]);
     params[7] = Math.round(m * N[0]),
     params[8] = Math.round(m * N[1]),
     params[9] = Math.round(m * N[2]),
     params[10] = Math.round(m * N[3])
   //begin mif2
   var paramNoiseM = Array(params.length).fill(null).map(() => Array(nvars))
   var states = Array(Np).fill(null).map(() => Array(nvars + 1))
   var s = (1 - 50 * timeLen * coolFrac) / (coolFrac - 1)
   var cmn = (s + 1)/ (s + k + (Nmif - 1) * timeLen)
     for (cnt = 0; cnt < params.length; cnt++){
       for (np = 0; np < Np; np++){
         if(rwIndex[cnt] === 1){
           paramNoiseM[cnt][np] = params[cnt] * ( 1 + cmn * rw[Nmif - 1][np] * rnorm(1))// weighted??????????????????????????
         } else {
           paramNoiseM[cnt][np] = params[cnt]
         } 
       }   
     }
     particles = mathjs.transpose(paramNoiseM)
     // particles = new Array(Np).fill(params)
   }//endif ( k === t0)
  var pop = interpolPop(k)
  var birthrate = interpolBirth(k) 
  var tt = (k - Math.floor(k)) * 365.25
  var loglik = 0
  var lik = new Array(Np)
  var weights = [], normalWeights = []

  //****************************************PARTICLE LOOP**************************************////////////////////////////////////////////////////////////////////////////
  for (np = 0; np < Np; np++){//calc for each particle
    var trans = new Array(6).fill(0)
    var R0 = particles[np][0], amplitude = particles[np][1], gamma = particles[np][2], mu = particles[np][3], sigma = particles[np][4] 
    var S = particles[np][7], E = particles[np][8], I = particles[np][9], R = particles[np][10], H = 0
    if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
      seas = 1 + amplitude * 0.2411 / 0.7589
    } else {
      seas = 1 - amplitude
    }                 
    var beta = R0 * (gamma + mu) * (sigma + mu) * seas / sigma//seasonal transmission rate
    var foi = beta * I / pop
    rate[0] = foi//force of infection
    rate[1] = mu// natural S death
    rate[2] = sigma// rate of ending of latent stage
    rate[3] = mu// natural E death
    rate[4] = gamma// recovery
    rate[5] = mu// natural I death       
    // transitions between classes
    steps = mathLib.numEulerSteps(k, k + deltaT, dt)//14
    var del_t = 1 / steps * deltaT
    for (let stp = 0; stp < steps; stp++) { // steps in each time interval
      var births = rpois(1, birthrate * (1 - va) * del_t )// Poisson births
      mathLib.reulermultinom(2, Math.round(S), 0, del_t, 0, rate, trans)
      mathLib.reulermultinom(2, Math.round(E), 2, del_t, 2, rate, trans)
      mathLib.reulermultinom(2, Math.round(I), 4, del_t, 4, rate, trans)
      S += (births - trans[0] - trans[1])
      E += (trans[0] - trans[2] - trans[3]) 
      I += (trans[2] - trans[4] - trans[5]) 
      // E = (E > 0) ? (E + trans[0] - trans[2] - trans[3]) : 0 
      // I = (I > 0) ? (I + trans[2] - trans[4] - trans[5]) :  0
      R = pop - S - E - I
      H += trans[4] 
    }
    particles[np][7] = S
    particles[np][8] = E
    particles[np][9] = I
    particles[np][10] = R
    states[np][0] = S || 0
    states[np][1] = E || 0
    states[np][2] = I || 0
    states[np][3] = R || 0
    states[np][4] = H || 0
    //***********RESAMPLE*************
    if (k >=  tdata) {
      var rho = particles[np][5], psi = particles[np][6]
      var mn = rho * H
      var v = mn * (1.0 - rho + psi * psi * mn)
      var tol = 1.0e-18
      var modelCases = Number(dataCases[timeCounter][1])
      var likvalue = snippet.dmeasure(rho, psi, H, modelCases, giveLog = 0)
      weights.push(likvalue)
    } 
  }//end particle loop/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (k >= tdata) {
    //normalize
    let wsum = 0
    for (let i = 0; i < Np; i++) {
      wsum += weights[i]
    }
    for (let i = 0; i < Np; i++) {
      normalWeights[i] = weights[i] / wsum
    };
    // check the weights and compute sum and sum of squares
   var  w = 0, ws = 0, nlost = 0
   for (let i = 0; i < Np; i++) {
     if (normalWeights[i] > toler) {
       w += normalWeights[i];
       ws += normalWeights[i]*normalWeights[i];
     } else {      // this particle is lost
       normalWeights[i] = 0;
       nlost++;
     }
   }
   if (nlost >= Np) allFail = 1 // all particles are lost
  if (allFail) {
    loglik = Math.log(toler) // minimum log-likelihood
    ess = 0  // zero effective sample size
  } else {
    loglik = Math.log(w / Np) // mean of weights is likelihood
    ess = w * w / ws  // effective sample size
  }





// compute prediction mean
    if (doPredictionMean || doPredictionVariance) {
      for (let j = 0; j< nvars + 1;j++) {
        var sum = 0, nlost = 0
        for (let nrow =0; nrow < Np; nrow++){
          if (states[nrow][j]) {
            sum += states[nrow][j]
          } else {
            nlost++
          }
        }
        predictionMean[timeCounter][j] = sum / Np
      }
    }
// compute prediction variance
    if (doPredictionVariance) {
      for (let j = 0; j< nvars + 1;j++) {
        for (let nrow =0; nrow < Np; nrow++){
          if (states[nrow][j]) {
            vsq = states[nrow][j] - sum
            sumsq += vsq * vsq
          }
        }
        predictionVariance[timeCounter][j] = sumsq / (Np - 1) 
      }
    }
//  compute filter mean
    if (doFilterMean) {
      if (allFail) {   // unweighted average
        ws = 0
        for (let j = 0; j< nvars + 1;j++) {
          for (let nrow =0; nrow < Np; nrow++){
            if (states[nrow][j]) {
              ws += states[nrow][j]
            }
          }
        } 
        filterMean[j] = ws / Np
      } else {      // weighted average
        ws = 0
        for (let j = 0; j< nvars + 1;j++) {
          for (let nrow =0; nrow < Np; nrow++){
            if (states[nrow][j]) {
              ws += states[nrow][j] * weights[k]
            }
          }
        }
        filterMean[j] = ws / w
      }
    }

    timeCounter++
    mathLib.nosortResamp(Np, normalWeights, Np, sampleNum, 0)
    for (np = 0; np < Np; np++) { // copy the particles
      particles[np] = [].concat(particles[sampleNum[np]])
    }
  }  
}//endTime

const createCsvWriter = require('csv-writer').createArrayCsvWriter;
const csvWriter = createCsvWriter({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: './mean.csv'
})
 
csvWriter.writeRecords(predictionMean)
  .then(() => {
  console.log('...predictionMean')
})

var createCsvWriter2 = require('csv-writer').createArrayCsvWriter;
var csvWriter2 = createCsvWriter2({
  header: [],
  path: './states.csv'
})
 
csvWriter2.writeRecords(states)
  .then(() => {
  console.log('...states')
})
  
const logMeanExp = function (x) {
  var mx = Math.max(...x)
  var s = x.map((x,i) => Math.exp(x - mx))
  var q = s.reduce((a, b) => a + b, 0)
  return mx + Math.log(q / x.length)
}

const numMapSteps = function (t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}

