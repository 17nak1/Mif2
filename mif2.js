
var trans = new Array(6).fill(3)
var rate = new Array(6)
fs = require('fs')
let fmin = require ('fmin')
var mathLib = require('./mathLib')
const libR = require('lib-r-math.js')
let rpois = require('./rpois')

var {
  Normal,
  rng: {
    LecuyerCMRG,
    normal: { BoxMuller }
  }
} = libR
const ad = new LecuyerCMRG(1234)
const { rnorm } = Normal(new BoxMuller(ad))

const R0Index = 0
const AMPLITUDE = 1
const GAMMA = 2
const MU = 3
const SIGMA = 4
const RHO = 5
const PSI = 6
const ALPHA = 7
const IOTA = 8
const SIGMASE = 9
const COHORT = 10
const S_0 = 11
const E_0 = 12
const R_0 = 13
const I_0 = 14
// data
var LondonBidata, LondonCovar, 
params = [49.05178717, 0.490067873611598, 73.05, 0.348209790270606, 45.66, 0.999999958052764, 1.07168303292765, 0.974365234, 4.835449219, 0.066357422, 0.927490234, 0.032927404, 2.46E-06, 0.967057826, 1.23E-05]
var times =[1940, 1944]
var Np = 10
var nvars = 4
var toler = 1e-17
var nlost = 0
//begin mif2
var Nmif = 30, nmif = 1
var rw_size = .05, rw = new Array(Nmif).fill(null).map(() => Array(Np).fill(rw_size))//change it to matrix
var coolFrac = 0.5 
var rwIndex = new Array(11).fill(0)
rwIndex = [R0Index, AMPLITUDE, GAMMA, MU, SIGMA, RHO, PSI, S_0, E_0, I_0, R_0] //determine by randomWalk.js



//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var dataCovar = []
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataCovar.push(lines[i].split(','))
}

//* 2nd data set
dataCases = []
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataCases.push(lines[i].split(','))
}
// console.log(dataCases[0],dataCases[dataCases.length - 2])
//* main function****************************************************************

var [R0, amplitude, gamma, mu, sigma, rho, psi, alpha, iota, sigmaSE, cohort, S0, E0, R0, I0] = params
// var estim = []
// var place = []
// // var index = new Array(12) // estimating params
// index[AMPLITUDE] = 1; index[MU] = 1; index[RHO] = 1; index[PSI] = 1;
var coolFrac = 0.5// rw2 = coolFrac * rw1
var rw_size = 0.05, delT = 0.03832991102// = 2/52 
var timeLen = dataCases.length;
// for (let i = 0; i < params.length - 1; i++) {
//   if (index[i] === 1) {
//     place.push(i)
//     estim.push(params[i])
//   }
// }
var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
let interpolPop = mathLib.interpolator (d1)
let interpolBirth = mathLib.interpolator (d2)
  
// nmif =1 
for (k = t0; k <= Number(dataCases[timeLen - 2][0]) + deltaT / 3 ; k += deltaT){
  
}  

  
  var dt = 0.1
  var va = 0, seas, dy = new Array(6).fill(0)
  // cooling type = hyperbolic; if it is geometric cmn = Math.pow(coolFrac, ((time - 1 + (nmif - 1) * timeLen) / (50 * timeLen)) 
  var s = (1 - 50 * timeLen * coolFrac) / ( coolFrac - 1)
  var cmn = (s + 1)/ (s + time + (nmif - 1) * timeLen)

  for (j = 0; j < rw.length; j++) {
    params[rwIndex[j]] += cmn * rw[j] * rnorm(1)
  }

  // console.log(params)
  var R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4],
  alpha = params[7], iota = params[8], sigmaSE = params[9], cohort = params[10] 
  var S = N[0], E = N[1], R = N[2], I = N[3], H = 0
  var pop = interpolPop(time)
  var birthrate = interpolBirth(time)
  // cohort effect
  if (Math.abs(time - Math.floor(time) - 251 / 365) < 0.5 * dt) {
    var br = cohort * birthrate / dt + (1 - cohort) * birthrate
  } else {
      br = (1 - cohort) * birthrate
  }
  // term-time seasonality
  var tt = (time - Math.floor(time)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - amplitude
  }                 
  
//   var beta = R0 * (gamma + mu) * seas// transmission rate
//   var foi = beta * Math.pow(I + iota, alpha) / pop// expected force of infection
//   var dw = mathLib.rgammawn(sigmaSE, dt)// white noise (extrademographic stochasticity)
//   rate[0] = foi * dw / dt// stochastic force of infection
//   rate[1] = mu// natural S death
//   rate[2] = sigma// rate of ending of latent stage
//   rate[3] = mu// natural E death
//   rate[4] = gamma// recovery
//   rate[5] = mu// natural I death       
//   var births = rpois(1, br * (1 - va) * dt)// Poisson births
//   // transitions between classes
//   var m = pop / (params[11] + params[12] + params[13] + params[14]),
//   S = Math.round(m * params[11]),
//   E = Math.round(m * params[12]),
//   R = Math.round(m * params[13]),
//   I = Math.round(m * params[14]),
//   W = 0,
//   C = 0,
//   N = [S, E, R, I, W, C]
// // console.log("ll",mathLib.reulermultinom(2,1,1,.1,1,[1,2], [0,0]))
//   mathLib.reulermultinom(2, S, 0, dt, 0, rate, trans)
//   mathLib.reulermultinom(2, E, 2, dt, 2, rate, trans)
//   mathLib.reulermultinom(2, I, 4, dt, 4, rate, trans)
//   dy[0] += births - trans[0] - trans[1]
//   dy[1] += trans[0] - trans[2] - trans[3]
//   dy[2] += trans[2] - trans[4] - trans[5]
//   dy[3] = pop - S - E - I
//   dy[4] += (dw - dt) / sigmaSE // standardized i.i.d. white noise
//   dy[5] += trans[4]
//   return dy
// }
// // console.log(trans,poly(params, 1944, N))


// // ODE solver
// function EulersMethod (params, covarData, delT) {
//   var rho = params[5], psi = params[6], t0 = params[15], tdata = params[16]
//   var steps = 1, arr2, arr = [], pop = interpolPop(t0)
//   var m = pop / (params[11] + params[12] + params[13] + params[14]),
//     S = Math.round(m * params[11]),
//     E = Math.round(m * params[12]),
//     R = Math.round(m * params[13]),
//     I = Math.round(m * params[14]),
//     W = 0,
//     C = 0,
//     N = [S, E, R, I, W, C]
//   for (let k = t0; k <= 1945; k += delT) {//Number(dataCases[dataCases.length - 2][0]) + delT / 3
//     N[4] = 0
//     N[5] = 0
//     if (k <= tdata && k > tdata - delT) {
//       k = tdata
//     }
//     // for (let stp = 0; stp < Np; stp++) { // steps in each time interval
//       arr2 = poly(params, k + delT, N)
//       N = N.map((a, i) => Math.round(a + arr2[i] * 1 / steps * delT))
//       // console.log(N)
//     // }
//     C = Math.round(N[5])
//     var mn = rho * C
//     var v = mn * (1.0 - rho + psi * psi * mn)
//     var tol = 1.0e-18
//     var cases = rnorm(1, mn, Math.sqrt(v) + tol)
//     if (cases > 0) {
//       cases = Math.round(cases)
//     } else {
//       cases = 0
//     }
//     if (k > tdata - 2 * delT) {
//       if (k <= tdata - delT ) {
//       k = tdata - delT
//     }
//       console.log(C)
//       arr.push([k + delT , N[0], N[1], N[2], N[3], C])
//     }
//   }
//   return arr
// }


    
      


