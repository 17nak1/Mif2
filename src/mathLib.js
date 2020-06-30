
/**
 *  @file        MathLib.js
 *               Library of depentent functions  
 *
 *  @autor       Nazila Akhavan, nazila@kingsds.network
 *  @date        March 2019
 */

let mathLib = {}

let erf = require('math-erf')
let rbinom = require('./rbinom')
let libUnif = require('lib-r-math.js');
const { arrayrify } = require('lib-r-math.js/dist/src/lib/r-func');
const {
    R: { numberPrecision },
    rng: { MersenneTwister, timeseed }
} = libUnif
let U = new MersenneTwister(0) 

// Normal  distribution function
mathLib.pnorm = function (x, mu = 0, sd = 1, lower_tail = true, give_log = false) {
  if (sd < 0) {
    return NaN
  }
  let ans = 1 / 2 * (1 + erf((x - mu) / sd / Math.sqrt(2)))
  if (!lower_tail) {
    ans = 1 - ans
  }
  if (give_log) {
    ans = Math.log(ans)
  }
  return ans
}


mathLib.numEulerSteps = function(t1, t2, dt) {
  let DOUBLE_EPS = 10e-8
  let tol = Math.sqrt(DOUBLE_EPS)
  let nstep
  if (t1 >= t2) {
    dt = 0
    nstep = 0
  } else if (t1 + dt >= t2) {
    dt = t2 - t1
    nstep = 1
  } else {
    nstep = Math.ceil((t2 - t1) / dt /(1 + tol))
    dt = (t2 - t1) / nstep
  }
  return nstep
}

mathLib.numMapSteps = function (t1, t2, dt) {
  let DOUBLE_EPS = 10e-8
  let tol = Math.sqrt(DOUBLE_EPS)
  let nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 + tol))
  return (nstep > 0) ? nstep : 0
}

// Resampling function
mathLib.nosortResamp = function (nw, w, np, p, offset) {
  for (j = 1; j < nw; j++) {
   w[j] += w[j-1];
 }
  if (w[nw - 1] <= 0) {
    throw "in 'systematic_resampling': non-positive sum of weight"
  }
  let du = w[nw - 1] / np
  let u = -du * U.unif_rand()//Math.random()
  let i = 0;
  for (let j = 0; j < np; j++) {
    u += du;
    while ((u > w[i]) && (i < nw - 1)) i++;//looking for the low weight
    p[j] = i;
  }
  
  if (offset){// add offset if needed
    for (j = 0; j < np; j++) p[j] += offset;
  }
}

// The Euler-multinomial distributions
mathLib.reulermultinom = function (m = 1, size, rateAdd, dt, transAdd, rate, trans) {
  let p = 0
  let j, k
  if ((size < 0) || (dt < 0) || (Math.floor(size + 0.5) !== size)) {
    for (k = 0; k < m; k++) trans[k + transAdd] = NaN
    return 0
  }
  for (k = 0; k < m; k++) {
    if (rate[k + rateAdd] < 0.0) {
      for (j = 0; j < m; j++) trans[j + transAdd] = NaN
      return 0
    }
    p += rate[k + rateAdd]// total event rate
  }
  if (p > 0) {
    size = rbinom.rbinomOne(size, 1 - Math.exp(-p * dt))// total number of events
    if (!(isFinite(size)))
      throw 'result of binomial draw is not finite.'
    m -= 1
    for (k = 0; k < m; k++) {
      if (rate[k + rateAdd] > p) p = rate[k + rateAdd]
      trans[k + transAdd] = ((size > 0) && (p > 0)) ? rbinom.rbinomOne(size, rate[k + rateAdd] / p) : 0
      if (!(isFinite(size) && isFinite(p) && isFinite(rate[k + rateAdd]) && isFinite(trans[k + transAdd]))) {
        throw 'result of binomial draw is not finite.'
      }
      size -= trans[k + transAdd]
      p -= rate[k + rateAdd]
    }
    trans[m + transAdd] = size
  } else {
    for (k = 0; k < m; k++) trans[k + transAdd] = 0
  }
}

mathLib.interpolator = function (points) {
  let first, n = points.length - 1,
    interpolated,
    leftExtrapolated,
    rightExtrapolated;

  if (points.length === 0) {
    return function () {
      return 0
    }
  }

  if (points.length === 1) {
    return function () {
      return points[0][1]
    }
  }

  points = points.sort(function (a, b) {
    return a[0] - b[0]
  })
  first = points[0]

  leftExtrapolated = function (x) {
    let a = points[0], b = points[1];
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  interpolated = function (x, a, b) {
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  rightExtrapolated = function (x) {
    let a = points[n - 1], b = points[n];
    return b[1] + (x - b[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  return function (x) {
    let i
    if (x <= first[0]) {
      return leftExtrapolated(x)
    }
    for (i = 0; i < n; i += 1) {
      if (x > points[i][0] && x <= points[i + 1][0]) {
        return interpolated(x, points[i], points[i + 1])
      }
    }
    return rightExtrapolated(x);
  }
}

mathLib.sign = function (x, signal) {
  if (isNaN(x))
      return x
  return signal ? Math.abs(x) : -Math.abs(x);
}


mathLib.expRand = function (uniformRand) {
    let q = [
        0.6931471805599453,
        0.9333736875190459,
        0.9888777961838675,
        0.9984959252914960,
        0.9998292811061389,
        0.9999833164100727,
        0.9999985691438767,
        0.9999998906925558,
        0.9999999924734159,
        0.9999999995283275,
        0.9999999999728814,
        0.9999999999985598,
        0.9999999999999289,
        0.9999999999999968,
        0.9999999999999999,
        1.0000000000000000
    ];
    let a = 0.;
    let u = uniformRand();
    while (u <= 0. || u >= 1.)
        u = uniformRand();
    while (true) {
        u += u;
        if (u > 1.)
            break;
        a += q[0];
    }
    u -= 1.;
    if (u <= q[0])
      return a + u
    let i = 0;
    let ustar = uniformRand();
    let umin = ustar;
    do {
        ustar = uniformRand();
        if (umin > ustar)
            umin = ustar;
        i++;
    } while (u > q[i]);
    return a + umin * q[0];
}
// Normal RNG
mathLib.normalRand =function() {
    let u = 0, v = 0;
    while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
    while(v === 0) v = Math.random();
    return Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v )
}

mathLib.fromLogBarycentric = function (xN) {
  var sum = 0;
  for (let i = 0; i < xN.length; i++) {
    xN[i] = Math.exp(xN[i]);
    sum += xN[i];
  }
  for (let i = 0; i < xN.length; i++) {
    xN[i] = xN[i] / sum;
  }
  return xN;
}

mathLib.toLogBarycentric = function (xN) {
  let sum = 0;
  for (let i = 0; i < xN.length; i++) {
    sum += xN[i];
  }
  for (let i = 0; i < xN.length; i++) {
    xN[i] = Math.log(xN[i] / sum);
  }
  return xN;
}

mathLib.logit = function (p) {
  return Math.log(p / (1 - p))
}

mathLib.expit = function (x) {
  return 1 / (1 + Math.exp(-x))
}

mathLib.rgammawn = function (sigma, dt) {
  var sigmasq
  sigmasq = Math.pow(sigma, 2)
  return (sigmasq > 0) ? rgamma(1, dt / sigmasq, sigmasq) : dt
}

mathLib.logMeanExp = function (x) {
  var mx = Math.max(...x)
  var s = x.map((x,i) => Math.exp(x - mx))
  var q = s.reduce((a, b) => a + b, 0)
  return mx + Math.log(q / x.length)
}

mathLib.index = function(paramnames, a) {
  let index = new Array(paramnames.length).fill(null);
  for ( let  i =0; i < paramnames.length; i++) {
    for (let j = 0; j < a.length; j++) {
      if(paramnames[i] === a[j]) {
        index[i] = j;
        j = a.length;
      }
    }
  }
  return index;
}

mathLib.mean = function(x , w = 0) {
  let nrow = x.length;
  let ncol = x[0].length;
  let mean = new Array(ncol).fill(0);
  if (w === 0)
    w = new Array(nrow).fill(1);
  let sumw = w.reduce((a,b) => a+b, 0);
  for (let i = 0; i < ncol; i ++) {
    for (let j = 0; j < nrow; j++) {
      mean[i] += x[j][i] * w[j];
    }
    mean[i] /= sumw; 
  }  
  return mean;
}

module.exports = mathLib;
