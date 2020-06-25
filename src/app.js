/**
 * 
 *
 */
let mif2 = require('./mif2.js');
let snippet = require('./modelSnippet.js');
let fs = require('fs');
let mathLib = require('./mathLib.js');


rootDir = '..'

let dataCases = [];
let dataCasesTimes = [];
let dataCovar = [];
let dataCovarTimes = [];
let currentParams = []; 

// 1st data set; read all rows and delete last one if it is ['']
let temp, file;
file = fs.readFileSync(rootDir+'/samples/London_covar.csv').toString()
let lines = file.split('\n');
for (let i = 1; i < lines.length; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(function (x) {return Number(x)});
    dataCovarTimes.push(temp[0]);
    dataCovar.push(temp.slice(1));
  }
}

//* 2nd data set
file = fs.readFileSync(rootDir+'/samples/London_BiDataMain.csv').toString()
lines = file.split('\n');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(function (x) {return Number(x)});
    dataCasesTimes.push(temp[0]);
    dataCases.push(temp[1]);
  }
}
//* 3nd data set
file = fs.readFileSync(rootDir+'/samples/initial_parameters.csv').toString()
lines = file.split('\n');
let currentParams_name = lines[0].split(',');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(function (x) {return Number(x)});
    currentParams.push(temp);
  }
}

// need to sort curren_param based on currentParams and paramnames
let st = currentParams_name[currentParams_name.length - 1];
currentParams_name[currentParams_name.length - 1] = st.substring(0, st.length - 1);
let initialParams = new Array(currentParams.length).fill(Array(currentParams[0].length))
for(let i = 0; i < snippet.paramnames.length; i++) {
  for( let j = 0; j < currentParams_name.length; j++) {
    if(snippet.paramnames[i] === currentParams_name[j]) {
      for (let k = 0; k < currentParams.length; k++) {
        initialParams[k][i] = currentParams[k][j];
      }
    }
  }       
}


current_params = initialParams[0];
let params_ic_fit = [];
let params_mod_fit = ["R0", "amplitude", "mu", "rho", "psi"];
let cool_fraction = 0.05;
let paramnames_rw = ["R0","amplitude","mu","rho","psi", "S_0", "E_0", "I_0", "R_0"];

let param_rwIndex = mathLib.index(snippet.paramnames,paramnames_rw);//index of params that are in rw;

const rw_sd_f = function(time) {
  let rwSize = 0.05;
  let R0 = time < 1944 ? 0 : rwSize;
  let amplitude = time < 1944 ? 0 : rwSize;
  let mu = time < 1944 ? 0 : rwSize;
  let rho = time < 1944 ? 0 : rwSize;
  let psi = time < 1944 ? 0 : rwSize;
  let S_0 = time < 1944 ? 0 : rwSize;
  let E_0 = time < 1944 ? 0 : rwSize;
  let I_0 = time < 1944 ? 0 : rwSize;
  let R_0 = time < 1944 ? 0 : rwSize;
  return [R0, amplitude, mu, rho, psi, S_0, E_0, I_0, R_0];
}
///////////////////////////////////////////////////

/////////////////////////////////////
const pomp = {
  data :  dataCases,
  times:  dataCasesTimes,
  t0: 1940,
  rprocess :  { type:"euler_sim", stepFunction: snippet.rproc, deltaT: 1/365.25 },
  rmeasure: snippet.rmeas,
  covar: dataCovar,
  tcovar: dataCovarTimes,
  dmeasure: snippet.dmeasure,
  zeronames: snippet.zeronames,
  initializer: snippet.initz,
  toEstimationScale: snippet.toEst,
  fromEstimationScale: snippet.fromEst,
  statenames: snippet.statenames,
  paramnames: snippet.paramnames,
  paramnamesRw: snippet.paramnames_rw,
  coef: current_params,
}
let d1 = [], d2 = [];
for (let i = 0; i < pomp.covar.length; i++) {
  d1.push([Number(pomp.covar[i][0]), Number(pomp.covar[i][1])])
  d2.push([Number(pomp.covar[i][0]), Number(pomp.covar[i][2])])
}

pomp.population = mathLib.interpolator(d1);
pomp.birth = mathLib.interpolator(d2);
pomp.pIndex = {
  "R0": 0,
  "amplitude": 1,
  "gamma": 2,
  "mu": 3,
  "sigma": 4,
  "rho": 5,
  "psi": 6,
  "S": 7,
  "E": 8,
  "I": 9,
  "R":10,
  "H": 11
};

mif2.mif2Internal(
  {pomp: pomp,
  Nmif: 1,
  start: current_params,
  transform: true,
  ivps: params_ic_fit,
  pars: params_mod_fit,
  rw_sd: rw_sd_f,
  param_rwIndex: param_rwIndex,
  Np: 2,
  varFactor: 2,
  coolingType: "hyperbolic",
  coolingFraction: cool_fraction
  }
)


