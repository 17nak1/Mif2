/**
 * 
 *
 */
const mif2 = require('./mif2.js');
const snippet = require('./modelSnippet.js');
const fs = require('fs');
const mathLib = require('./mathLib.js');
const { coef } = require("./mif2Helpers.js");


rootDir = '..'

let dataCases = [];
let dataCasesTimes = [];
let dataCovar = [];
let dataCovarTimes = [];
let currentParams = []; 

// 1st data set;
let temp, file;
file = fs.readFileSync(rootDir+'/samples/London_covar.csv').toString();
lines = file.split(/\r\n|\n/);
let dataCovar_name = lines[0].split(',');
for (let i = 1; i < lines.length; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCovarTimes.push(temp[0]);
    dataCovar.push(temp.slice(1));
  }
}

//* 2nd data set
file = fs.readFileSync(rootDir+'/samples/London_BiDataMainsh.csv').toString()
lines = file.split(/\r\n|\n/);
let dataCases_name = lines[0].split(',');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCasesTimes.push(temp[0]);
    dataCases.push(temp.slice(1));
  }
}

//* 3nd data set and names
file = fs.readFileSync(rootDir+'/samples/initial_parameters.csv').toString()
lines = file.split(/\r\n|\n/);
let currentParams_name = lines[0].split(',');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(function (x) {return Number(x)});
    currentParams.push(temp);
  }
}


let sortedCurrentParams = new Array(currentParams.length).fill(Array(currentParams[0].length));
// sortedCurrentParams is sorted currentParams based on snippet.paramnames.
temp = [...snippet.paramsMod, ...snippet.paramsIc];
for(let i = 0; i < temp.length; i++) {
  for( let j = 0; j < currentParams_name.length; j++) {
    if(temp[i] === currentParams_name[j]) {
      for (let k = 0; k < currentParams.length; k++) {
        sortedCurrentParams[k][i] = currentParams[k][j];
      }
    }
  }       
}

let params_ic_fit = [];
let params_mod_fit = ["R0", "amplitude", "mu", "rho", "psi"];
let cool_fraction = 0.05;
let paramnames_rw = ["R0","amplitude","mu","rho","psi", "S_0", "E_0", "I_0", "R_0"];

let current_params = sortedCurrentParams[0];//only for this example:we need "for" loop

let param_rwIndex = mathLib.index([...snippet.paramsMod, ...snippet.paramsIc], paramnames_rw);//index of params that are in rw;

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
  rprocess :  { type:"euler_sim", stepFunction: snippet.rprocess, deltaT: 1/365.25 },
  rmeasure: snippet.rmeas,
  covar: dataCovar,
  tcovar: dataCovarTimes,
  dmeasure: snippet.dmeasure,
  zeronames: snippet.zeronames,
  initializer: snippet.initz,
  toEstimationScale: snippet.toEst,
  fromEstimationScale: snippet.fromEst,
  statenames: snippet.statenames,
  paramnames: [snippet.paramsMod, snippet.paramsIc],
  paramnamesRw: snippet.paramnames_rw,
  // coef: current_params,
}

let d1 = [], d2 = [];
for (let i = 0; i < pomp.covar.length; i++) {
  d1.push([Number(pomp.tcovar[i]), Number(pomp.covar[i][0])])
  d2.push([Number(pomp.tcovar[i]), Number(pomp.covar[i][1])])
}

pomp.population = mathLib.interpolator(d1);
pomp.birthrate = mathLib.interpolator(d2);

let t = new Date()
pomp.params = current_params;//coef

let mf = mif2.mif2Internal(
  {pomp: pomp,
  Nmif: 5,
  start: current_params,
  transform: true,
  ivps: params_ic_fit,
  pars: params_mod_fit,
  rw_sd: rw_sd_f,
  param_rwIndex: param_rwIndex,
  Np: 10,
  varFactor: 2,
  coolingType: "hyperbolic",
  coolingFraction: cool_fraction
  }
)

console.log((new Date() - t)/1000, mf.loglik, coef(mf)[0])
// let createCsvWriter = require('csv-writer').createArrayCsvWriter;
//   let csvWriter = createCsvWriter({
//     header: [],
//     path: '../samples/convRec.csv',
//     append : true
//   })
//   csvWriter.writeRecords( mf.convRec)

