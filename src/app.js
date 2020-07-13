/**
 * 
 *
 */
const mif2 = require('./mif2.js');
const snippet = require('./modelSnippet.js');
const fs = require('fs');
const { coef } = require("./mif2Helpers.js");
let pomp = require('./pomp.js');


rootDir = '..'

let dataCases = [];
let dataCasesTimes = [];
let dataCovar = [];
let dataCovarTimes = [];
let currentParams = []; 

// 1st data set; read all rows and delete last one if it is ['']
let temp, file;
file = fs.readFileSync(rootDir+'/samples/London_covar.csv').toString();
let lines = file.split(/\r\n|\n/);
let dataCovar_name = lines[0].replace(/['"]+/g, '').split(',');
dataCovar_name.shift();
for (let i = 1; i < lines.length; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCovarTimes.push(temp[0]);
    dataCovar.push(temp.slice(1));
  }
}

//* 2nd data set
let data;
file = fs.readFileSync(rootDir+'/samples/London_BiDataMain.csv').toString()
lines = file.split(/\r\n|\n/);
let dataCases_name = lines[0].replace(/['"]+/g, '').split(',');
dataCases_name.shift();
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCasesTimes.push(temp[0]);
    data = {};
    for(let j = 0; j < temp.length - 1; j++){
      data[dataCases_name[j]] = temp[j + 1];
    }
    dataCases.push(data)
  }
}

//* 3nd data set and names
file = fs.readFileSync(rootDir+'/samples/initial_parameters.csv').toString()
lines = file.split(/\r\n|\n/);
let currentParams_name = lines[0].replace(/['"]+/g, '').split(',');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    data = {};
    for(let j = 0; j < temp.length; j++){
      data[currentParams_name[j]] = temp[j];
    }
    currentParams.push(data)
  }
}

currentParams = currentParams[0]//Only for this example, we need loop over currentParams
let params_ic_fit = [];
let params_mod_fit = ["R0", "amplitude", "mu", "rho", "psi"];
let cool_fraction = 0.05;
// let paramnames_rw = ["R0","amplitude","mu","rho","psi", "S_0", "E_0", "I_0", "R_0"];

// let param_rwIndex = mathLib.index([...snippet.paramsMod, ...snippet.paramsIc], paramnames_rw);//index of params that are in rw;

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
  return {R0: R0, amplitude: amplitude, mu: mu, rho: rho, psi: psi, S_0: S_0, E_0: E_0, I_0: I_0, R_0};
}

const mypomp = new pomp({
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
  paramnames: snippet.paramnames,
  coef: currentParams,
  covarnames: dataCovar_name,
  obsnames: dataCases_name,
});

// const mypomp = new pomp({
//   data :  dataCases,
//   times:  dataCasesTimes,
//   t0: 1940,
//   rprocess :  { type:"euler_sim", stepFunction: snippet.rprocess, deltaT: 1/365.25 },
//   rmeasure: snippet.rmeas,
//   covar: dataCovar,
//   tcovar: dataCovarTimes,
//   dmeasure: snippet.dmeasure,
//   zeronames: snippet.zeronames,
//   initializer: snippet.initz,
//   toEstimationScale: snippet.toEst,
//   fromEstimationScale: snippet.fromEst,
//   statenames: snippet.statenames,
//   paramnames: snippet.paramnames,
//   coef: currentParams,
//   covarnames: dataCovar_name,
//   obsnames: dataCases_name,
// });

let t = new Date()
mypomp.params = currentParams;//coef

let mf = mif2.mif2Internal(
  {object: mypomp,
  Nmif: 5,
  start: currentParams,
  transform: true,
  ivps: params_ic_fit,
  pars: params_mod_fit,
  rw_sd: rw_sd_f,
  Np: 2000,
  varFactor: 2,
  coolingType: "hyperbolic",
  coolingFraction: cool_fraction
  }
)

console.log((new Date() - t)/1000, mf.loglik);

const createCsvWriter = require('csv-writer').createObjectCsvWriter;
// let header = [];
// for (let i = 0; i < Object.keys(mf.convRec[0]).length; i++) {
//   header.push({id: Object.keys(mf.convRec[0])[i]})
// }

const csvWriter = createCsvWriter({
    path: './oo.csv',
    header: [ { id: 'loglik' },
    { id: 'nfail', title: 'nfail'},
    { id: 'R0', title: 'R0'},
    { id: 'amplitude', title: 'amplitude'},
    { id: 'gamma', title: 'gamma'},
    { id: 'mu', title: 'mu'},
    { id: 'sigma', title: 'sigma'},
    { id: 'rho', title: 'rho'},
    { id: 'psi', title: 'psi'},
    { id: 'S_0', title: 'S_0'},
    { id: 'E_0', title: 'E_0'},
    { id: 'I_0', title: 'I_0'},
    { id: 'R_0', title: 'R_0'} ]
});
 
 
csvWriter.writeRecords(mf.convRec)       // returns a promise
    .then(() => {
        console.log('...Done');
    });

