

/** randwalk_perturbation adds random normal value to the parameters.
 * 
 * @param {matrix} params  
 * @param {array} rw_sd    random walk corspond with parametrs in params. 
 * @param {array} pidx     indices of parameters undergoing random walk
 */

const { qnorm } = require('lib-r-math.js/dist/src/lib/normal/qnorm');
exports.randwalk_perturbation = function (params, rw_sd) { 
  for(let i = 0; i < params.length; i++) {
    Object.keys(rw_sd).forEach(key => {
      if (rw_sd[key] !== null) {
        params[i][key] += rw_sd[key] * normalRand();
      } ;
    });
  }

  return params;
}

// const { qnorm5 } = require("./qnorm.js")
// INVERSION method
const normalRand = function () {

  let BIG = 134217728; /* 2^27 */
	u1 = Math.random();
  u1 = BIG * u1 + Math.random();
  
	return qnorm(u1 / BIG, 0.0, 1.0, 1, 0);
}
// let records = []
// for (let i = 0; i < 1000; i++) {
//   records.push([normalRand()])
// }

// let createCsvWriter = require('csv-writer').createArrayCsvWriter;
//   let csvWriter = createCsvWriter({
//     header: [],
//     path: '../samples/rngi.csv',
//   })
//   csvWriter.writeRecords(records)       // returns a promise
//   .then(() => {
//       console.log('...Done');
//   });
  
