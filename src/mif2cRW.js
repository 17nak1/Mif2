
// let mathLib = require('./mathLib');
/** randwalk_perturbation adds random normal value to the parameters.
 * 
 * @param {matrix} params  
 * @param {array} rw_sd    random walk corspond with parametrs in params. 
 * @param {array} pidx     indices of parameters undergoing random walk
 */
let mathLib = require('./mathLib.js');
const { qnorm } = require('lib-r-math.js/dist/src/lib/normal/qnorm');
exports.randwalk_perturbation = function (params, rw_sd, pidx) { 
  let nreps = params.length;
    
  for (let i = 0; i < nreps; i++) {
    for (let j = 0; j < pidx.length; j++) {
      if(pidx[j] !== null) {
        params[i][j] += rw_sd[pidx[j]] * normalRand();//TODO
      }
    }
  }

  return(params);
}

const { qnorm5 } = require("./qnorm.js")
// INVERSION method
const normalRand = function () {

  let BIG = 134217728; /* 2^27 */
	/* unif_rand() alone is not of high enough precision */
	u1 = Math.random();
	u1 = BIG * u1 + Math.random();
	return qnorm5(u1 / BIG, 0.0, 1.0, 1, 0);
}
