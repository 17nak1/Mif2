
// let mathLib = require('./mathLib');
/** randwalk_perturbation adds random normal value to the parameters.
 * 
 * @param {matrix} params  
 * @param {array} rw_sd    random walk corspond with parametrs in params. 
 * @param {array} pidx     indices of parameters undergoing random walk
 */
let mathLib = require('./mathLib.js');
exports.randwalk_perturbation = function (params, rw_sd, pidx) { 
  let nreps = params.length;
    
  for (let i = 0; i < nreps; i++) {
    for (let j = 0; j < pidx.length; j++) {
      if(pidx[j] !== null) {
        params[i][j] += rw_sd[pidx[j]] * mathLib.normalRand();
      }
    }
  }

  return(params);
}
