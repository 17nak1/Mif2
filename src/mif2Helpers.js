

const cooling = function(type, fraction, ntimes) {
  switch(type){
  case "geometric":
    let factor = Math.pow(fraction, (1/50));
    return function (nt, m) {
      let alpha = Math.pow(factor, (nt / ntimes + m - 1));
      return {alpha: alpha, gamma: alpha ** 2};
    }
    
  case "hyperbolic":
    if (fraction < 1) {
      let scal = (50 * ntimes * fraction - 1) / (1 - fraction);
      return function (nt, m) {
        let alpha = (1 + scal) / (scal + nt + ntimes * (m - 1));
        return {alpha: alpha, gamma: alpha ** 2};
      }
    } else {
      return function (nt, m) {
        return {alpha: 1, gamma: 1};
      }
    }
  }
}

/**
 * Rescale the parameters.
 * @param {object} pomp 
 * @param {array} params 
 * @param {string} dir 
 */
const partrans = function (pomp, params, dir = ["fromEstimationScale","toEstimationScale"]) {
  if (!Array.isArray(params[0])) params = new Array(params);
  let transParam = [].concat(params);
  switch(dir){
  case "fromEstimationScale":
    for (let i = 0 ; i < params.length; i++) {
      transParam[i] = pomp.fromEstimationScale(params[i]);
    }
    break;
    
  case "toEstimationScale":
    for (let i = 0 ; i < params.length; i++) {
      transParam[i] = pomp.toEstimationScale(params[i]);
    }
    break;
  } 
  
  if (transParam.length === 1) return [...transParam];
  return transParam;
}


module.exports = {cooling, partrans}
