// randomWalk = {}
// randomWalk.rw = function(time) {
//   var rwArray = []
// R0=ifelse(time<1944,0,rw_size),
//              amplitude=ifelse(time<1944,0,rw_size),
//              mu=ifelse(time<1944,0,rw_size),
//              rho=ifelse(time<1944,0,rw_size),
//              psi=ifelse(time<1944,0,rw_size),
//              S_0=ifelse(time<1944,0,rw_size),
//              E_0=ifelse(time<1944,0,rw_size),
//              I_0=ifelse(time<1944,0,rw_size),
//              R_0=ifelse(time<1944,0,rw_size))
//   
let start = new Date()           
let a=[]
for (i=0;i<1000000;i++) {
  a.push([Math.random()])
}
console.log(new Date() - start)
var createCsvWriter = require('csv-writer').createArrayCsvWriter;
  var csvWriter = createCsvWriter({
    header: ['S'],
    path:'res.csv'
  })
  
  csvWriter.writeRecords(a)
    .then(() => {
    console.log('...predictionvar')
  })
   