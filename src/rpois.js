
/**
 *  @file             rpois.js   
 *                    Random number generator using the poisson distribution.    
 *                         
 *  @references       Jacob K.F. Bogers
 *                    Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com
 *                    https://github.com/R-js/libRmath.js/blob/master/src/lib/poisson/rpois.ts
 *
 *  @autor            Nazila Akhavan, nazila@kingsds.network
 *  @date             March 2019
 */

let mathLib = require('./mathLib')

rpois = {}

let M_1_SQRT_2PI = 1 / Math.sqrt(2 * Math.PI)
let a0 = -0.5;
let a1 = 0.3333333;
let a2 = -0.2500068;
let a3 = 0.2000118;
let a4 = -0.1661269;
let a5 = 0.1421878;
let a6 = -0.1384794;
let a7 = 0.125006;
let one_7 = 0.1428571428571428571;
let one_12 = 0.0833333333333333333;
let one_24 = 0.0416666666666666667;

rpois.rpoisOne = function (mu) {
    let fact = [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880];
    let l = 0;
    let m = 0;
    let pp = new Array(36);
    let b1 = 0;
    let b2 = 0;
    let c = 0;
    let c0 = 0;
    let c1 = 0;
    let c2 = 0;
    let c3 = 0;
    let p0 = 0;
    let p = 0;
    let q = 0;
    let s = 0;
    let d = 0;
    let omega = 0;
    let big_l = 0;
    let muprev = 0;
    let muprev2 = 0;
    let del;
    let difmuk = 0;
    let E = 0;
    let fk = 0;
    let fx;
    let fy;
    let g;
    let px;
    let py;
    let t = 0;
    let u = 0;
    let v;
    let x;
    let pois = -1;
    let k;
    let kflag = 0;
    let big_mu;
    let new_big_mu = false;
    if (!isFinite(mu) || mu < 0) {
        throw "error"
    }
    if (mu <= 0)
        return 0;
    big_mu = mu >= 10;
    if (big_mu) {
        new_big_mu = false;
    }
    if (!(big_mu && mu === muprev)) {
        if (big_mu) {
            new_big_mu = true;
            muprev = mu;
            s = Math.sqrt(mu);
            d = 6 * mu * mu;
            big_l = Math.floor(mu - 1.1484);
        }
        else {
            if (mu !== muprev) {
                muprev = mu;
                m = Math.trunc(Math.max(1, Math.trunc(mu)));
                l = 0;
                q = p0 = p = Math.exp(-mu);
            }
            while (true) {
                u = Math.random();
                if (u <= p0)
                    return 0;
                if (l !== 0) {
                    for (k = u <= 0.458 ? 1 : Math.trunc(Math.min(l, m)); k <= l; k++)
                        if (u <= pp[k])
                            return k;
                    if (l === 35)
                        continue;
                }
                l++;
                for (k = l; k <= 35; k++) {
                    p *= mu / k;
                    q += p;
                    pp[k] = q;
                    if (u <= q) {
                        l = k;
                        return k;
                    }
                }
                l = 35;
            }
        }
    }
    function randNormal() {
      var u = 0, v = 0;
      while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
      while(v === 0) v = Math.random();
      return Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
    }
    g =
        mu + s * randNormal(); //mathLib.normalRand()//rng.norm_randOne();
    if (g >= 0) {
        pois = Math.floor(g);
        if (pois >= big_l)
            return pois;
        fk = pois;
        difmuk = mu - fk;
        u = Math.random();
        if (d * u >= difmuk * difmuk * difmuk)
            return pois;
    }
    if (new_big_mu || mu !== muprev2) {
        muprev2 = mu;
        omega = M_1_SQRT_2PI / s;
        b1 = one_24 / mu;
        b2 = 0.3 * b1 * b1;
        c3 = one_7 * b1 * b2;
        c2 = b2 - 15 * c3;
        c1 = b1 - 6 * b2 + 45 * c3;
        c0 = 1 - b1 + 3 * b2 - 15 * c3;
        c = 0.1069 / mu;
    }
    let gotoStepF = false;
    let once = true;
    while (true) {
        if (once) {
            once = false;
            if (g >= 0) {
                kflag = 0;
                gotoStepF = true;
            }
        }
        if (!gotoStepF) {
            E = mathLib.expRand(Math.random)//sMath.exp_1.Math.exp_rand(rng.unif_rand);
            u = 2 * Math.random() - 1;
            t = 1.8 + mathLib.sign(E, u >= 0);
        }
        if (t > -0.6744 || gotoStepF) {
            if (!gotoStepF) {
                pois = Math.floor(mu + s * t);
                fk = pois;
                difmuk = mu - fk;
                kflag = 1;
            }
            gotoStepF = false;
            if (pois < 10) {
                px = -mu;
                py = Math.pow(mu, pois) / fact[Math.trunc(pois)];
            }
            else {
                del = one_12 / fk;
                del = del * (1 - 4.8 * del * del);
                v = difmuk / fk;
                if (Math.abs(v) <= 0.25)
                    px =
                        fk *
                            v *
                            v *
                            (((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v +
                                a1) *
                                v +
                                a0) -
                            del;
                else
                    px = fk * Math.log(1 + v) - difmuk - del;
                py = M_1_SQRT_2PI / Math.sqrt(fk);
            }
            x = (0.5 - difmuk) / s;
            x *= x;
            fx = -0.5 * x;
            fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
            if (kflag > 0) {
                if (c * Math.abs(u) <= py * Math.exp(px + E) - fy * Math.exp(fx + E)) {
                    break;
                }
            }
            else if (fy - u * fy <= py * Math.exp(px - fx)) {
                break;
            }
        }
    }
    return pois;
}

module.exports = rpois;

