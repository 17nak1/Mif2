
var mathLib = require('./mathLib')

rpois = {}



rpois.rpoisOne = function (mu) {
    var M_1_SQRT_2PI = 1 / Math.sqrt(2 * Math.PI)
    var a0 = -0.5;
    var a1 = 0.3333333;
    var a2 = -0.2500068;
    var a3 = 0.2000118;
    var a4 = -0.1661269;
    var a5 = 0.1421878;
    var a6 = -0.1384794;
    var a7 = 0.125006;
    var one_7 = 0.1428571428571428571;
    var one_12 = 0.0833333333333333333;
    var one_24 = 0.0416666666666666667; 

    var fact = [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880];
    var l = 0;
    var m = 0;
    var pp = new Array(36);
    var b1 = 0;
    var b2 = 0;
    var c = 0;
    var c0 = 0;
    var c1 = 0;
    var c2 = 0;
    var c3 = 0;
    var p0 = 0;
    var p = 0;
    var q = 0;
    var s = 0;
    var d = 0;
    var omega = 0;
    var big_l = 0;
    var muprev = 0;
    var muprev2 = 0;
    var del;
    var difmuk = 0;
    var E = 0;
    var fk = 0;
    var fx;
    var fy;
    var g;
    var px;
    var py;
    var t = 0;
    var u = 0;
    var v;
    var x;
    var pois = -1;
    var k;
    var kflag = 0;
    var big_mu;
    var new_big_mu = false;
    var gotoStepF, once
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
                u = Math.random();//rng.unif_rand()
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
    g = mu + s * Math.random() //mathLib.normalRand()//rng.norm_randOne();
    if (g >= 0) {
        pois = Math.floor(g);
        if (pois >= big_l)
            return pois;
        fk = pois;
        difmuk = mu - fk;
        u = Math.random();//rng.unif_rand()
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
    gotoStepF = false;
    once = true;
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
            u = 2 * Math.random() - 1;//rng.unif_rand()
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
// console.log(rpois.rpoisOne(10))