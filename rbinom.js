/*
 * file        rbinom.js     Random variates from the binomial distribution.
 *                         
 * References                Jacob K.F. Bogers  info@mail.jacob-bogers.com
 *                           https://github.com/R-js/libRmath.js/blob/master/src/lib/binomial/rbinom.ts
 *                           https://github.com/SurajGupta/r-source/blob/master/src/nmath/rbinom.c  
 */

rbinom = {}

rbinom. rbinomOne = function (size, pp) {
    var c = 0;
    var fm = 0;
    var npq = 0;
    var p1 = 0;
    var p2 = 0;
    var p3 = 0;
    var p4 = 0;
    var qn = 0;
    var xl = 0;
    var xll = 0;
    var xlr = 0;
    var xm = 0;
    var xr = 0;
    var psave = -1;
    var nsave = -1;
    var m = 0;
    var f;
    var f1;
    var f2;
    var u;
    var v;
    var w;
    var w2;
    var x;
    var x1;
    var x2;
    var z;
    var z2;
    var p;
    var q;
    var np;
    var g;
    var r;
    var al;
    var alv;
    var amaxp;
    var ffm;
    var ynorm;
    var i;
    var ix = 0;
    var k;
    var n;
    if (!isFinite(size)){
        throw "Input values should be finite"
    }
    r = Math.round(size)
    if (r !== size)
        return NaN
    if (!isFinite(pp) ||
        r < 0 ||
        pp < 0 ||
        pp > 1) {
        return NaN
    }
    if (r === 0 || pp === 0)
        return 0;
    if (pp === 1)
        return r;
    if (r >= Number.MAX_SAFE_INTEGER) {
        throw 'Evade overflow' + r + '> MAX_SAFE_INTEGER'
    }
    n = Math.trunc(r);
    p = Math.min(pp, 1 - pp);
    q = 1 - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);
    var gotoL_np_small = false;
    if (pp !== psave || n !== nsave) {
        psave = pp;
        nsave = n;
        if (np < 30.0) {
            qn = rbinom.R_pow_di(q, n);
            gotoL_np_small = true;
        } else {
            ffm = np + p;
            m = Math.trunc(ffm);
            fm = m;
            npq = np * q;
            p1 = Math.trunc(2.195 * Math.sqrt(npq) - 4.6 * q) + 0.5;
            xm = fm + 0.5;
            xl = xm - p1;
            xr = xm + p1;
            c = 0.134 + 20.5 / (15.3 + fm);
            al = (ffm - xl) / (ffm - xl * p);
            xll = al * (1 + 0.5 * al);
            al = (xr - ffm) / (xr * q);
            xlr = al * (1 + 0.5 * al);
            p2 = p1 * (1 + c + c);
            p3 = p2 + c / xll;
            p4 = p3 + c / xlr;
        }
    } else if (n === nsave) {
        if (np < 30.0)
            gotoL_np_small = true;
    }
    var gotoFinis = false;
    while (true && !gotoL_np_small) {
        u = Math.random() * p4; // rng.unif_rand()
        v = Math.random(); // rng.unif_rand
        if (u <= p1) {
            ix = Math.trunc(xm - p1 * v + u);
            gotoFinis = true;
            break;
        }
        if (u <= p2) {
            x = xl + (u - p1) / c;
            v = v * c + 1 - Math.abs(xm - x) / p1;
            if (v > 1 || v <= 0)
                continue;
            ix = Math.trunc(x);
        } else {
            if (u > p3) {
                ix = Math.trunc(xr - Math.log(v) / xlr);
                if (ix > n)
                    continue;
                v = v * (u - p3) * xlr;
            } else {
                ix = Math.trunc(xl + Math.log(v) / xll);
                if (ix < 0)
                    continue;
                v = v * (u - p2) * xll;
            }
        }
        k = Math.abs(ix - m);
        if (k <= 20 || k >= npq / 2 - 1) {
            f = 1;
            if (m < ix) {
                for (i = m + 1; i <= ix; i++)
                    f *= g / i - r;
            } else if (m !== ix) {
                for (i = ix + 1; i <= m; i++)
                    f /= g / i - r;
            }
            if (v <= f) {
                gotoFinis = true;
                break
            }
        } else {
            amaxp = k / npq * ((k * (k / 3 + 0.625) + 0.1666666666666) / npq + 0.5);
            ynorm = -k * k / (2.0 * npq);
            alv = Math.log(v);
            if (alv < ynorm - amaxp) {
                gotoFinis = true;
                break;
            }
            if (alv <= ynorm + amaxp) {
                x1 = ix + 1;
                f1 = fm + 1;
                z = n + 1 - fm;
                w = n - ix + 1;
                z2 = z * z;
                x2 = x1 * x1;
                f2 = f1 * f1;
                w2 = w * w;
                if (alv <=
                    xm * Math.log(f1 / x1) +
                        (n - m + 0.5) * Math.log(z / w) +
                        (ix - m) * Math.log(w * p / (x1 * q)) +
                        (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) /
                            f1 /
                            166320.0 +
                        (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) /
                            z /
                            166320.0 +
                        (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) /
                            x1 /
                            166320.0 +
                        (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) /
                            w /
                            166320) {
                    gotoFinis = true;
                    break;
                }
            }
        }
    }
    if (!gotoFinis) {
        while (true) {
            ix = 0;
            f = qn;
            u = Math.random() // rng.unif_rand
            while (true) {
                if (u < f) {
                    gotoFinis = true;
                    break;
                }
                if (ix > 110)
                    break;
                u -= f;
                ix++;
                f *= g / ix - r;
            }
            if (gotoFinis) {
                break;
            }
        }
    }
    if (psave > 0.5) {
        ix = n - ix;
    }
    return ix;
}

rbinom.R_pow_di = function (x, n) {
    var pow = 1;
    if (Number.isNaN(x))
        return x;
    if (n !== 0) {
        if (!Number.isFinite(x))
            return R_pow(x, n);
        if (n < 0) {
            n = -n;
            x = 1 / x;
        }
        while (true) {
            if (n & 1)
                pow *= x;
            if ((n >>= 1))
                x *= x;
            else
                break;
        }
    }
    return pow;
}

module.exports = rbinom
