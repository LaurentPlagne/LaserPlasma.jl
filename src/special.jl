# Numerical Recipes special functions, ported from decker.c (Bessel J and K1)
# and dcisi.c (sine/cosine integrals). Same polynomial approximations as the
# oracle, so the historical curves can be reproduced.

"""
    bessel_j0(x)

Bessel function `J0`, NR polynomial approximation (`decker.c`).
"""
function bessel_j0(x::Float64)
    ax = abs(x)
    if ax < 8.0
        y = x * x
        ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7 +
               y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))))
        ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 +
               y * (59272.64853 + y * (267.8532712 + y * 1.0))))
        return ans1 / ans2
    end
    z = 8.0 / ax
    y = z * z
    xx = ax - 0.785398164
    ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 +
           y * (-0.2073370639e-5 + y * 0.2093887211e-6)))
    ans2 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 +
           y * (0.7621095161e-6 - y * 0.934935152e-7)))
    return sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2)
end

"""
    bessel_j1(x)

Bessel function `J1`, NR polynomial approximation (`decker.c`).
"""
function bessel_j1(x::Float64)
    ax = abs(x)
    if ax < 8.0
        y = x * x
        ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 +
               y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))))
        ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74 +
               y * (99447.43394 + y * (376.9991397 + y * 1.0))))
        return ans1 / ans2
    end
    z = 8.0 / ax
    y = z * z
    xx = ax - 2.356194491
    ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 +
           y * (0.2457520174e-5 + y * (-0.240337019e-6))))
    ans2 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 +
           y * (-0.88228987e-6 + y * 0.105787412e-6)))
    ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2)
    return x < 0.0 ? -ans : ans
end

const BESS_ACC = 40.0
const BESS_BIGNO = 1.0e10
const BESS_BIGNI = 1.0e-10

"""
    bessel_jn(n, x)

Bessel function `Jn`, NR `bessj` (requires `n >= 2`).
"""
function bessel_jn(n::Int, x::Float64)
    n >= 2 || error("bessel_jn: n < 2")
    ax = abs(x)
    ax == 0.0 && return 0.0
    if ax > n
        tox = 2.0 / ax
        bjm = bessel_j0(ax)
        bj = bessel_j1(ax)
        for j in 1:n-1
            bjp = j * tox * bj - bjm
            bjm = bj
            bj = bjp
        end
        ans = bj
    else
        tox = 2.0 / ax
        m = 2 * ((n + trunc(Int, sqrt(BESS_ACC * n))) ÷ 2)
        jsum = false
        bjp = ans = 0.0
        sum = 0.0
        bj = 1.0
        for j in m:-1:1
            bjm = j * tox * bj - bjp
            bjp = bj
            bj = bjm
            if abs(bj) > BESS_BIGNO
                bj *= BESS_BIGNI
                bjp *= BESS_BIGNI
                ans *= BESS_BIGNI
                sum *= BESS_BIGNI
            end
            jsum && (sum += bj)
            jsum = !jsum
            j == n && (ans = bjp)
        end
        sum = 2.0 * sum - bj
        ans /= sum
    end
    return x < 0.0 && isodd(n) ? -ans : ans
end

"""
    bessel_jn_safe(n, x)

`Jn(x)` for any `n >= 0`, dispatching to `J0`/`J1` for the low orders.
"""
bessel_jn_safe(n::Int, x::Float64) = n == 0 ? bessel_j0(x) :
                                     n == 1 ? bessel_j1(x) : bessel_jn(n, x)

"""
    bessel_i1(x)

Modified Bessel function `I1`, NR polynomial approximation (`decker.c`).
"""
function bessel_i1(x::Float64)
    ax = abs(x)
    if ax < 3.75
        y = x / 3.75
        y *= y
        ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934 +
              y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))))
    else
        y = 3.75 / ax
        ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 -
              y * 0.420059e-2))
        ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 +
              y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))))
        ans *= exp(ax) / sqrt(ax)
    end
    return x < 0.0 ? -ans : ans
end

"""
    bessel_k1(x)

Modified Bessel function `K1`, NR polynomial approximation (`decker.c`).
"""
function bessel_k1(x::Float64)
    if x <= 2.0
        y = x * x / 4.0
        return log(x / 2.0) * bessel_i1(x) + (1.0 / x) *
               (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 +
                y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))))
    end
    y = 2.0 / x
    return (exp(-x) / sqrt(x)) * (1.25331414 + y * (0.23498619 +
           y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2 +
           y * (0.325614e-2 + y * (-0.68245e-3)))))))
end

"""
    cisi(x)

Sine and cosine integrals `(Ci(x), Si(x))`, NR `cisi` (port of `dcisi.c`,
verbatim algorithm).
"""
function cisi(x::Float64)
    eps = 6.0e-8
    euler = 0.57721566
    maxit = 100
    piby2 = 1.5707963
    fpmin = 1.0e-30
    tmin = 2.0
    t = abs(x)
    if t == 0.0
        return -1.0 / fpmin, 0.0
    end
    if t > tmin
        b = ComplexF64(1.0, t)
        c = ComplexF64(1.0 / fpmin, 0.0)
        h = 1.0 / b
        d = h
        i = 0
        for i2 in 2:maxit
            i = i2
            a = -(i - 1) * (i - 1)
            b += 2.0
            d = 1.0 / (a * d + b)
            c = b + a / c
            del = c * d
            h *= del
            abs(real(del) - 1.0) + abs(imag(del)) < eps && break
        end
        i > maxit && error("cf failed in cisi")
        h *= ComplexF64(cos(t), -sin(t))
        ci = -real(h)
        si = piby2 + imag(h)
    elseif t < sqrt(fpmin)
        ci = 0.0 + log(t) + euler
        si = t
    else
        sum = 0.0
        sums = 0.0
        sumc = 0.0
        sign = 1.0
        fact = 1.0
        odd = true
        k = 0
        for k2 in 1:maxit
            k = k2
            fact *= t / k
            term = fact / k
            sum += sign * term
            err = term / abs(sum)
            if odd
                sign = -sign
                sums = sum
                sum = sumc
            else
                sumc = sum
                sum = sums
            end
            err < eps && break
            odd = !odd
        end
        k > maxit && error("maxits exceeded in cisi")
        ci = sumc + log(t) + euler
        si = sums
    end
    return ci, x < 0.0 ? -si : si
end

"""
    si(x)

Sine integral `Si(x)`, the quantity of the historical `Si` potential.
"""
si(x::Float64) = cisi(x)[2]
