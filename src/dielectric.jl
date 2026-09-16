# Collision frequency from the dielectric theory (Dawson-Oberman / Decker),
# ported from decker.c. Two models:
#  - bare Coulomb with the η cutoff `ηmax = kT^1.5/ωp` (the "Decker model" of
#    the Diel/expew figures),
#  - soft core, with the extra `(γ K1(γ))²` factor, `γ = η bmin ωp/vth`; the
#    cutoff is carried by ε and the η integral runs over all k.
#
# The double integral over (y, η) is reduced to a 1D integral over the Bessel
# argument `a = k η` using the window `W_n(a) = ∫_{-1}^1 J_n(a y)² dy`, which is
# evaluated once per order on a phase-resolved grid. The n-sum convergence
# criterion is the C one.

"""
    DielectricTheory(kT, omega, omega_p, bmin)

Parameters of the Dawson-Oberman/Decker dielectric theory.
"""
struct DielectricTheory
    kT::Float64
    omega::Float64
    omega_p::Float64
    bmin::Float64
end

"""
    DeckerCoefficients(th, v0; softcore = false)

Precomputed coefficients of the `decker.c` integrals for quiver velocity `v0`.
"""
struct DeckerCoefficients
    frac_omega::Float64
    frac_veloc::Float64
    etamax::Float64
    bmin_omegap_vth::Float64
    softcore::Bool
end

function DeckerCoefficients(th::DielectricTheory, v0::Float64; softcore::Bool = false)
    vth = sqrt(th.kT)
    return DeckerCoefficients(th.omega_p / th.omega, v0 / vth, th.kT^1.5 / th.omega_p,
                              th.bmin * th.omega_p / vth, softcore)
end

"""
    gauss_legendre(n)

Nodes and weights of Gauss-Legendre on `[-1, 1]` (NR `gauleg`). Kept for the 2D
reference integrand used in the tests.
"""
function gauss_legendre(n::Int)
    x = zeros(n)
    w = zeros(n)
    for i in 1:(n+1)÷2
        z = cos(π * (i - 0.25) / (n + 0.5))
        pp = 1.0
        for _ in 1:100
            p1 = 1.0
            p2 = 0.0
            for j in 1:n
                p3 = p2
                p2 = p1
                p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j
            end
            pp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / pp
            abs(z - z1) < 1.0e-15 && break
        end
        x[i] = -z
        x[n+1-i] = z
        w[i] = 2.0 / ((1.0 - z * z) * pp * pp)
        w[n+1-i] = w[i]
    end
    return x, w
end

"""
    decker_f1(n, eta, y, c)

Integrand of `function_1` / `function_1_k1` (`decker.c`), for one Bessel order.
"""
@inline function decker_f1(n::Int, eta::Float64, y::Float64, c::DeckerCoefficients)
    dn = Float64(n)
    etam1 = 1.0 / eta
    arg_bess = eta * c.frac_omega * c.frac_veloc * y
    arg_exp1 = dn * etam1 / c.frac_omega
    gauss = exp(-0.5 * arg_exp1 * arg_exp1)
    bess = bessel_jn_safe(n, arg_bess)
    v = 2.0 * dn * dn * bess * bess * gauss * etam1 * etam1 * etam1
    if c.softcore
        gamma = eta * c.bmin_omegap_vth
        k1 = bessel_k1(gamma)
        v *= k1 * k1 * gamma * gamma
    end
    return v
end

"""
    decker_term_2d(n, c; ngl = 384)

2D reference integral `∫_{-1}^{1} dy ∫ dη f_n`, with tensor Gauss-Legendre. Used
only in the tests to validate the window reduction.
"""
function decker_term_2d(n::Int, c::DeckerCoefficients; ngl::Int = 384)
    xgl, wgl = gauss_legendre(ngl)
    eta_hi = c.softcore ? 8.0 / c.bmin_omegap_vth : c.etamax
    half_eta = 0.5 * eta_hi
    s = 0.0
    for (y, wy) in zip(xgl, wgl)
        inner = 0.0
        for (u, wu) in zip(xgl, wgl)
            eta = half_eta * (u + 1.0)
            inner += wu * decker_f1(n, eta, y, c)
        end
        s += wy * (half_eta * inner)
    end
    return s
end

"""
    decker_window_integrand(n, a, Ccum, k, c)

1D integrand on the Bessel-argument grid: the `y` integral is the window
`W_n(a) = 2 C_n(a)/a` with `C_n(a) = ∫_0^a J_n(t)² dt`, and `dη = da/k`.
"""
@inline function decker_window_integrand(n::Int, a::Float64, Ccum::Float64, k::Float64,
                                         c::DeckerCoefficients)
    dn = Float64(n)
    eta = a / k
    etam1 = 1.0 / eta
    arg_exp1 = dn * etam1 / c.frac_omega
    gauss = exp(-0.5 * arg_exp1 * arg_exp1)
    W = 2.0 * Ccum / a
    v = 2.0 * dn * dn * W * gauss * etam1 * etam1 * etam1 / k
    if c.softcore
        gamma = eta * c.bmin_omegap_vth
        k1 = bessel_k1(gamma)
        v *= k1 * k1 * gamma * gamma
    end
    return v
end

"""
    decker_term(n, c; hdiv = 32.0)

`∫_{-1}^{1} dy ∫ dη f_n(η, y)` through the window reduction. Composite Simpson
on a phase-resolved grid (16 points per oscillation of `J_n²`, asymptotic period
`π`); at least `4n` points so that the rise of `J_n` near `a ~ n` is resolved at
small `a`.
"""
function decker_term(n::Int, c::DeckerCoefficients; hdiv::Float64 = 32.0)
    kfac = c.frac_omega * c.frac_veloc
    kfac == 0.0 && return 0.0
    eta_hi = c.softcore ? 8.0 / c.bmin_omegap_vth : c.etamax
    a_hi = kfac * eta_hi
    # resolve both the oscillation (period π) and the Gaussian rise at small a
    # (scale n kfac/frac_omega), whichever is finer
    hres = min(π / hdiv, n * kfac / (hdiv / 2.0 * c.frac_omega))
    m = max(2 * ceil(Int, a_hi / (2.0 * hres)), 2n + 1)
    isodd(m) && (m += 1)
    H = a_hi / m          # Simpson step; C is stored at the nodes of this grid
    hf = 0.5 * H          # fine grid used only for the cumulative integral
    jf = [bessel_jn_safe(n, i * hf)^2 for i in 0:2m]
    C = zeros(m + 1)
    for i in 1:m
        C[i+1] = C[i] + H / 6.0 * (jf[2i-1] + 4.0 * jf[2i] + jf[2i+1])
    end

    # Simpson on F(a) = integrand; F(0) = 0 (the Gaussian kills the endpoint)
    s = 0.0
    for i in 1:m
        F = decker_window_integrand(n, i * H, C[i+1], kfac, c)
        s += (i == m ? 1.0 : (isodd(i) ? 4.0 : 2.0)) * F
    end
    return s * H / 3.0
end

"""
    decker_i_tilde(v0, th; softcore = false, reltol = 1e-6, nmax = 400)

Sum over the Bessel order until the C relative criterion
`|terme/somme| <= 1e-6`.
"""
function decker_i_tilde(v0::Float64, th::DielectricTheory; softcore::Bool = false,
                        reltol::Float64 = 1.0e-6, nmax::Int = 400)
    c = DeckerCoefficients(th, v0; softcore = softcore)
    n = 1
    terme = decker_term(n, c)
    somme = terme
    while n < nmax
        n += 1
        terme = decker_term(n, c)
        somme += terme
        diff = sqrt(terme * terme / (somme * somme))
        diff > reltol || break
    end
    return somme
end

"""
    decker_coll_freq(v0, th; softcore = false, kwargs...)

Collision frequency `ν_ei/ωp` for a quiver velocity `v0` (`decker.c` `coll_freq`).
"""
function decker_coll_freq(v0::Float64, th::DielectricTheory; softcore::Bool = false,
                          kwargs...)
    vth = sqrt(th.kT)
    coef1 = 2.0 * π
    coef2 = coef1^(-1.5)
    coef3 = (th.omega / th.omega_p) * (vth / v0)
    coef4 = 4.0 * π * th.omega_p / (vth * vth * vth)
    coef = coef2 * coef3 * coef3 * coef4
    return coef * decker_i_tilde(v0, th; softcore = softcore, kwargs...)
end
