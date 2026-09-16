# Oracle random generators, ported verbatim.
# `dran2` is the Numerical Recipes ran2 (dodeint.c); `gasdev` is the double
# version from ewald.c (Box-Muller over `dran2`). Both share their state, as
# the C `static` variables do.

const IM1 = 2147483563
const IM2 = 2147483399
const AM = 1.0 / IM1
const IMM1 = IM1 - 1
const IA1 = 40014
const IA2 = 40692
const IQ1 = 53668
const IQ2 = 52774
const IR1 = 12211
const IR2 = 3791
const NRAN2_NTAB = 32
const NDIV = 1 + IMM1 ÷ NRAN2_NTAB
const RNMX = 1.0 - 1.2e-7

"""
    Ran2(seed = -1)

State of the Numerical Recipes `ran2` generator (`dodeint.c`).
"""
mutable struct Ran2
    idum::Int64
    idum2::Int64
    iy::Int64
    iv::Vector{Int64}
end

Ran2(seed::Integer = -1) = Ran2(Int64(seed), 0, 0, zeros(Int64, NRAN2_NTAB))

"""
    dran2!(r)

Draw one uniform deviate and advance `r`; bit-identical to the C `dran2`.
"""
function dran2!(r::Ran2)
    if r.idum <= 0
        r.idum = abs(r.idum) < 1 ? 1 : abs(r.idum)
        r.idum2 = r.idum
        for j in NRAN2_NTAB+7:-1:0
            k = r.idum ÷ IQ1
            r.idum = IA1 * (r.idum - k * IQ1) - k * IR1
            r.idum < 0 && (r.idum += IM1)
            j < NRAN2_NTAB && (r.iv[j+1] = r.idum)
        end
        r.iy = r.iv[1]
    end
    k = r.idum ÷ IQ1
    r.idum = IA1 * (r.idum - k * IQ1) - k * IR1
    r.idum < 0 && (r.idum += IM1)
    k = r.idum2 ÷ IQ2
    r.idum2 = IA2 * (r.idum2 - k * IQ2) - k * IR2
    r.idum2 < 0 && (r.idum2 += IM2)
    j = r.iy ÷ NDIV
    r.iy = r.iv[j+1] - r.idum2
    r.iv[j+1] = r.idum
    r.iy < 1 && (r.iy += IMM1)
    return min(AM * r.iy, RNMX)
end

"""
    Gasdev()

State of the Box-Muller normal generator (`gasdev` in ewald.c).
"""
mutable struct Gasdev
    iset::Bool
    gset::Float64
end

Gasdev() = Gasdev(false, 0.0)

"""
    gasdev!(g, r)

Draw one standard normal deviate from the shared `dran2` stream; paired draws
are cached in `g` exactly as in the C.
"""
function gasdev!(g::Gasdev, r::Ran2)
    if !g.iset
        local v1, v2, rsq
        while true
            v1 = 2.0 * dran2!(r) - 1.0
            v2 = 2.0 * dran2!(r) - 1.0
            rsq = v1 * v1 + v2 * v2
            (rsq < 1.0 && rsq != 0.0) && break
        end
        fac = sqrt(-2.0 * log(rsq) / rsq)
        g.gset = v1 * fac
        g.iset = true
        return v2 * fac
    end
    g.iset = false
    return g.gset
end
