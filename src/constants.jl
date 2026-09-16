# Atomic units, identical to ppbs_const.h.
const C_LIGHT = 137.0
const M_ELEC = 1.0
const Q_ELEC = -1.0
const M_ION = 1836.0 * 23.0

"""
    erfc(x)

Complementary error function. Not exported by Base; Julia embeds OpenLibm, and
we call its function directly (same implementation family as the C libm).
"""
erfc(x::Real) = ccall(:erfc, Cdouble, (Cdouble,), Float64(x))
