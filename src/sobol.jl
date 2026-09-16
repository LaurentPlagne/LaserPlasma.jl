# Numerical Recipes Sobol quasi-random sequence (sobol.c), used for the initial
# ion positions. Ported verbatim: same tables, same draw order, and Float32
# arithmetic like the C `float x[]`.

const SOBOL_MAXBIT = 30
const SOBOL_MAXDIM = 6
const SOBOL_FAC = 1.0f0 / Float32(1 << SOBOL_MAXBIT)

"""
    Sobol()

State of the Numerical Recipes Sobol quasi-random sequence (`sobol.c`).
"""
mutable struct Sobol
    n::UInt64
    iu::Matrix{UInt64}
    ix::Vector{UInt64}
end

function Sobol()
    iu = zeros(UInt64, SOBOL_MAXBIT, SOBOL_MAXDIM)
    iu[1:6, :] .= UInt64[
        1 1 1 1 1 1
        0 1 3 3 1 1
        0 0 7 3 3 5
        0 0 0 11 13 9
        0 0 0 0 25 29
        0 0 0 0 0 53
    ]
    s = Sobol(0, iu, zeros(UInt64, SOBOL_MAXDIM))
    init_sobol!(s)
    return s
end

"""
    init_sobol!(s)

Build the direction numbers of `s` (`sobseq` initialization).
"""
function init_sobol!(s::Sobol)
    mdeg = (1, 2, 3, 4, 5, 6)
    ip = (0, 1, 1, 4, 7, 19)
    for k in 1:SOBOL_MAXDIM
        for j in 1:mdeg[k]
            s.iu[j, k] <<= SOBOL_MAXBIT - j
        end
        for j in mdeg[k]+1:SOBOL_MAXBIT
            ipp = UInt64(ip[k])
            ind = s.iu[j-mdeg[k], k]
            ind ⊻= ind >> mdeg[k]
            for i in mdeg[k]:-1:2
                isodd(ipp) && (ind ⊻= s.iu[j-1, k])
                ipp >>= 1
            end
            s.iu[j, k] = ind
        end
    end
    return s
end

"""
    sobol_next!(s, Val(N) = Val(3))

Return the first `N` coordinates of the next point, like the C `x[0..n-1]`;
`Float32` arithmetic as in the C `float x[]`.
"""
function sobol_next!(s::Sobol, ::Val{N} = Val(3)) where {N}
    N <= SOBOL_MAXDIM || error("Sobol: N > MAXDIM")
    ind = s.n
    s.n += 1
    j = 0
    while j < SOBOL_MAXBIT
        iseven(ind) && break
        ind >>= 1
        j += 1
    end
    j == SOBOL_MAXBIT && error("Sobol: MAXBIT too small")
    return ntuple(Val(N)) do k
        s.ix[k] ⊻= s.iu[j+1, k]
        Float32(s.ix[k]) * SOBOL_FAC
    end
end
