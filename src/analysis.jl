# Analysis of the `deltaU` output, port of freq.c (`analyse_deltaU`).
# The first line of the file (the reference cycle, ΔU = 0) is skipped; the
# collision frequency is `ν_ei = 2 (d⟨u⟩/dt)/v0²` with `v0 = F0/ω`.

"""
    CollisionFrequency

Result of the `freq.c` analysis: collision frequency `nu`, its dispersion
`sigma`, the standard error `uncertainty`, the number of cycles and the mean
`ΔU` per cycle.
"""
struct CollisionFrequency
    nu::Float64
    sigma::Float64
    uncertainty::Float64
    ncycles::Int
    dt::Float64
    mean_du::Float64
end

"""
    collision_frequency(delta_u; nb_elec, F0, omega)

Collision frequency `ν_ei = 2 (d⟨u⟩/dt)/v0²` with `v0 = F0/ω`, from a `deltaU`
matrix `(t, ΔU, cumul)`. The first line (reference cycle, `ΔU = 0`) is skipped.
"""
function collision_frequency(delta_u::AbstractMatrix; nb_elec::Int, F0::Float64,
                             omega::Float64)
    size(delta_u, 1) >= 2 || error("deltaU needs at least two lines")
    du = delta_u[2:end, 2]
    t = delta_u[2:end, 1]
    nb = length(du)
    dtm = nb == 1 ? (t[1] - delta_u[1, 1]) : (t[end] - t[1]) / (nb - 1.0)
    dum = sum(du) / nb
    du2m = sum(abs2, du) / nb
    sigma = sqrt(du2m - dum * dum)
    v0 = F0 / omega
    nu = 2.0 * (dum / (nb_elec * dtm)) / (v0 * v0)
    sig = 2.0 * (sigma / (nb_elec * dtm)) / (v0 * v0)
    return CollisionFrequency(nu, sig, sig / sqrt(nb), nb, dtm, dum)
end
