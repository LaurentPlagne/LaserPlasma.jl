# Dipole laser field, plasma branch of field_laser_i (ewald.c/ppbs.c).
# Zero until `n_relax` periods, linear ramp over `n_increase` periods, then
# E_x = F0 cos(ωt). The magnetic field vanishes in plasma mode.

"""
    laser_half_period(p)

Half laser period `π/ω`.
"""
laser_half_period(p::PlasmaParams) = 0.5 * (2.0 * π) / p.ω

"""
    laser_electric_field(p, t)

`E_x` of the dipole field: zero until `n_relax` periods, linear ramp over
`n_increase` periods, then `F0 cos(ωt)`.
"""
function laser_electric_field(p::PlasmaParams, t::Float64)
    t_relax = 2.0 * laser_half_period(p) * p.n_relax
    t_full = 2.0 * laser_half_period(p) * p.n_increase + t_relax
    t <= t_relax && return 0.0
    amp = t > t_full ? 1.0 : (t - t_relax) / (t_full - t_relax)
    return amp * p.F0 * cos(t * p.ω)
end

"""
    laser_fields(p, t, x, y, z)

`(Ex, Ey, Ez, Bx, By, Bz)`; B is zero in plasma mode.
"""
laser_fields(p::PlasmaParams, t::Float64, x::Float64, y::Float64, z::Float64) =
    (laser_electric_field(p, t), 0.0, 0.0, 0.0, 0.0, 0.0)
