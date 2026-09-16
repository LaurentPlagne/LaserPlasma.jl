# Ewald summation on a table, ported from ewald.c.
#
# The table holds only the periodic *correction* V_corr = real erfc + reciprocal
# − 1/r; the softened pair potential (1/√(r²+ε²)) is added separately by the
# caller. The table is sampled on (size+1)³ nodes and read by trilinear
# interpolation. Expressions reproduce the C operation order (multiplications
# by reciprocals, `sqrt(π)` in the denominator).

"""
    FieldJacobian

Nine grids `∂F_i/∂x_j` of the Ewald correction field, in the C table order
`(xx, xy, xz, yx, yy, yz, zx, zy, zz)`. The sums are symmetric, so the triplet
`(xx, xy, xz)` is also the column `∂F/∂x` used for the dipole terms of the tree.

`packed` is the same table in a node-major layout — 16 rows per grid node
`(pot, fx, fy, fz, then the nine derivatives, padded to 16)` — which makes the
reader a rank-1 update that LLVM vectorizes (`<2 x double>` on NEON).
"""
struct FieldJacobian
    xx::Array{Float64,3}
    xy::Array{Float64,3}
    xz::Array{Float64,3}
    yx::Array{Float64,3}
    yy::Array{Float64,3}
    yz::Array{Float64,3}
    zx::Array{Float64,3}
    zy::Array{Float64,3}
    zz::Array{Float64,3}
    packed::Matrix{Float64}
end

"""
    EwaldTable

Trilinear-interpolation table of the periodic Ewald correction (potential and
field on a `(size+1)³` grid). `hess` is `nothing` for the monopole-only table
used by the direct forces, and a [`FieldJacobian`](@ref) when the tree
quadrupole/dipole corrections need the field derivatives.
"""
struct EwaldTable{H}
    h::Float64
    invh::Float64
    coef::Float64
    xmin::Float64
    size::Int
    pot::Array{Float64,3}
    fx::Array{Float64,3}
    fy::Array{Float64,3}
    fz::Array{Float64,3}
    packed::Matrix{Float64}
    hess::H
end

"""
    pack_nodes(pot, fx, fy, fz)

Node-major copy of the table: four contiguous lanes `(pot, fx, fy, fz)` per grid
node, linear node index `i + j n + k n²`.

The trilinear read combines eight nodes, and in the `(i,j,k)` layout that is 24
loads strided across three separate arrays — LLVM emits no vector instruction for it
at all (measured: 429 lines of scalar LLVM). Node-major, the same combination is a
rank-1 update on four lanes, which vectorizes to `<2 x double>` pairs: 44 vector ops,
273 lines, **×1.35** on the reader and bit-identical, since each component still
accumulates its eight terms in the same order.

This is the layout `ewald_force_jacobian_packed` already uses for the tree (16 lanes);
this is the same trick for the three-component field the direct scheme reads.
"""
function pack_nodes(pot::Array{Float64,3}, fx::Array{Float64,3}, fy::Array{Float64,3},
                    fz::Array{Float64,3})
    n = Base.size(pot, 1)
    packed = Matrix{Float64}(undef, 4, n^3)
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        idx = i + (j - 1) * n + (k - 1) * n * n
        packed[1, idx] = pot[i, j, k]
        packed[2, idx] = fx[i, j, k]
        packed[3, idx] = fy[i, j, k]
        packed[4, idx] = fz[i, j, k]
    end
    return packed
end

"""
    build_ewald_table(L; size = 16, nmax = 5, hmax = 5, alpha = 2/L)

Sample the correction `V_corr = real erfc + reciprocal − 1/r` on the grid
(potential and field only, no field derivatives).
"""
function build_ewald_table(L::Float64; size::Int = 16, nmax::Int = 5, hmax::Int = 5,
                           alpha::Float64 = 2.0 / L)
    h = L / size
    coef = 1.0 / (h * h * h)
    xmin = -L / 2
    n = size + 1
    grid = [ewald_corr_point(L, alpha, nmax, hmax, xmin + (i - 1) * h, xmin + (j - 1) * h,
                             xmin + (k - 1) * h)
            for i in 1:n, j in 1:n, k in 1:n]
    pot = getindex.(grid, 1)
    fx = getindex.(grid, 2)
    fy = getindex.(grid, 3)
    fz = getindex.(grid, 4)
    return EwaldTable(h, 1.0 / h, coef, xmin, size, pot, fx, fy, fz,
                      pack_nodes(pot, fx, fy, fz), nothing)
end

"""
    build_ewald_hessian_table(L; size = 16, nmax = 5, hmax = 5, alpha = 2/L)

Same grid as [`build_ewald_table`](@ref) plus the nine field derivatives
`∂F_i/∂x_j` (`Sum_Ewald_corr` of the C), needed by the tree dipole/quadrupole
terms. This is the C table built by `build_tab_ewald`.
"""
function build_ewald_hessian_table(L::Float64; size::Int = 16, nmax::Int = 5,
                                   hmax::Int = 5, alpha::Float64 = 2.0 / L)
    h = L / size
    coef = 1.0 / (h * h * h)
    xmin = -L / 2
    n = size + 1
    pot = Array{Float64,3}(undef, n, n, n)
    fx = similar(pot)
    fy = similar(pot)
    fz = similar(pot)
    jac = FieldJacobian(similar(pot), similar(pot), similar(pot), similar(pot),
                        similar(pot), similar(pot), similar(pot), similar(pot),
                        similar(pot), zeros(16, n * n * n))
    packed = jac.packed
    @inbounds for i in 1:n, j in 1:n, k in 1:n
        g = ewald_corr_point_hessian(L, alpha, nmax, hmax, xmin + (i - 1) * h,
                                     xmin + (j - 1) * h, xmin + (k - 1) * h)
        pot[i, j, k] = g.pot
        fx[i, j, k], fy[i, j, k], fz[i, j, k] = g.f
        jac.xx[i, j, k], jac.xy[i, j, k], jac.xz[i, j, k] = g.j[1], g.j[2], g.j[3]
        jac.yx[i, j, k], jac.yy[i, j, k], jac.yz[i, j, k] = g.j[4], g.j[5], g.j[6]
        jac.zx[i, j, k], jac.zy[i, j, k], jac.zz[i, j, k] = g.j[7], g.j[8], g.j[9]
        idx = i + (j - 1) * n + (k - 1) * n * n
        packed[1, idx] = g.pot
        packed[2, idx], packed[3, idx], packed[4, idx] = g.f
        for q in 1:9
            packed[4+q, idx] = g.j[q]
        end
    end
    return EwaldTable(h, 1.0 / h, coef, xmin, size, pot, fx, fy, fz,
                      pack_nodes(pot, fx, fy, fz), jac)
end

"""
    ewald_corr_point(L, alpha, nmax, hmax, xj, yj, zj)

Real (erfc) sum + reciprocal sum − direct `1/r` interaction.
"""
function ewald_corr_point(L, alpha, nmax, hmax, xj, yj, zj)
    p1, f1 = ewald_real(L, alpha, nmax, xj, yj, zj)
    p2, f2 = ewald_recip(L, alpha, hmax, xj, yj, zj)
    pd, fd = ewald_direct(xj, yj, zj)
    p = p1 + p2 - pd
    f = ntuple(c -> f1[c] + f2[c] - fd[c], 3)
    return (p, f...)
end

"""
    ewald_real(L, alpha, nmax, xj, yj, zj)

Real-space `erfc` sum (potential and field).
"""
function ewald_real(L, alpha, nmax, xj, yj, zj)
    epsilon2 = (1.0e-5)^2
    alpha2 = alpha * alpha
    coef1 = 2.0 * alpha / (sqrt(π))
    coef2 = (alpha * alpha * alpha) / (sqrt(π))
    pot = 0.0
    fx = 0.0
    fy = 0.0
    fz = 0.0
    for nx in -nmax:nmax
        x = xj - nx * L
        for ny in -nmax:nmax
            y = yj - ny * L
            for nz in -nmax:nmax
                z = zj - nz * L
                radius2 = max(x * x + y * y + z * z, epsilon2)
                radius = sqrt(radius2)
                erfc1 = erfc(alpha * radius)
                gauss1 = exp(-alpha2 * radius2)
                b11 = (erfc1 + coef1 * radius * gauss1) / (radius * radius2)
                pot += erfc1 / radius
                fx += x * b11
                fy += y * b11
                fz += z * b11
            end
        end
    end
    return pot, (fx, fy, fz)
end

"""
    ewald_recip(L, alpha, hmax, xj, yj, zj)

Reciprocal-space sum (potential and field), `h` vectors up to `hmax`.
"""
function ewald_recip(L, alpha, hmax, xj, yj, zj)
    lenghtm1 = 1.0 / L
    alpham1 = 1.0 / alpha
    coef0 = π * lenghtm1 * alpham1
    coef1 = -1.0 * coef0 * coef0
    coef2 = 2.0 * π * lenghtm1
    coef3 = 2.0 * lenghtm1 * lenghtm1
    coef4 = lenghtm1 / π
    pot = 0.0
    fx = 0.0
    fy = 0.0
    fz = 0.0
    for hx in -hmax:hmax
        x = xj * hx
        shx = hx * hx
        for hy in -hmax:hmax
            y = yj * hy
            shy = hy * hy
            for hz in -hmax:hmax
                hx * hx + hy * hy + hz * hz == 0 && continue
                z = zj * hz
                shxyz = shx + shy + hz * hz
                sinxyz = sin(coef2 * (x + y + z))
                cosxyz = cos(coef2 * (x + y + z))
                ah = exp(coef1 * shxyz) / shxyz
                pot += ah * cosxyz
                fx += hx * ah * sinxyz
                fy += hy * ah * sinxyz
                fz += hz * ah * sinxyz
            end
        end
    end
    return coef4 * pot, (coef3 * fx, coef3 * fy, coef3 * fz)
end

"""
    ewald_direct(xj, yj, zj)

Direct `1/r` interaction subtracted from the Ewald correction.
"""
function ewald_direct(xj, yj, zj)
    epsilon2 = (1.0e-5)^2
    radius2 = max(xj * xj + yj * yj + zj * zj, epsilon2)
    radius = sqrt(radius2)
    radiusm3 = 1.0 / (radius * radius2)
    return 1.0 / radius, (xj * radiusm3, yj * radiusm3, zj * radiusm3)
end

# --- Hessian variants (C `Sum_Ewald_*`), used by the tree table only. ---

"""
    ewald_real_hessian(L, alpha, nmax, xj, yj, zj)

Real-space `erfc` sum with the field Jacobian (`Sum_Ewald_1`).
"""
function ewald_real_hessian(L, alpha, nmax, xj, yj, zj)
    epsilon2 = (1.0e-5)^2
    alpha2 = alpha * alpha
    coef1 = 2.0 * alpha / (sqrt(π))
    coef2 = (alpha * alpha * alpha) / (sqrt(π))
    pot = 0.0
    fx = fy = fz = 0.0
    fxx = fxy = fxz = fyx = fyy = fyz = fzx = fzy = fzz = 0.0
    for nx in -nmax:nmax
        x = xj - nx * L
        for ny in -nmax:nmax
            y = yj - ny * L
            for nz in -nmax:nmax
                z = zj - nz * L
                radius2 = max(x * x + y * y + z * z, epsilon2)
                radius = sqrt(radius2)
                radius3 = radius * radius2
                radius4 = radius * radius3
                radius5 = radius * radius4
                erfc1 = erfc(alpha * radius)
                gauss1 = exp(-alpha2 * radius2)
                b1 = erfc1 + coef1 * radius * gauss1
                b2 = coef2 * gauss1
                b11 = b1 / radius3
                b12 = b1 / radius5
                b21 = b2 / radius2
                pot += erfc1 / radius
                fx += x * b11
                fy += y * b11
                fz += z * b11
                fxx += (3.0 * x * x - radius2) * b12 + 4.0 * x * x * b21
                fxy += (3.0 * x * y) * b12 + 4.0 * x * y * b21
                fxz += (3.0 * x * z) * b12 + 4.0 * x * z * b21
                fyx += (3.0 * y * x) * b12 + 4.0 * y * x * b21
                fyy += (3.0 * y * y - radius2) * b12 + 4.0 * y * y * b21
                fyz += (3.0 * y * z) * b12 + 4.0 * y * z * b21
                fzx += (3.0 * z * x) * b12 + 4.0 * z * x * b21
                fzy += (3.0 * z * y) * b12 + 4.0 * z * y * b21
                fzz += (3.0 * z * z - radius2) * b12 + 4.0 * z * z * b21
            end
        end
    end
    return pot, (fx, fy, fz),
           (fxx, fxy, fxz, fyx, fyy, fyz, fzx, fzy, fzz)
end

"""
    ewald_recip_hessian(L, alpha, hmax, xj, yj, zj)

Reciprocal-space sum with the field Jacobian (`Sum_Ewald_2`).
"""
function ewald_recip_hessian(L, alpha, hmax, xj, yj, zj)
    lenghtm1 = 1.0 / L
    alpham1 = 1.0 / alpha
    coef0 = π * lenghtm1 * alpham1
    coef1 = -1.0 * coef0 * coef0
    coef2 = 2.0 * π * lenghtm1
    coef3 = 2.0 * lenghtm1 * lenghtm1
    coef4 = lenghtm1 / π
    coef5 = -4.0 * π * lenghtm1 * lenghtm1 * lenghtm1
    pot = 0.0
    fx = fy = fz = 0.0
    fxx = fxy = fxz = fyx = fyy = fyz = fzx = fzy = fzz = 0.0
    for hx in -hmax:hmax
        x = xj * hx
        shx = hx * hx
        for hy in -hmax:hmax
            y = yj * hy
            shy = hy * hy
            for hz in -hmax:hmax
                hx * hx + hy * hy + hz * hz == 0 && continue
                z = zj * hz
                shxyz = shx + shy + hz * hz
                sinxyz = sin(coef2 * (x + y + z))
                cosxyz = cos(coef2 * (x + y + z))
                ah = exp(coef1 * shxyz) / shxyz
                pot += ah * cosxyz
                fx += hx * ah * sinxyz
                fy += hy * ah * sinxyz
                fz += hz * ah * sinxyz
                fxx += hx * hx * ah * cosxyz
                fxy += hx * hy * ah * cosxyz
                fxz += hx * hz * ah * cosxyz
                fyx += hy * hx * ah * cosxyz
                fyy += hy * hy * ah * cosxyz
                fyz += hy * hz * ah * cosxyz
                fzx += hz * hx * ah * cosxyz
                fzy += hz * hy * ah * cosxyz
                fzz += hz * hz * ah * cosxyz
            end
        end
    end
    return coef4 * pot, (coef3 * fx, coef3 * fy, coef3 * fz),
           (coef5 * fxx, coef5 * fxy, coef5 * fxz, coef5 * fyx, coef5 * fyy,
            coef5 * fyz, coef5 * fzx, coef5 * fzy, coef5 * fzz)
end

"""
    ewald_direct_hessian(xj, yj, zj)

Direct `1/r` interaction with the field Jacobian (`Sum_Ewald_direct`).
"""
function ewald_direct_hessian(xj, yj, zj)
    epsilon2 = (1.0e-5)^2
    radius2 = max(xj * xj + yj * yj + zj * zj, epsilon2)
    radius = sqrt(radius2)
    radius3 = radius * radius2
    radiusm3 = 1.0 / radius3
    radiusm5 = 1.0 / (radius3 * radius2)
    pot = 1.0 / radius
    fdx = xj * radiusm3
    fdy = yj * radiusm3
    fdz = zj * radiusm3
    return pot, (fdx, fdy, fdz),
           (3.0 * xj * xj * radiusm5 - 1.0 * radiusm3, 3.0 * xj * yj * radiusm5,
            3.0 * xj * zj * radiusm5, 3.0 * yj * xj * radiusm5,
            3.0 * yj * yj * radiusm5 - 1.0 * radiusm3, 3.0 * yj * zj * radiusm5,
            3.0 * zj * xj * radiusm5, 3.0 * zj * yj * radiusm5,
            3.0 * zj * zj * radiusm5 - 1.0 * radiusm3)
end

"""
    ewald_corr_point_hessian(L, alpha, nmax, hmax, xj, yj, zj)

Real (erfc) sum + reciprocal sum − direct `1/r`, with potential, field and field
Jacobian (`Sum_Ewald_corr`).
"""
function ewald_corr_point_hessian(L, alpha, nmax, hmax, xj, yj, zj)
    p1, f1, j1 = ewald_real_hessian(L, alpha, nmax, xj, yj, zj)
    p2, f2, j2 = ewald_recip_hessian(L, alpha, hmax, xj, yj, zj)
    pd, fd, jd = ewald_direct_hessian(xj, yj, zj)
    p = p1 + p2 - pd
    f = ntuple(c -> f1[c] + f2[c] - fd[c], 3)
    j = ntuple(c -> j1[c] + j2[c] - jd[c], 9)
    return (pot = p, f = f, j = j)
end

"""
    check_cell(tab, i, j, k)

The trilinear readers index `i+1 … i+2` in a `(size+1)³` table, so the cell index
must lie in `0:size-1`. That holds for any displacement folded by [`min_image`](@ref)
into `[-L/2, L/2)`, which is what every caller passes. Checking it once here is three
comparisons instead of the twenty-four bounds checks of the reads, and it turns a
silent out-of-range read into an error.
"""
@inline function check_cell(tab::EwaldTable, i::Int, j::Int, k::Int)
    sz = tab.size
    (0 <= i < sz && 0 <= j < sz && 0 <= k < sz) || throw(ArgumentError(
        "ewald table read outside the box at cell ($i, $j, $k): pass a minimum image"))
    return nothing
end

@inline ewald_load4(packed::Matrix{Float64}, idx::Int) =
    @inbounds (packed[1, idx], packed[2, idx], packed[3, idx], packed[4, idx])

@inline ewald_addmul4(a::NTuple{4,Float64}, v::NTuple{4,Float64}, c::Float64) =
    (a[1] + c * v[1], a[2] + c * v[2], a[3] + c * v[3], a[4] + c * v[4])

"""
    ewald_force(tab, x, y, z)

Trilinear read of the correction field, from the node-major [`pack_nodes`](@ref)
layout. The eight weights sum to `h³`; `tab.coef = 1/h³` normalizes, like
`coef_ewald` in the C.

The node order is the one the `(i,j,k)` reader used — `(i,j,k)`, `(i+1,j,k)`,
`(i,j+1,k)`, … — so every component still accumulates its eight terms in the same
order and the result is bit-identical to a strided read. This is the hot path of the
direct scheme: 92 % of the cost of one pair force is here.
"""
function ewald_force(tab::EwaldTable, x::Float64, y::Float64, z::Float64)
    h = tab.h
    invh = tab.invh
    i = trunc(Int, (x - tab.xmin) * invh)
    j = trunc(Int, (y - tab.xmin) * invh)
    k = trunc(Int, (z - tab.xmin) * invh)
    check_cell(tab, i, j, k)
    bx = x - (tab.xmin + i * h)
    by = y - (tab.xmin + j * h)
    bz = z - (tab.xmin + k * h)
    ax = h - bx
    ay = h - by
    az = h - bz
    axay = ax * ay
    bxay = bx * ay
    axby = ax * by
    bxby = bx * by
    cs = (axay * az, bxay * az, axby * az, bxby * az,
          axay * bz, bxay * bz, axby * bz, bxby * bz)
    n = tab.size + 1
    n2 = n * n
    i0 = i + j * n + k * n2 + 1
    idxs = (i0, i0 + 1, i0 + n, i0 + 1 + n,
            i0 + n2, i0 + 1 + n2, i0 + n + n2, i0 + 1 + n + n2)
    packed = tab.packed
    acc = (0.0, 0.0, 0.0, 0.0)
    @inbounds for l in 1:8
        acc = ewald_addmul4(acc, ewald_load4(packed, idxs[l]), cs[l])
    end
    coef = tab.coef
    return coef * acc[2], coef * acc[3], coef * acc[4]
end

"""
    ewald_pot(tab, x, y, z)

Trilinear read of the correction potential, same weights as [`ewald_force`](@ref).
"""
function ewald_pot(tab::EwaldTable, x::Float64, y::Float64, z::Float64)
    h = tab.h
    invh = tab.invh
    i = trunc(Int, (x - tab.xmin) * invh)
    j = trunc(Int, (y - tab.xmin) * invh)
    k = trunc(Int, (z - tab.xmin) * invh)
    check_cell(tab, i, j, k)
    bx = x - (tab.xmin + i * h)
    by = y - (tab.xmin + j * h)
    bz = z - (tab.xmin + k * h)
    ax = h - bx
    ay = h - by
    az = h - bz
    ip = i + 1
    jp = j + 1
    kp = k + 1
    @inbounds p = ax * ay * az * tab.pot[ip, jp, kp] +
                  bx * ay * az * tab.pot[ip+1, jp, kp] +
                  ax * by * az * tab.pot[ip, jp+1, kp] +
                  bx * by * az * tab.pot[ip+1, jp+1, kp] +
                  ax * ay * bz * tab.pot[ip, jp, kp+1] +
                  bx * ay * bz * tab.pot[ip+1, jp, kp+1] +
                  ax * by * bz * tab.pot[ip, jp+1, kp+1] +
                  bx * by * bz * tab.pot[ip+1, jp+1, kp+1]
    return p * tab.coef
end

"""
    ewald_force_jacobian(tab, x, y, z)

Trilinear read of the correction potential, field and field Jacobian
(`calc_ewald_sum` of the C). Same weights as [`ewald_force`](@ref); `j` is in
the C table order `(xx, xy, xz, yx, yy, yz, zx, zy, zz)`. Requires a table built
built by [`build_ewald_hessian_table`](@ref).

Coordinates must lie in `[-L/2, L/2[`: at `x = +L/2` the top index would read
past the table (the C does not check either; live callers fold with `min_image`
or `BICV` first).
"""
function ewald_force_jacobian(tab::EwaldTable{<:FieldJacobian}, x::Float64, y::Float64,
                              z::Float64)
    h = tab.h
    invh = tab.invh
    i = trunc(Int, (x - tab.xmin) * invh)
    j = trunc(Int, (y - tab.xmin) * invh)
    k = trunc(Int, (z - tab.xmin) * invh)
    check_cell(tab, i, j, k)
    bx = x - (tab.xmin + i * h)
    by = y - (tab.xmin + j * h)
    bz = z - (tab.xmin + k * h)
    ax = h - bx
    ay = h - by
    az = h - bz
    axay = ax * ay
    axby = ax * by
    bxay = bx * ay
    bxby = bx * by
    coefs = (axay * az, axay * bz, axby * az, axby * bz, bxay * az, bxay * bz,
             bxby * az, bxby * bz)
    # linear indexing and manual unrolling of the cube of 8 nodes: this is the
    # tree hot loop, and the index arithmetic / tuple indexing of a loop shows in
    # the profile
    n = tab.size + 1
    n2 = n * n
    i0 = i + j * n + k * n2 + 1
    i1 = i0 + n2
    i2 = i0 + n
    i3 = i2 + n2
    i4 = i0 + 1
    i5 = i4 + n2
    i6 = i4 + n
    i7 = i6 + n2
    c0, c1, c2, c3, c4, c5, c6, c7 = coefs
    hess = tab.hess
    @inbounds begin
        pot = c0 * tab.pot[i0] + c1 * tab.pot[i1] + c2 * tab.pot[i2] + c3 * tab.pot[i3] +
              c4 * tab.pot[i4] + c5 * tab.pot[i5] + c6 * tab.pot[i6] + c7 * tab.pot[i7]
        fx = c0 * tab.fx[i0] + c1 * tab.fx[i1] + c2 * tab.fx[i2] + c3 * tab.fx[i3] +
             c4 * tab.fx[i4] + c5 * tab.fx[i5] + c6 * tab.fx[i6] + c7 * tab.fx[i7]
        fy = c0 * tab.fy[i0] + c1 * tab.fy[i1] + c2 * tab.fy[i2] + c3 * tab.fy[i3] +
             c4 * tab.fy[i4] + c5 * tab.fy[i5] + c6 * tab.fy[i6] + c7 * tab.fy[i7]
        fz = c0 * tab.fz[i0] + c1 * tab.fz[i1] + c2 * tab.fz[i2] + c3 * tab.fz[i3] +
             c4 * tab.fz[i4] + c5 * tab.fz[i5] + c6 * tab.fz[i6] + c7 * tab.fz[i7]
        fxx = c0 * hess.xx[i0] + c1 * hess.xx[i1] + c2 * hess.xx[i2] + c3 * hess.xx[i3] +
              c4 * hess.xx[i4] + c5 * hess.xx[i5] + c6 * hess.xx[i6] + c7 * hess.xx[i7]
        fxy = c0 * hess.xy[i0] + c1 * hess.xy[i1] + c2 * hess.xy[i2] + c3 * hess.xy[i3] +
              c4 * hess.xy[i4] + c5 * hess.xy[i5] + c6 * hess.xy[i6] + c7 * hess.xy[i7]
        fxz = c0 * hess.xz[i0] + c1 * hess.xz[i1] + c2 * hess.xz[i2] + c3 * hess.xz[i3] +
              c4 * hess.xz[i4] + c5 * hess.xz[i5] + c6 * hess.xz[i6] + c7 * hess.xz[i7]
        fyx = c0 * hess.yx[i0] + c1 * hess.yx[i1] + c2 * hess.yx[i2] + c3 * hess.yx[i3] +
              c4 * hess.yx[i4] + c5 * hess.yx[i5] + c6 * hess.yx[i6] + c7 * hess.yx[i7]
        fyy = c0 * hess.yy[i0] + c1 * hess.yy[i1] + c2 * hess.yy[i2] + c3 * hess.yy[i3] +
              c4 * hess.yy[i4] + c5 * hess.yy[i5] + c6 * hess.yy[i6] + c7 * hess.yy[i7]
        fyz = c0 * hess.yz[i0] + c1 * hess.yz[i1] + c2 * hess.yz[i2] + c3 * hess.yz[i3] +
              c4 * hess.yz[i4] + c5 * hess.yz[i5] + c6 * hess.yz[i6] + c7 * hess.yz[i7]
        fzx = c0 * hess.zx[i0] + c1 * hess.zx[i1] + c2 * hess.zx[i2] + c3 * hess.zx[i3] +
              c4 * hess.zx[i4] + c5 * hess.zx[i5] + c6 * hess.zx[i6] + c7 * hess.zx[i7]
        fzy = c0 * hess.zy[i0] + c1 * hess.zy[i1] + c2 * hess.zy[i2] + c3 * hess.zy[i3] +
              c4 * hess.zy[i4] + c5 * hess.zy[i5] + c6 * hess.zy[i6] + c7 * hess.zy[i7]
        fzz = c0 * hess.zz[i0] + c1 * hess.zz[i1] + c2 * hess.zz[i2] + c3 * hess.zz[i3] +
              c4 * hess.zz[i4] + c5 * hess.zz[i5] + c6 * hess.zz[i6] + c7 * hess.zz[i7]
    end
    coef = tab.coef
    return (pot = coef * pot, f = (coef * fx, coef * fy, coef * fz),
            j = (coef * fxx, coef * fxy, coef * fxz, coef * fyx, coef * fyy,
                 coef * fyz, coef * fzx, coef * fzy, coef * fzz))
end

"""
    ewald_force_jacobian(tab, x, y, z)

Guard: the table must have been built with `hessian = true`.
"""
ewald_force_jacobian(tab::EwaldTable{Nothing}, x, y, z) =
    error("Ewald table has no Hessian; build it with build_ewald_hessian_table(L)")

@inline ewald_load16(packed::Matrix{Float64}, idx::Int) =
    @inbounds @ntuple 16 q -> packed[q, idx]

@inline ewald_addmul16(acc::NTuple{16,Float64}, v::NTuple{16,Float64}, c::Float64) =
    @ntuple 16 q -> acc[q] + c * v[q]

"""
    ewald_force_jacobian_packed(tab, x, y, z)

Same as [`ewald_force_jacobian`](@ref), read from the node-major `packed` layout:
each cube node is 16 contiguous lanes (13 components padded), so the eight-node
combination is a rank-1 update LLVM vectorizes (`<2 x double>` on NEON). Verified
bit-identical to the reference reader and ~2x faster; used by the tree forces.
"""
function ewald_force_jacobian_packed(tab::EwaldTable{<:FieldJacobian}, x::Float64,
                                     y::Float64, z::Float64)
    h = tab.h
    invh = tab.invh
    i = trunc(Int, (x - tab.xmin) * invh)
    j = trunc(Int, (y - tab.xmin) * invh)
    k = trunc(Int, (z - tab.xmin) * invh)
    check_cell(tab, i, j, k)
    bx = x - (tab.xmin + i * h)
    by = y - (tab.xmin + j * h)
    bz = z - (tab.xmin + k * h)
    ax = h - bx
    ay = h - by
    az = h - bz
    axay = ax * ay
    axby = ax * by
    bxay = bx * ay
    bxby = bx * by
    cs = (axay * az, axay * bz, axby * az, axby * bz, bxay * az, bxay * bz, bxby * az,
          bxby * bz)
    n = tab.size + 1
    n2 = n * n
    i0 = i + j * n + k * n2 + 1
    idxs = (i0, i0 + n2, i0 + n, i0 + n + n2, i0 + 1, i0 + 1 + n2, i0 + 1 + n,
            i0 + 1 + n + n2)
    acc = ntuple(_ -> 0.0, Val(16))
    packed = tab.hess.packed
    @inbounds for l in 0:7
        acc = ewald_addmul16(acc, ewald_load16(packed, idxs[l+1]), cs[l+1])
    end
    coef = tab.coef
    return (pot = coef * acc[1], f = (coef * acc[2], coef * acc[3], coef * acc[4]),
            j = (coef * acc[5], coef * acc[6], coef * acc[7], coef * acc[8],
                 coef * acc[9], coef * acc[10], coef * acc[11], coef * acc[12],
                 coef * acc[13]))
end

"""
    ewald_force_jacobian_packed(tab, x, y, z)

Guard: the table must have been built with `hessian = true`.
"""
ewald_force_jacobian_packed(tab::EwaldTable{Nothing}, x, y, z) =
    error("Ewald table has no Hessian; build it with build_ewald_hessian_table(L)")
