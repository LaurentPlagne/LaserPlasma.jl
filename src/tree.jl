# Barnes-Hut octree, ported from load.c. Node storage is a structure of arrays
# (`pos[:, i]`, `charge[i]`, ...): the hot loops (walk, moments, forces) read
# contiguous per-field columns instead of chasing `Vector{Node}` references, and
# nodes are reused across rebuilds (C free list via `newtree!`/`makecell!`).
# Node identifiers are columns: bodies 1..nbody (electrons then ions, as the C
# `mkbody_tab_rec`), cells appended during construction; 0 means "none".
# `Q` is the C `usequad` (0, 1 or 2) as a type parameter, so the kernel branches
# disappear at compile time.

# Node type tags, as the C `BODY`/`CELL` (exposed for the tests and the walk).
const BODY_T = UInt8(1)
const CELL_T = UInt8(2)

"""
    Octree{Q}

Barnes-Hut octree workspace. `Q` is the C `usequad` (0: monopole, 1: +dipole,
2: +quadrupole and Ewald Hessian). Node fields are separate arrays; capacities
are the second dimension / length of the arrays, and `cellused` counts the live
cells, so `nbody + cellused` is the number of live nodes.
"""
mutable struct Octree{Q}
    nbody::Int
    ntype::Vector{UInt8}
    mass::Vector{Float64}
    charge::Vector{Float64}
    pos::Matrix{Float64}
    next::Vector{Int}
    more::Vector{Int}
    subp::Matrix{Int}
    rcrit2::Vector{Float64}
    dip::Matrix{Float64}
    quad::Matrix{Float64}
    root::Int
    rsize::Float64
    theta::Float64
    sw93::Bool
    maxlevel::Int
    cellused::Int
end

function _octree(::Val{Q}, nbody::Int, rsize::Float64, theta::Float64,
                 sw93::Bool) where {Q}
    cap = max(4 * nbody, 64)
    return Octree{Q}(nbody, fill(BODY_T, cap), zeros(cap), zeros(cap), zeros(3, cap),
                     zeros(Int, cap), zeros(Int, cap), zeros(Int, 8, cap), zeros(cap),
                     zeros(3, cap), zeros(9, cap), 0, rsize, theta, sw93, 0, 0)
end

"""
    Octree(nbody; rsize, theta = 0.5, sw93 = true, usequad = 2)

`rsize` is the side of the root cell; as in the C it is carried by the tree and
only grows across rebuilds (`expandbox`), so the caller must initialise it to the
C value for plasma runs: `lenght_cell = L` (set by `init_plasma`; the
`4 N_ion^{1/3} r_ws` of `init_part` is the cluster path, never used here).
"""
function Octree(nbody::Int; rsize::Float64, theta::Float64 = 0.5, sw93::Bool = true,
                usequad::Int = 2)
    (usequad == 0 || usequad == 1 || usequad == 2) ||
        error("Octree: usequad must be 0, 1 or 2")
    rsize > 0.0 || error("Octree: rsize must be positive (expandbox doubles it)")
    if usequad == 2
        return _octree(Val(2), nbody, rsize, theta, sw93)
    elseif usequad == 1
        return _octree(Val(1), nbody, rsize, theta, sw93)
    end
    return _octree(Val(0), nbody, rsize, theta, sw93)
end

# Grow every backing array when the cell count exceeds the allocated capacity.
function ensure_capacity!(t::Octree, n::Int)
    cap = size(t.pos, 2)
    n <= cap && return t
    newcap = max(n, 2 * cap)
    extra = newcap - cap
    resize!(t.ntype, newcap)
    resize!(t.mass, newcap)
    resize!(t.charge, newcap)
    t.pos = hcat(t.pos, zeros(3, extra))
    resize!(t.next, newcap)
    resize!(t.more, newcap)
    t.subp = hcat(t.subp, zeros(Int, 8, extra))
    resize!(t.rcrit2, newcap)
    t.dip = hcat(t.dip, zeros(3, extra))
    t.quad = hcat(t.quad, zeros(9, extra))
    return t
end

"""
    set_bodies!(tree, mass, charge, pos)

Fill the body nodes (no tree is built): `mass` and `charge` are length-`nbody`
vectors, `pos` a `3×nbody` matrix. Mirrors `mkbody_tab_rec`: bodies are the
electrons then the ions.
"""
function set_bodies!(t::Octree, mass::AbstractVector{Float64},
                     charge::AbstractVector{Float64}, pos::AbstractMatrix{Float64})
    (length(mass) == t.nbody && length(charge) == t.nbody &&
     size(pos, 2) == t.nbody) || error("set_bodies!: expected $(t.nbody) bodies")
    @inbounds for i in 1:t.nbody
        t.ntype[i] = BODY_T
        t.mass[i] = mass[i]
        t.charge[i] = charge[i]
        t.pos[1, i] = pos[1, i]
        t.pos[2, i] = pos[2, i]
        t.pos[3, i] = pos[3, i]
        t.next[i] = 0
        t.more[i] = 0
    end
    return t
end

"""
    set_bodies!(tree, mass, charge, epos, ipos)

Fill the body nodes from separate electron and ion position matrices (electrons
first, ions next), as `mkbody_tab_rec` does at every stage.
"""
function set_bodies!(t::Octree, mass::AbstractVector{Float64},
                     charge::AbstractVector{Float64}, epos::AbstractMatrix{Float64},
                     ipos::AbstractMatrix{Float64})
    ne = size(epos, 2)
    ni = size(ipos, 2)
    ne + ni == t.nbody || error("set_bodies!: expected $(t.nbody) bodies")
    @inbounds for i in 1:ne
        t.ntype[i] = BODY_T
        t.mass[i] = mass[i]
        t.charge[i] = charge[i]
        t.pos[1, i] = epos[1, i]
        t.pos[2, i] = epos[2, i]
        t.pos[3, i] = epos[3, i]
        t.next[i] = 0
        t.more[i] = 0
    end
    @inbounds for j in 1:ni
        i = ne + j
        t.ntype[i] = BODY_T
        t.mass[i] = mass[i]
        t.charge[i] = charge[i]
        t.pos[1, i] = ipos[1, j]
        t.pos[2, i] = ipos[2, j]
        t.pos[3, i] = ipos[3, j]
        t.next[i] = 0
        t.more[i] = 0
    end
    return t
end

"""
    maketree!(tree)

Build the tree from the current body positions (C `maketree`): root cell, box
expansion, insertion, centre of mass and SW93 radii, threading, then dipole and
quadrupole moments according to the `Q` parameter. Cells are dropped first
(`newtree`).
"""
function maketree!(t::Octree{Q}) where {Q}
    newtree!(t)
    t.root = makecell!(t)
    expandbox!(t)
    t.maxlevel = 0
    @inbounds for i in 1:t.nbody
        t.mass[i] == 0.0 && continue       # exclude test particles, as the C
        loadbody!(t, i)
    end
    hackcofm!(t, t.root, t.rsize)
    threadtree!(t, t.root, 0)
    Q > 0 && hackdip!(t, t.root)
    Q > 1 && hackquad!(t, t.root)
    return t
end

"""
    build_tree!(tree, mass, charge, pos)

`set_bodies!` followed by `maketree!`.
"""
function build_tree!(t::Octree, mass::AbstractVector{Float64},
                     charge::AbstractVector{Float64}, pos::AbstractMatrix{Float64})
    set_bodies!(t, mass, charge, pos)
    return maketree!(t)
end

"""
    newtree!(tree)

Flush the cells: keep the body nodes, reset the root and the cell count. Like the
C free list (`newtree`), the cell slots are kept and reused by `makecell!`, so a
rebuild allocates nothing once the high-water mark of cells is reached.
"""
function newtree!(t::Octree)
    t.root = 0
    t.cellused = 0
    return t
end

"""
    makecell!(tree)

Return a reset cell, reusing the existing slot when possible (C `makecell`).
"""
function makecell!(t::Octree)
    t.cellused += 1
    idx = t.nbody + t.cellused
    ensure_capacity!(t, idx)
    @inbounds begin
        t.ntype[idx] = CELL_T
        t.mass[idx] = 0.0
        t.charge[idx] = 0.0
        t.pos[1, idx] = 0.0
        t.pos[2, idx] = 0.0
        t.pos[3, idx] = 0.0
        t.next[idx] = 0
        t.more[idx] = 0
        for k in 1:8
            t.subp[k, idx] = 0
        end
        t.rcrit2[idx] = 0.0
        t.dip[1, idx] = 0.0
        t.dip[2, idx] = 0.0
        t.dip[3, idx] = 0.0
        for k in 1:9
            t.quad[k, idx] = 0.0
        end
    end
    return idx
end

"""
    expandbox!(tree)

Double the root cell until it fits every body (exact powers of two, as the C);
`tree.rsize` only grows and persists across rebuilds.
"""
function expandbox!(t::Octree)
    xyzmax = 0.0
    @inbounds for i in 1:t.nbody
        xyzmax = max(xyzmax, abs(t.pos[1, i]), abs(t.pos[2, i]), abs(t.pos[3, i]))
    end
    while t.rsize < 2.0 * xyzmax
        t.rsize = 2.0 * t.rsize
    end
    return t
end

"""
    bicv1(d, L)

Fold one component into `(-L/2, L/2]`, exactly as the C macro `BICV` (`> L/2`
subtracts, `<= -L/2` adds).
"""
@inline function bicv1(d::Float64, L::Float64)
    half = 0.5 * L
    mhalf = -half
    while d > half
        d -= L
    end
    while d <= mhalf
        d += L
    end
    return d
end

"""
    subindex(p, c)

Octant of `p` in the cell centred at `c`: bit `k` set when `p[k] >= c[k]`
(returns `0..7`).
"""
@inline function subindex(px::Float64, py::Float64, pz::Float64, cx::Float64,
                          cy::Float64, cz::Float64)
    ind = 0
    ind += (cx <= px) ? 4 : 0
    ind += (cy <= py) ? 2 : 0
    ind += (cz <= pz) ? 1 : 0
    return ind
end

"""
    loadbody!(tree, ip)

Descend the tree and insert body `ip`; split a cell when it lands on a body
(C `loadbody`).
"""
function loadbody!(t::Octree, ip::Int)
    px = t.pos[1, ip]
    py = t.pos[2, ip]
    pz = t.pos[3, ip]
    q = t.root
    qind = subindex(px, py, pz, t.pos[1, q], t.pos[2, q], t.pos[3, q])
    qsize = t.rsize
    lev = 0
    @inbounds while t.subp[qind+1, q] != 0
        qs = t.subp[qind+1, q]
        if t.ntype[qs] == BODY_T
            c = makecell!(t)
            t.pos[1, c] = t.pos[1, q] + (px < t.pos[1, q] ? -qsize : qsize) / 4
            t.pos[2, c] = t.pos[2, q] + (py < t.pos[2, q] ? -qsize : qsize) / 4
            t.pos[3, c] = t.pos[3, q] + (pz < t.pos[3, q] ? -qsize : qsize) / 4
            sub = subindex(t.pos[1, qs], t.pos[2, qs], t.pos[3, qs], t.pos[1, c],
                           t.pos[2, c], t.pos[3, c]) + 1
            t.subp[sub, c] = qs
            t.subp[qind+1, q] = c
        end
        q = t.subp[qind+1, q]
        qind = subindex(px, py, pz, t.pos[1, q], t.pos[2, q], t.pos[3, q])
        qsize = qsize / 2
        lev += 1
    end
    t.subp[qind+1, q] = ip
    t.maxlevel = max(t.maxlevel, lev)
    return t
end

"""
    setrcrit!(tree, ip, cx, cy, cz, psize)

Critical squared radius of cell `ip` from its centre of mass `(cx, cy, cz)`, its
geometric centre and size (SW93 criterion of C `setrcrit`); the C also supports
bh86, unused by the archived runs.
"""
function setrcrit!(t::Octree, ip::Int, cx::Float64, cy::Float64, cz::Float64,
                   psize::Float64)
    if t.theta == 0.0
        rc = 2.0 * t.rsize
    elseif t.sw93
        dmin = cx - (t.pos[1, ip] - psize / 2)
        bmax2 = max(dmin, psize - dmin)^2
        dmin = cy - (t.pos[2, ip] - psize / 2)
        bmax2 += max(dmin, psize - dmin)^2
        dmin = cz - (t.pos[3, ip] - psize / 2)
        bmax2 += max(dmin, psize - dmin)^2
        rc = sqrt(bmax2) / t.theta
    else
        error("setrcrit!: only sw93 and theta = 0 are ported (bh86 is unused)")
    end
    t.rcrit2[ip] = rc * rc
    return t
end

"""
    hackcofm!(tree, ip, psize)

Accumulate mass, charge and centre of mass, then set the critical radius and
overwrite the cell position with the centre of mass (C `hackcofm`).
"""
function hackcofm!(t::Octree, ip::Int, psize::Float64)
    @inbounds begin
        t.mass[ip] = 0.0
        t.charge[ip] = 0.0
        cx = cy = cz = 0.0
        for i in 1:8
            q = t.subp[i, ip]
            q == 0 && continue
            t.ntype[q] == CELL_T && hackcofm!(t, q, psize / 2)
            m = t.mass[q]
            t.charge[ip] += t.charge[q]
            t.mass[ip] += m
            cx += m * t.pos[1, q]
            cy += m * t.pos[2, q]
            cz += m * t.pos[3, q]
        end
        m = t.mass[ip]
        cx /= m
        cy /= m
        cz /= m
        for k in 1:3
            c = k == 1 ? cx : k == 2 ? cy : cz
            if c < t.pos[k, ip] - psize / 2 || t.pos[k, ip] + psize / 2 <= c
                error("hackcofm!: tree structure error")
            end
        end
        setrcrit!(t, ip, cx, cy, cz, psize)
        t.pos[1, ip] = cx
        t.pos[2, ip] = cy
        t.pos[3, ip] = cz
    end
    return t
end

"""
    hackcofm_nl!(tree, ip, psize)

Refresh mass, charge, centre of mass and critical radius of cell `ip` from the
current subcells and body positions, without rebuilding the topology and without
the bounds check (C `hackcofm_nl`, used at every RK4 stage in tree mode).
"""
function hackcofm_nl!(t::Octree, ip::Int, psize::Float64)
    @inbounds begin
        t.mass[ip] = 0.0
        t.charge[ip] = 0.0
        cx = cy = cz = 0.0
        for i in 1:8
            q = t.subp[i, ip]
            q == 0 && continue
            t.ntype[q] == CELL_T && hackcofm_nl!(t, q, psize / 2)
            m = t.mass[q]
            t.charge[ip] += t.charge[q]
            t.mass[ip] += m
            cx += m * t.pos[1, q]
            cy += m * t.pos[2, q]
            cz += m * t.pos[3, q]
        end
        m = t.mass[ip]
        cx /= m
        cy /= m
        cz /= m
        setrcrit!(t, ip, cx, cy, cz, psize)
        t.pos[1, ip] = cx
        t.pos[2, ip] = cy
        t.pos[3, ip] = cz
    end
    return t
end

"""
    refresh_moments!(tree)

Recompute the cell moments from the current body positions and the frozen
topology (C `hackcofm_nl` + `hackdip_nl` + `hackquad_nl`), as `euler1_rec` does
before each RK4 stage in tree mode.
"""
function refresh_moments!(t::Octree{Q}) where {Q}
    hackcofm_nl!(t, t.root, t.rsize)
    Q > 0 && hackdip!(t, t.root)
    Q > 1 && hackquad!(t, t.root)
    return t
end

"""
    threadtree!(tree, ip, nxt)

Thread the subtree rooted at `ip` with `more`/`next` links (C `threadtree`),
without the C stack array: the children are linked in subcell order, the last
one to `nxt`.
"""
function threadtree!(t::Octree, ip::Int, nxt::Int)
    t.next[ip] = nxt
    t.ntype[ip] == CELL_T || return t
    first_child = 0
    prev = 0
    @inbounds for i in 1:8
        q = t.subp[i, ip]
        q == 0 && continue
        first_child == 0 && (first_child = q)
        prev != 0 && (t.next[prev] = q)
        prev = q
    end
    prev != 0 && (t.next[prev] = nxt)
    t.more[ip] = first_child
    @inbounds for i in 1:8
        q = t.subp[i, ip]
        q == 0 && continue
        threadtree!(t, q, t.next[q])
    end
    return t
end

"""
    hackdip!(tree, ip)

Dipole moments, accumulated bottom-up (C `hackdip`).
"""
function hackdip!(t::Octree, ip::Int)
    @inbounds begin
        t.dip[1, ip] = 0.0
        t.dip[2, ip] = 0.0
        t.dip[3, ip] = 0.0
        for i in 1:8
            q = t.subp[i, ip]
            q == 0 && continue
            t.ntype[q] == CELL_T && hackdip!(t, q)
            dx = t.pos[1, q] - t.pos[1, ip]
            dy = t.pos[2, q] - t.pos[2, ip]
            dz = t.pos[3, q] - t.pos[3, ip]
            if t.ntype[q] == CELL_T
                t.dip[1, ip] += t.dip[1, q]
                t.dip[2, ip] += t.dip[2, q]
                t.dip[3, ip] += t.dip[3, q]
            end
            c = t.charge[q]
            t.dip[1, ip] += c * dx
            t.dip[2, ip] += c * dy
            t.dip[3, ip] += c * dz
        end
    end
    return t
end

# Identity of the quadrupole stencil, as a 9-tuple (row-major).
const QUAD_ID9 = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

"""
    hackquad!(tree, ip)

Quadrupole moments, accumulated bottom-up (C `hackquad`): for each child,
`Q_p += q (3 dr⊗dr − dr² I) + Q_q`. The nine-component update is built with
`@ntuple` so LLVM can vectorize it.
"""
function hackquad!(t::Octree, ip::Int)
    @inbounds begin
        for k in 1:9
            t.quad[k, ip] = 0.0
        end
        for i in 1:8
            q = t.subp[i, ip]
            q == 0 && continue
            t.ntype[q] == CELL_T && hackquad!(t, q)
            dx = t.pos[1, q] - t.pos[1, ip]
            dy = t.pos[2, q] - t.pos[2, ip]
            dz = t.pos[3, q] - t.pos[3, ip]
            drsq = dx * dx + dy * dy + dz * dz
            c = t.charge[q]
            drr = (dx * dx, dx * dy, dx * dz, dy * dx, dy * dy, dy * dz, dz * dx,
                   dz * dy, dz * dz)
            tm = @ntuple 9 ab -> (3.0 * drr[ab] - QUAD_ID9[ab] * drsq) * c
            if t.ntype[q] == CELL_T
                tm = @ntuple 9 ab -> tm[ab] + t.quad[ab, q]
            end
            for ab in 1:9
                t.quad[ab, ip] += tm[ab]
            end
        end
    end
    return t
end

"""
    walk_interaction_list!(arena, n, tree, ip, L)

Walk the threaded tree from the root and write into `arena`, starting after the
first `n` entries, the nodes accepted for body `ip` (C `hackgrav_s`/`treescan_s`
with `subdivp_s`): a cell is opened when `|BICV(Pos(q) - pos0)|² < rcrit2(q)`,
otherwise it is accepted; the body itself is skipped. Returns the new count;
grows `arena` (by doubling) only when the cursor would overrun it.
"""
function walk_interaction_list!(arena::Vector{Int}, n::Int, t::Octree, ip::Int,
                                L::Float64)
    px = t.pos[1, ip]
    py = t.pos[2, ip]
    pz = t.pos[3, ip]
    q = t.root
    @inbounds while q != 0
        if t.ntype[q] == CELL_T
            # the C `subdivp_s` takes Pos(q) − pos0 (only the square is used)
            dx = bicv1(t.pos[1, q] - px, L)
            dy = bicv1(t.pos[2, q] - py, L)
            dz = bicv1(t.pos[3, q] - pz, L)
            if dx * dx + dy * dy + dz * dz < t.rcrit2[q]
                q = t.more[q]
                continue
            end
        end
        if q != ip
            n += 1
            if n > length(arena)
                resize!(arena, max(n, 2 * length(arena)))
            end
            arena[n] = q
        end
        q = t.next[q]
    end
    return n
end

"""
    append_interaction_list!(arena, tree, ip, L)

Append the walk of [`walk_interaction_list!`](@ref) to `arena` and return it.
"""
function append_interaction_list!(arena::Vector{Int}, t::Octree, ip::Int, L::Float64)
    n = walk_interaction_list!(arena, length(arena), t, ip, L)
    resize!(arena, n)
    return arena
end

"""
    fill_interaction_list!(list, tree, ip, L)

Empty `list`, then [`append_interaction_list!`](@ref) into it.
"""
fill_interaction_list!(list::Vector{Int}, t::Octree, ip::Int, L::Float64) =
    append_interaction_list!(empty!(list), t, ip, L)

"""
    build_interaction_list!(tree, ip, L)

Same as [`fill_interaction_list!`](@ref) into a fresh vector.
"""
build_interaction_list!(t::Octree, ip::Int, L::Float64) =
    append_interaction_list!(Int[], t, ip, L)

"""
    gravsub_plasma!(ax, ay, az, phi0, tree, q, px, py, pz, tab, L, eps_tree)

Scalar accumulation of the contribution of node `q` to the tree acceleration and
potential of a body at `(px, py, pz)` (C `gravsub_i_plasma`): softened monopole,
Ewald correction table read with its field Jacobian, then the cell dipole (and
quadrupole) terms. The `usequad` branches are compile-time (`Octree{Q}`).
"""
@inline function gravsub_plasma!(ax::Float64, ay::Float64, az::Float64, phi0::Float64,
                                 t::Octree{Q}, q::Int, px::Float64, py::Float64,
                                 pz::Float64, tab::EwaldTable, L::Float64,
                                 eps_tree::Float64) where {Q}
    # the C `gravsub_i_plasma` takes dr = pos0 − Pos(q) (opposite of `subdivp_s`)
    dx = bicv1(px - t.pos[1, q], L)
    dy = bicv1(py - t.pos[2, q], L)
    dz = bicv1(pz - t.pos[3, q], L)
    drsq = dx * dx + dy * dy + dz * dz + eps_tree * eps_tree
    drab = sqrt(drsq)            # C `rsqrt` is a plain square root (misnomer)
    charge = t.charge[q]
    phii = charge / drab
    mor3 = phii / drsq
    ax += mor3 * dx
    ay += mor3 * dy
    az += mor3 * dz
    phi0 += phii
    ew = ewald_force_jacobian_packed(tab, dx, dy, dz)
    f1, f2, f3 = ew.f
    qf1 = charge * f1
    qf2 = charge * f2
    qf3 = charge * f3
    ax += qf1
    ay += qf2
    az += qf3
    phi0 += charge * ew.pot
    if Q > 0 && t.ntype[q] == CELL_T
        dip1 = t.dip[1, q]
        dip2 = t.dip[2, q]
        dip3 = t.dip[3, q]
        dr3inv = 1.0 / (drsq * drab)
        phidip = (dip1 * dx + dip2 * dy + dip3 * dz) * dr3inv
        phi0 += phidip
        t3 = 3.0 * drab * phidip
        ax += (t3 * dx - dip1) * dr3inv
        ay += (t3 * dy - dip2) * dr3inv
        az += (t3 * dz - dip3) * dr3inv
        phi0 += dip1 * qf1 + dip2 * qf2 + dip3 * qf3
    end
    if Q > 1 && t.ntype[q] == CELL_T
        dip1 = t.dip[1, q]
        dip2 = t.dip[2, q]
        dip3 = t.dip[3, q]
        j1, j2, j3, j4, j5, j6, j7, j8, j9 = ew.j
        ax += dip1 * j1 + dip2 * j2 + dip3 * j3
        ay += dip1 * j4 + dip2 * j5 + dip3 * j6
        az += dip1 * j7 + dip2 * j8 + dip3 * j9
        dr5inv = 1.0 / (drsq * drsq * drab)
        qx = t.quad[1, q] * dx + t.quad[2, q] * dy + t.quad[3, q] * dz
        qy = t.quad[4, q] * dx + t.quad[5, q] * dy + t.quad[6, q] * dz
        qz = t.quad[7, q] * dx + t.quad[8, q] * dy + t.quad[9, q] * dz
        phiquad = -0.5 * dr5inv * (dx * qx + dy * qy + dz * qz)
        phi0 -= phiquad
        phiquad = 5.0 * phiquad / drsq
        ax -= phiquad * dx
        ay -= phiquad * dy
        az -= phiquad * dz
        ax -= dr5inv * qx
        ay -= dr5inv * qy
        az -= dr5inv * qz
    end
    return ax, ay, az, phi0
end

"""
    tree_force_i(tree, ip, list, tab, L, eps_tree)

Tree acceleration and potential of body `ip` from its frozen interaction list
(C `hackgrav_i_plasma`); the force on an electron is `-acc`.
"""
function tree_force_i(t::Octree, ip::Int, list::AbstractVector{Int}, tab::EwaldTable,
                      L::Float64, eps_tree::Float64)
    px = t.pos[1, ip]
    py = t.pos[2, ip]
    pz = t.pos[3, ip]
    ax = ay = az = 0.0
    phi0 = 0.0
    for q in list
        ax, ay, az, phi0 = gravsub_plasma!(ax, ay, az, phi0, t, q, px, py, pz, tab, L,
                                           eps_tree)
    end
    return Vec3(ax, ay, az), phi0
end
