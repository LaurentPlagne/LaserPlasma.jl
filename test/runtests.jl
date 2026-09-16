using LaserPlasma
using Test

const PC = LaserPlasma

# Frozen reference values produced by the recompiled 2000 C oracle
# (eelinux04/C++/new). See AGENTS.md ("Chemins numériques") and
# docs/reprise/2026-09-13.md.

@testset "RNG identical to the oracle" begin
    r = PC.Ran2(-1)
    @test [PC.dran2!(r) for _ in 1:8] == [
        0.28538089909468611,
        0.25335818926591708,
        0.093468531009194042,
        0.60849689073964752,
        0.90342026007860998,
        0.19587319281381618,
        0.46295354298830549,
        0.93902133024149248,
    ]
    gd = PC.Gasdev()
    @test [PC.gasdev!(gd, r) for _ in 1:8] == [
        -0.22811536231958363,
        -1.0115252677965521,
        -0.97264016936075881,
        0.083595842060255054,
        1.678086777821336,
        -1.119635269354307,
        1.0152410766438908,
        1.9075967347134744,
    ]
    so = PC.Sobol()
    @test [Float64(x) for _ in 1:5 for x in PC.sobol_next!(so)] == [
        0.5, 0.5, 0.5,
        0.25, 0.75, 0.25,
        0.75, 0.25, 0.75,
        0.375, 0.625, 0.125,
        0.875, 0.125, 0.625,
    ]
    # consumption identical to init_plasma (sample=2, 3 electrons)
    r2 = PC.Ran2(-1)
    gd2 = PC.Gasdev()
    for _ in 1:3
        foreach(_ -> PC.dran2!(r2), 1:3)
        foreach(_ -> PC.gasdev!(gd2, r2), 1:3)
        foreach(_ -> PC.dran2!(r2), 1:3)
    end
    for _ in 1:3
        foreach(_ -> PC.dran2!(r2), 1:3)
        foreach(_ -> PC.gasdev!(gd2, r2), 1:3)
    end
    @test r2.idum == 1389597872
end

# Validation case: N=50, Γ=0.1 (kT=6.935), ε=0.144, F0=1.57.
const PC_PARAMS = PC.PlasmaParams(50, 1, 1.442, 6.935, 0.144, 0.144, 1.57, 3.0, 10, 2, 2,
                                  true)

@testset "Ewald table" begin
    tab = PC.build_ewald_table(PC.box_length(PC_PARAMS))
    @test tab.pot[1, 1, 1] ≈ -0.13677093545597976 atol = 1e-14
    @test tab.fx[1, 1, 1] ≈ 0.010497241274580248 atol = 1e-14
    @test tab.pot[1, 1, 2] ≈ -0.14142398592990479 atol = 1e-14
    @test tab.fz[1, 1, 2] ≈ 0.0068894781828234726 atol = 1e-14
    @test tab.pot[9, 9, 9] ≈ -0.23960980284027755 atol = 1e-14
    @test tab.pot[17, 17, 17] ≈ -0.13677093545597976 atol = 1e-14
end

@testset "Ewald table with field Jacobian (tree)" begin
    # Frozen values from the recompiled oracle (`calc_ewald_sum` after
    # build_tab_ewald(16)); see docs/architecture/barnes-hut.md, phase 2.
    # Each row is (x, y, z, pot, fx, fy, fz, fxx..fzz).
    ref13 = [
        -4.2817516057232599 -3.7465326550078526 -3.2113137042924449 -0.1502774780974985 0.015354922182881791 0.0094741793174999284 0.0047950790385521986 0.0068490672447747833 -0.0040434228765116112 -0.0034657910370099521 -0.0040434228765116112 0.0071159934874873432 -0.0021762626286082449 -0.0034657910370099521 -0.0021762626286082449 0.006045325405197634
        -2.71828 3.14159 -0.57721 -0.18690162472710353 0.01209425020892824 -0.01878067751178112 0.00073606497927815122 0.0083027476148237033 0.0065604268396436518 -0.00058076685230953786 0.0065604268396436509 0.010214058490660238 0.0010587335233004011 -0.00058076685230953786 0.0010587335233004011 0.0014935800319758271
        2.0000001 -2.0000001 4.0 -0.16874313987742087 -0.0014884571644713882 0.0014884571644713847 -0.028270269614131257 0.0040819250394283678 -0.00020392259874155988 -0.0062523647396144884 -0.00020392259874155988 0.004081925039428367 0.0062523647396144893 -0.0062523647396144884 0.0062523647396144893 0.011846536058603025
    ]
    tab = PC.build_ewald_hessian_table(8.56350321144652)
    for row in eachrow(ref13)
        q = PC.ewald_force_jacobian(tab, row[1], row[2], row[3])
        @test [q.pot; q.f...; q.j...] ≈ row[4:16] atol = 1e-14
    end
    # the Jacobian is symmetric (fxy ≈ fyx, ...), and the monopole readers of the
    # same table agree to the summation-order roundoff of `calc_ewald`
    q = PC.ewald_force_jacobian(tab, -2.71828, 3.14159, -0.57721)
    @test q.j[2] ≈ q.j[4] atol = 1e-15
    @test q.j[3] ≈ q.j[7] atol = 1e-15
    @test q.j[6] ≈ q.j[8] atol = 1e-15
    @test q.pot ≈ PC.ewald_pot(tab, -2.71828, 3.14159, -0.57721) atol = 1e-15
    @test collect(q.f) ≈ collect(PC.ewald_force(tab, -2.71828, 3.14159, -0.57721)) atol =
        1e-15
    # second box, grid node (3,6,9) and an interior point
    tab8 = PC.build_ewald_hessian_table(8.0)
    ref8_node = [-2.5, -1.0, 0.5, -0.22433909798906509, 0.024543314211333697,
                 0.0046805180030534549, -0.0020789008333079659, 0.014604589161439585,
                 -0.0032803512997566994, 0.0015912829454432301, -0.0032803512997566994,
                 0.0055591306976579546, 0.00025601440509775107, 0.0015912829454432301,
                 0.00025601440509775107, 0.0043799727470727259]
    q8 = PC.ewald_force_jacobian(tab8, ref8_node[1], ref8_node[2], ref8_node[3])
    @test [q8.pot; q8.f...; q8.j...] ≈ ref8_node[4:16] atol = 1e-14
    ref8_pt = [0.1234567, -0.7654321, 2.3456789, -0.22890692264088994,
               -0.00057770537228805984, 0.0038824998813650588, -0.023327539460903799,
               0.0046508210656266528, 5.4839184336724942e-5, -0.00037175229392393336,
               5.4839184336724942e-5, 0.0055759791009162406, 0.0023725057000087311,
               -0.00037175229392393336, 0.0023725057000087311, 0.014316892439627377]
    q8b = PC.ewald_force_jacobian(tab8, ref8_pt[1], ref8_pt[2], ref8_pt[3])
    @test [q8b.pot; q8b.f...; q8b.j...] ≈ ref8_pt[4:16] atol = 1e-14
    # explicit guard on a monopole-only table
    @test_throws ErrorException PC.ewald_force_jacobian(PC.build_ewald_table(8.0), 0.0,
                                                       0.0, 0.0)
    # the packed (vectorized) reader is bit-identical to the reference
    for row in eachrow(ref13)
        qp = PC.ewald_force_jacobian_packed(tab, row[1], row[2], row[3])
        qr = PC.ewald_force_jacobian(tab, row[1], row[2], row[3])
        @test [qp.pot; qp.f...; qp.j...] == [qr.pot; qr.f...; qr.j...]
    end
    @test_throws ErrorException PC.ewald_force_jacobian_packed(PC.build_ewald_table(8.0),
                                                               0.0, 0.0, 0.0)
end

@testset "Barnes-Hut octree (oracle dump)" begin
    # Frozen dump of the recompiled oracle `maketree` with 10 electrons and 10
    # ions (test/data/tree_ref.txt); see docs/architecture/barnes-hut.md, phase 3.
    rows = [split(l) for l in eachline(joinpath(@__DIR__, "data", "tree_ref.txt"))
            if !startswith(l, "#")]
    head = rows[1]
    rsize = parse(Float64, head[3])
    n_elec = parse(Int, head[4])
    n_ion = parse(Int, head[5])
    maxlevel = parse(Int, head[6])
    cellused = parse(Int, head[7])
    body = rows[2:1+n_elec+n_ion]
    noderows = rows[2+n_elec+n_ion:end]
    nbody = n_elec + n_ion
    mass = [parse(Float64, r[2]) for r in body]
    charge = [parse(Float64, r[3]) for r in body]
    pos = zeros(3, nbody)
    for (i, r) in enumerate(body)
        pos[:, i] .= parse.(Float64, r[4:6])
    end

    tree = PC.Octree(nbody; rsize = rsize, theta = 0.5, sw93 = true, usequad = 2)
    PC.build_tree!(tree, mass, charge, pos)
    @test tree.maxlevel == maxlevel
    @test tree.cellused == cellused
    @test tree.nbody + tree.cellused == nbody + cellused
    @test tree.mass[tree.root] == Float64(nbody)
    @test tree.charge[tree.root] == 0.0   # 10 electrons, 10 ions

    # pre-order traversal, children in subcell order, as the dump
    order = Int[]
    depths = Int[]
    function visit(id, d)
        push!(order, id)
        push!(depths, d)
        if tree.ntype[id] == PC.CELL_T
            for k in 1:8
                q = tree.subp[k, id]
                q == 0 || visit(q, d + 1)
            end
        end
    end
    visit(tree.root, 0)
    @test length(order) == length(noderows)
    index_of = Dict(order[i] => i for i in eachindex(order))
    @test depths == parse.(Int, getindex.(noderows, 2))
    @test Int[tree.ntype[id] for id in order] == parse.(Int, getindex.(noderows, 3))
    @test collect(0:length(order)-1) == parse.(Int, getindex.(noderows, 4))
    @test [(tree.more[id] == 0 ? -1 : index_of[tree.more[id]] - 1)
           for id in order] == parse.(Int, getindex.(noderows, 5))
    @test [(tree.next[id] == 0 ? -1 : index_of[tree.next[id]] - 1)
           for id in order] == parse.(Int, getindex.(noderows, 6))

    # positions, masses, charges (all nodes), then cell radii and moments
    refnum = permutedims(reduce(hcat, [parse.(Float64, r[7:24]) for r in noderows]))
    got = zeros(length(order), 18)
    for (i, id) in enumerate(order)
        got[i, 1:3] .= tree.pos[:, id]
        got[i, 4] = tree.mass[id]
        got[i, 5] = tree.charge[id]
        if tree.ntype[id] == PC.CELL_T
            got[i, 6] = tree.rcrit2[id]
            got[i, 7:9] .= tree.dip[:, id]
            got[i, 10:18] .= tree.quad[:, id]
        end
    end
    @test got[:, 1:5] ≈ refnum[:, 1:5] atol = 1e-13
    cells = findall(==(2), parse.(Int, getindex.(noderows, 3)))
    @test got[cells, 6:18] ≈ refnum[cells, 6:18] atol = 1e-13

    # rebuilding is idempotent and reuses the grown root cell
    rsize0 = tree.rsize
    PC.maketree!(tree)
    @test tree.rsize == rsize0
    @test tree.cellused == cellused
    @test tree.maxlevel == maxlevel
end

@testset "Barnes-Hut walk, lists and forces (oracle dump)" begin
    # Frozen dumps of init_plasma(sample=2) -> maketree -> make_interaction_list
    # -> hackgrav_i_plasma (see docs/architecture/barnes-hut.md, phase 4). The
    # theta=2 dump accepts cells with non-zero dipoles, exercising the Jacobian.
    function check_interl(path, theta)
        rows = [split(l) for l in eachline(path) if !startswith(l, "#")]
        head = [r for r in rows if r[1] == "state"][1]
        L = parse(Float64, head[2])
        rsize = parse(Float64, head[3])
        n_elec = parse(Int, head[4])
        n_ion = parse(Int, head[5])
        bodies = [r for r in rows if r[1] == "body"]
        pos = zeros(3, n_elec + n_ion)
        mass = [parse(Float64, b[2]) for b in bodies]
        charge = [parse(Float64, b[3]) for b in bodies]
        for (i, b) in enumerate(bodies)
            pos[:, i] .= parse.(Float64, b[4:6])
        end
        # the oracle initial state is reproduced bit-for-bit (sample = 2)
        p = PC.PlasmaParams(n_ion, 1, 1.442, 6.935, 0.05, 0.05, 79.0, 3.0, 10, 2, 2,
                            true)
        st = PC.PlasmaState(p)
        PC.init_plasma!(st, p; sample = 2)
        @test st.pos == pos[:, 1:n_elec]
        @test st.ions == pos[:, n_elec+1:end]

        tree = PC.Octree(n_elec + n_ion; rsize = rsize, theta = theta, sw93 = true,
                         usequad = 2)
        PC.build_tree!(tree, mass, charge, pos)
        order = Int[]
        function visit(id)
            push!(order, id)
            for k in 1:8
                q = tree.subp[k, id]
                q == 0 || visit(q)
            end
        end
        visit(tree.root)
        index_of = Dict(order[i] => i - 1 for i in eachindex(order))
        tab = PC.build_ewald_hessian_table(L)
        for row in [r for r in rows if r[1] == "list"]
            l = parse(Int, row[2]) + 1
            got = [index_of[q] for q in PC.build_interaction_list!(tree, l, L)]
            @test got == parse.(Int, row[4:end])
        end
        for row in [r for r in rows if r[1] == "force"]
            l = parse(Int, row[2]) + 1
            list = PC.build_interaction_list!(tree, l, L)
            acc, phi = PC.tree_force_i(tree, l, list, tab, L, 0.05)
            @test [acc[1], acc[2], acc[3], phi] ≈ parse.(Float64, row[3:6]) atol = 1e-13
        end
    end
    check_interl(joinpath(@__DIR__, "data", "tree_interl_ref.txt"), 0.5)
    check_interl(joinpath(@__DIR__, "data", "tree_interl_coarse_ref.txt"), 2.0)

    # the tree force is close to the direct one at theta = 0.5 (approximation)
    p = PC.PlasmaParams(10, 1, 1.442, 6.935, 0.05, 0.05, 79.0, 3.0, 10, 2, 2, true)
    st = PC.PlasmaState(p)
    PC.init_plasma!(st, p; sample = 2)
    L = PC.box_length(p)
    mass = vcat(ones(10), ones(10))
    charge = vcat(fill(-1.0, 10), ones(10))
    pos = hcat(st.pos, st.ions)
    tree = PC.Octree(20; rsize = L, theta = 0.5, sw93 = true, usequad = 2)
    PC.build_tree!(tree, mass, charge, pos)
    tab = PC.build_ewald_hessian_table(L)
    fdir = zeros(3, 10)
    PC.background_force!(fdir, st, PC.Forces(p), tab, L)
    PC.ee_force!(fdir, st, PC.Forces(p), tab, L)
    for l in 1:10
        acc, _ = PC.tree_force_i(tree, l, PC.build_interaction_list!(tree, l, L), tab, L,
                                 0.05)
        @test collect(acc.data) ≈ -fdir[:, l] atol = 0.05
    end
end

@testset "Fixed-step dynamics (verlet)" begin
    outdir = mktempdir()
    PC.run_verlet(PC_PARAMS, PC.RunConfig(1.0, 8, outdir))
    energy = reduce(vcat,
                    [permutedims(parse.(Float64, split(l)))
                     for l in eachline(joinpath(outdir, "energy.dat"))])
    ref_energy = [
        0.1308996938995747 491.26045093638646 -115.6090700831156 375.65138085327089
        0.91629785729702284 491.25099266093127 -115.71748017956102 375.53351248137028
    ]
    @test energy[[1, 7], :] ≈ ref_energy atol = 1e-11
    pos = reduce(vcat,
                 [permutedims(parse.(Float64, split(l)))
                  for l in eachline(joinpath(outdir, "pos_0001"))])
    ref_pos = [
        0.91629785729702284 2.3061283395536165 2.629710558928438 1.1080529662171683
        0.91629785729702284 -0.63257163971790231 3.5408347028002449 -0.34193235952491752
        0.91629785729702284 -0.27562132352363816 2.3866701604530576 1.4848796385108336
    ]
    @test pos[[1, 2, 50], :] ≈ ref_pos atol = 1e-12
    mom = reduce(vcat,
                 [permutedims(parse.(Float64, split(l)))
                  for l in eachline(joinpath(outdir, "mome_0001"))])
    ref_mom = [
        0.91629785729702284 0.39772882552066413 -1.4716498084714909 3.3685941920386178
        0.91629785729702284 1.9770359305700966 0.75826904087277058 0.48852123722950869
        0.91629785729702284 -5.1566291797481991 -1.5690953469608822 4.5544704043583195
    ]
    @test mom[[1, 2, 50], :] ≈ ref_mom atol = 1e-12
end

@testset "Adaptive individual step (multistep)" begin
    function multistep_energies(eps; nsteps = 6, interp = :quad)
        p = PC_PARAMS
        L = PC.box_length(p)
        tab = PC.build_ewald_table(L)
        fm = PC.Forces(p)
        st = PC.PlasmaState(p)
        PC.init_plasma!(st, p; sample = 2)
        dx = PC.laser_half_period(p) / 8
        es = [PC.energies(st, p, fm, tab, L)[3]]
        t = 0.0
        for _ in 1:nsteps
            PC.multistep_step!(st, p, fm, tab, L, t, dx, eps; interp = interp)
            t += dx
            push!(es, PC.energies(st, p, fm, tab, L)[3])
        end
        return es, st
    end
    # eps=0.5: no recursion, pure RK4 step doubling
    es05, _ = multistep_energies(0.5)
    ref05 = [381.62265712869691, 381.48884998600124, 381.49364270033635,
             381.42474472454671, 381.43302969839408, 381.40100669481041,
             381.44413984013551]
    @test es05 ≈ ref05 atol = 1e-11
    # eps=0.01: recursive subdivision (populations 4, 2, 0 at levels 1..3)
    es01, st01 = multistep_energies(0.01)
    ref01 = [381.62265712869691, 381.61224807772652, 381.62275050330618,
             381.61537981323715, 381.62856813873225, 381.58693934566247,
             381.62928513349141]
    @test es01 ≈ ref01 atol = 1e-11
    refpos = [2.2725012953161121 -0.86546373369805141
              2.824133209463632 3.39593853247337
              0.59591741097249507 -0.36556681726472062]
    @test [st01.pos[1, 1] st01.pos[1, 2]
           st01.pos[2, 1] st01.pos[2, 2]
           st01.pos[3, 1] st01.pos[3, 2]] ≈ refpos atol = 1e-12
    # interpolation stencil of the frozen particles: the lower-order variants must
    # run and change the trajectory, and the quadratic default must be untouched
    esq, _ = multistep_energies(0.5; interp = :quad)
    @test esq == es05
    esp, _ = multistep_energies(0.01; interp = :pwise)
    @test all(isfinite, esp) && any(esp .!= es01)
    ess, _ = multistep_energies(0.01; interp = :simple)
    @test all(isfinite, ess) && any(ess .!= es01)
    # same through the full run loop (`energy.dat` holds start-of-step states)
    outdir = mktempdir()
    PC.run_multistep(PC_PARAMS, PC.RunConfig(1.0, 8, outdir), 0.01)
    energy = reduce(vcat,
                    [permutedims(parse.(Float64, split(l)))
                     for l in eachline(joinpath(outdir, "energy.dat"))])
    @test energy[:, 4] ≈ ref01 atol = 1e-11
end

@testset "Frozen stencil: closure under dyadic refinement" begin
    # A particle frozen at some level is handed to its children as the triple
    # (deb, value at the quarter point, value at the half point); the child
    # re-interpolates from that triple. The quadratic stencil is a fixed point of
    # this nesting -- the three nodes it passes down all lie on its own parabola --
    # so a frozen trajectory is ONE polynomial over the whole interval whatever the
    # depth. The linear stencils re-linearize at every level and compound their
    # error: this is the mechanism behind the accuracy floor measured in
    # docs/notes/interp-order.md §4.2, so it is worth pinning.
    n = 1
    mk(v) = (s = PC.particle_set(n); s.pos[1, 1] = v; s.mom[1, 1] = v; s)
    stable = trues(n)
    # nodes of the freezing level, taken on a parabola through t = 0, 1/2, 1
    P(t) = 0.7 + 1.3t - 2.1t^2
    deb, mid, fin = mk(P(0.0)), mk(P(0.5)), mk(P(1.0))
    for (interp, exact) in ((:quad, true), (:pwise, false), (:simple, false))
        out = PC.particle_set(n)
        # child covering [0, 1/2]: its known triple is (deb, quarter, half)
        PC.interp_rec_fh!(stable, deb, mid, fin, out; interp = interp)
        cdeb, cmid, cfin = deb, out, mid
        nested = PC.particle_set(n)
        PC.interp_rec_fh!(stable, cdeb, cmid, cfin, nested; interp = interp)
        # the child's quarter point is t = 1/8 of the parent interval
        d = abs(nested.pos[1, 1] - P(0.125))
        exact ? (@test d < 1e-15) : (@test d > 1e-3)
    end
end

@testset "Recursion profile (oracle levels)" begin
    # Frozen from the recompiled C printout (`level = k ; population = n`) for the
    # first step of r14i0500e0500q001a00/F0157W0300N010 (N=500, eps=0.005). The
    # unstable population must collapse with the level: this is the cost argument
    # of the method, and it validates the stability masks, not only trajectories.
    p = PC.PlasmaParams(500, 1, 1.442, 6.935, 0.05, 0.05, 1.57, 3.0, 10, 2, 2, true)
    L = PC.box_length(p)
    tab = PC.build_ewald_table(L)
    fm = PC.Forces(p)
    st = PC.PlasmaState(p)
    PC.init_plasma!(st, p; sample = 2)
    w = PC.MultistepWork(500, 15)
    PC.multistep_step!(st, p, fm, tab, L, 0.0, PC.laser_half_period(p) / 8, 0.005;
                       work = w)
    @test w.visits[1:6] == [1, 2, 4, 8, 14, 2]
    @test w.populations[1:6] == [98, 58, 27, 12, 1, 0]
end

@testset "Adaptive multistep with Barnes-Hut tree (oracle)" begin
    # Frozen values from the recompiled oracle with use_tree=1 (ppbs_min run):
    # N=50, Γ=0.1, ε=0.144, F0=1.57, ω=3, eps=0.01, theta=0.5, 7 steps.
    p = PC.PlasmaParams(50, 1, 1.442, 6.935, 0.144, 0.144, 1.57, 3.0, 10, 2, 2, true)
    outdir = mktempdir()
    PC.run_multistep(p, PC.RunConfig(1.0, 8, outdir), 0.01; use_tree = true)
    readmat(path) = reduce(vcat,
                           [permutedims(parse.(Float64, split(l)))
                            for l in eachline(path) if !isempty(strip(l))])
    ref_energy = [
        0.0 495.45440778206813 -113.8317506533712 381.62265712869691
        0.1308996938995747 496.78326649987599 -115.16728136417235 381.61598513570362
        0.26179938779914941 496.04010247067038 -114.42043218909538 381.619670281575
        0.39269908169872414 500.27616816799275 -118.64322053382867 381.63294763416411
        0.52359877559829882 500.09617296073446 -118.44282736891681 381.65334559181764
        0.65449846949787349 496.48437200318199 -114.84715486713851 381.63721713604349
        0.78539816339744817 495.3199991189241 -113.61697961843734 381.70301950048679
    ]
    energy = readmat(joinpath(outdir, "energy.dat"))
    @test energy ≈ ref_energy atol = 1e-11
    ref_pos = [
        0.78539816339744817 2.2718001670803072 2.8238563303838697 0.5957238849981501
        0.78539816339744817 -0.86509700870275907 3.3952224473237744 -0.3655747531307702
        0.78539816339744817 0.36294655568566586 2.5142826113109402 1.0142855225178344
    ]
    pos = readmat(joinpath(outdir, "pos_0001"))
    @test pos[[1, 2, 50], :] ≈ ref_pos atol = 1e-12
    ref_mom = [
        0.78539816339744817 0.36829521399640508 -1.4453459373963184 3.1678994862977174
        0.78539816339744817 2.0614474089593791 0.71141997807789148 0.5481283079090139
        0.78539816339744817 -5.1631360396058446 -1.4208992071789646 4.7668632348052355
    ]
    mom = readmat(joinpath(outdir, "mome_0001"))
    @test mom[[1, 2, 50], :] ≈ ref_mom atol = 1e-12
end

@testset "Constant-temperature scheme (thermostat)" begin
    # archived configuration r14i0100e0100q001a01/F0078W0300N010/DIRE_RL_ELE:
    # N=100, Γ=0.1, ε=0.144, F0=0.785, ω=3, eps=0.005, thermostat on
    p = PC.PlasmaParams(100, 1, 1.442, 6.935, 0.144, 0.144, 0.785, 3.0, 10, 2, 2, true)
    outdir = mktempdir()
    PC.run_multistep(p, PC.RunConfig(10.0, 8, outdir), 0.005; thermostat = true)
    energy = reduce(vcat,
                    [permutedims(parse.(Float64, split(l)))
                     for l in eachline(joinpath(outdir, "energy.dat"))])
    ref_energy = [
        0.0 1205.3467842411364 -382.32067139911578 823.02611284202067
        9.8174770424681075 1067.2148262506512 -387.22818839542373 679.9866378552274
    ]
    @test energy[[1, 76], :] ≈ ref_energy atol = 1e-7
    du = reduce(vcat,
                [permutedims(parse.(Float64, split(l)))
                 for l in eachline(joinpath(outdir, "deltaU"))])
    @test du[1, :] ≈ [8.3775804095727811, 0.0, 683.43139323813648] atol = 1e-7
end

@testset "Checkpoint and restart" begin
    p = PC.PlasmaParams(50, 1, 1.442, 6.935, 0.144, 0.144, 1.57, 3.0, 10, 2, 2, true)
    out_cont = mktempdir()
    out_chunk = mktempdir()
    PC.run_multistep(p, PC.RunConfig(12.0, 8, out_cont), 0.01)
    PC.run_multistep(p, PC.RunConfig(6.0, 8, out_chunk), 0.01; checkpoint_every = 5)
    PC.run_multistep(p, PC.RunConfig(12.0, 8, out_chunk), 0.01; checkpoint_every = 5,
                     restart = true)
    readmat(path) = reduce(vcat,
                           [permutedims(parse.(Float64, split(l)))
                            for l in eachline(path) if !isempty(strip(l))])
    for f in ("energy.dat", "power.dat", "deltaU")
        @test readmat(joinpath(out_cont, f)) == readmat(joinpath(out_chunk, f))
    end
end

@testset "Special functions" begin
    @test PC.si(1.0) ≈ 0.946083070367183 atol = 1e-9
    @test PC.si(10.0) ≈ 1.658347594218874 atol = 1e-7
    @test PC.bessel_j0(2.5) ≈ -0.0483837764681977 atol = 1e-8
    @test PC.bessel_j1(2.5) ≈ 0.4970941024642741 atol = 1e-10
    @test PC.bessel_jn_safe(2, 1.0) ≈ 0.1149034849319005 rtol = 1e-6
    @test PC.bessel_jn_safe(3, 1.0) ≈ 0.0195633539826684 rtol = 1e-6
    @test PC.bessel_k1(1.0) ≈ 0.6019072301972346 rtol = 1e-6
end

@testset "Dielectric theory (decker.c)" begin
    vth = sqrt(6.935)
    # archived values, Diel.w10.xmgr (ω=10)
    th10 = PC.DielectricTheory(6.935, 10.0, 1.0, 0.05)
    @test PC.decker_coll_freq(0.01 * vth, th10) ≈ 0.0106671 rtol = 1e-3
    @test PC.decker_coll_freq(1.08251 * vth, th10) ≈ 0.0092457 rtol = 1e-3
    @test PC.decker_coll_freq(9.10282 * vth, th10) ≈ 0.000335768 rtol = 1e-3
    # soft core ε=0.05 (our converged value differs by 2.4% from the archived
    # one at the last point; see docs/reprise/2026-09-13.md)
    @test PC.decker_coll_freq(0.01 * vth, th10; softcore = true) ≈ 0.0206258 rtol = 1e-3
    @test PC.decker_coll_freq(1.08251 * vth, th10; softcore = true) ≈ 0.0170588 rtol = 1e-3
    @test PC.decker_coll_freq(100.0 * vth, th10; softcore = true) ≈ 7.66633e-7 rtol = 3e-2
    # soft core ε=0.144
    th144 = PC.DielectricTheory(6.935, 10.0, 1.0, 0.144)
    @test PC.decker_coll_freq(0.01 * vth, th144; softcore = true) ≈ 0.00856309 rtol = 1e-3
    @test PC.decker_coll_freq(1.08251 * vth, th144; softcore = true) ≈ 0.0073389 rtol = 1e-3
end

@testset "Collision frequency analysis (freq.c)" begin
    # constant ΔU per cycle: ν = 2 (ΔU/(N Δt))/v0²
    du = zeros(101, 3)
    du[:, 1] .= 0:100
    du[2:end, 2] .= 2.0
    du[:, 3] .= cumsum(du[:, 2])
    cf = PC.collision_frequency(du; nb_elec = 100, F0 = 1.0, omega = 1.0)
    @test cf.nu ≈ 0.04
    @test cf.sigma ≈ 0.0 atol = 1e-12
    @test cf.ncycles == 100
    @test cf.dt ≈ 1.0
end

@testset "Grace/xmgr reader" begin
    include(joinpath(@__DIR__, "..", "figures", "common.jl"))
    gs = read_xmgr(joinpath(@__DIR__, "data", "sample.xmgr"))
    @test length(gs) == 2
    g0 = gs[1]
    @test g0.name == "g0"
    @test g0.title == "sample"
    @test g0.xlabel == "x" && g0.ylabel == "y"
    @test g0.logx && g0.logy
    @test g0.world["xmax"] == 10.0
    @test [s.name for s in g0.sets] == ["S0", "S1"]
    @test size(g0.sets[1].data) == (3, 2)
    @test g0.sets[1].data[3, :] == [10.0, 1.0]
    @test size(g0.sets[2].data) == (2, 3)
    @test g0.sets[2].legend == "points"
    @test g0.sets[2].data[2, :] == [1.5, 0.15, 0.02]
    @test gs[2].name == "g1"
    @test [s.name for s in gs[2].sets] == ["S2"]
end

@testset "Type stability" begin
    @inferred PC.box_length(PC_PARAMS)
    @inferred PC.nb_elec(PC_PARAMS)
    @inferred PC.laser_electric_field(PC_PARAMS, 1.0)
    @inferred PC.build_ewald_table(8.0)
    @inferred PC.build_ewald_hessian_table(8.0)
    @inferred PC.ewald_force_jacobian(PC.build_ewald_hessian_table(8.0), 0.1, 0.2, 0.3)
    @inferred PC.ewald_force_jacobian_packed(PC.build_ewald_hessian_table(8.0), 0.1, 0.2,
                                             0.3)
    @inferred PC.build_tree!(PC.Octree(2; rsize = 8.0), [1.0, 1.0], [1.0, -1.0],
                             [0.1 -0.2; 0.3 0.4; -0.5 0.6])
    tree2 = PC.Octree(2; rsize = 8.0)
    PC.build_tree!(tree2, [1.0, 1.0], [1.0, -1.0], [0.1 -0.2; 0.3 0.4; -0.5 0.6])
    @inferred PC.build_interaction_list!(tree2, 1, 8.0)
    @inferred PC.tree_force_i(tree2, 1, PC.build_interaction_list!(tree2, 1, 8.0),
                              PC.build_ewald_hessian_table(8.0), 8.0, 0.05)
    @inferred PC.PlasmaState(PC_PARAMS)
    @inferred PC.Forces(PC_PARAMS)
    @inferred PC.energies(PC.PlasmaState(PC_PARAMS), PC_PARAMS, PC.Forces(PC_PARAMS),
                          PC.build_ewald_table(8.56350321144652), 8.56350321144652)
end

# The recursion is not symplectic and freezing a trajectory by interpolation breaks
# Newton's third law, so the energy budget has to be watched. Measured (session 5,
# `benchmarks/energy_drift.jl`): with DIRECT forces the drift is -4.8e-7 of E_kin per
# laser cycle at the production tolerance, independent of N, and it falls below the
# measurement floor at eps = 1e-5 for +0.5% of cost. This test is the guard: it is
# loose enough not to be brittle, tight enough that the regression which motivated it
# -- a force scheme that heats by 1e-2 of E_kin per cycle -- cannot slip through.
@testset "Adaptive recursion conserves energy (direct forces)" begin
    # a cold, dense plasma so that a handful of cycles is already representative
    p = PC.PlasmaParams(64, 1, 1.442, 6.935, 0.05, 0.05, 0.0, 3.0, 10, 2, 2, true)
    L = PC.box_length(p)
    tab = PC.build_ewald_table(L)
    fm = PC.Forces(p)
    st = PC.PlasmaState(p)
    PC.init_plasma!(st, p; sample = 2)
    w = PC.MultistepWork(64, 15)
    h = PC.laser_half_period(p) / 32
    ekin0 = 1.5 * 64 * p.kT
    e0 = PC.energies(st, p, fm, tab, L)[3]
    t = 0.0
    for _ in 1:(64 * 4)          # four laser cycles
        PC.multistep_step!(st, p, fm, tab, L, t, h, 1.748e-3; work = w,
                           scheme = PC.DirectScheme())
        t += h
        PC.wrap_in_box!(st, L)
    end
    drift = (PC.energies(st, p, fm, tab, L)[3] - e0) / ekin0
    @test abs(drift) < 1e-4      # measured ~1e-6 over four cycles
end

# The direct scheme reads the Ewald correction from the node-major `packed` layout
# (`pack_nodes`), which LLVM vectorizes; the `(i,j,k)` arrays are kept as the
# reference. The two must agree to the last bit -- each component accumulates its
# eight nodes in the same order -- and that is what makes the layout change free.
@testset "Ewald table: packed layout is bit-identical to the strided read" begin
    L = 8.5635032114465
    tab = PC.build_ewald_table(L)
    function strided(tab, x, y, z)
        h, invh = tab.h, tab.invh
        i = trunc(Int, (x - tab.xmin) * invh)
        j = trunc(Int, (y - tab.xmin) * invh)
        k = trunc(Int, (z - tab.xmin) * invh)
        bx = x - (tab.xmin + i * h); by = y - (tab.xmin + j * h); bz = z - (tab.xmin + k * h)
        ax = h - bx; ay = h - by; az = h - bz
        a = (ax * ay * az, bx * ay * az, ax * by * az, bx * by * az,
             ax * ay * bz, bx * ay * bz, ax * by * bz, bx * by * bz)
        ip, jp, kp = i + 1, j + 1, k + 1
        idx = ((ip, jp, kp), (ip + 1, jp, kp), (ip, jp + 1, kp), (ip + 1, jp + 1, kp),
               (ip, jp, kp + 1), (ip + 1, jp, kp + 1), (ip, jp + 1, kp + 1),
               (ip + 1, jp + 1, kp + 1))
        f(t) = a[1] * t[idx[1]...] + a[2] * t[idx[2]...] + a[3] * t[idx[3]...] +
               a[4] * t[idx[4]...] + a[5] * t[idx[5]...] + a[6] * t[idx[6]...] +
               a[7] * t[idx[7]...] + a[8] * t[idx[8]...]
        return tab.coef * f(tab.fx), tab.coef * f(tab.fy), tab.coef * f(tab.fz)
    end
    rng = PC.Ran2(-7)
    for _ in 1:200
        x = PC.min_image(L * (PC.dran2!(rng) - 0.5), L)
        y = PC.min_image(L * (PC.dran2!(rng) - 0.5), L)
        z = PC.min_image(L * (PC.dran2!(rng) - 0.5), L)
        @test PC.ewald_force(tab, x, y, z) === strided(tab, x, y, z)
    end
    # and the guard that turns an out-of-box read into an error rather than garbage
    @test_throws ArgumentError PC.ewald_force(tab, 0.6 * L, 0.0, 0.0)
end

# `PerParticleScheme` assembles the same physics as the oracle's symmetric pair sum
# but per particle, so it threads. It is deliberately NOT bit-identical -- Newton's
# third law changes the summation order -- and the two must agree to rounding.
@testset "Per-particle force assembly matches the symmetric pair sum" begin
    p = PC.PlasmaParams(200, 1, 1.442, 6.935, 0.05, 0.05, 1.57, 3.0, 10, 2, 2, true)
    L = PC.box_length(p)
    tab = PC.build_ewald_table(L)
    fm = PC.Forces(p)
    st = PC.PlasmaState(p)
    PC.init_plasma!(st, p; sample = 2)
    fa = zeros(3, PC.nb_elec(p))
    fb = similar(fa)
    t = 3.0                       # past the ramp, so the laser term is live
    PC.verlet_forces!(fa, st, p, fm, tab, L, t, PC.DirectScheme())
    PC.verlet_forces!(fb, st, p, fm, tab, L, t, PC.PerParticleScheme())
    @test maximum(abs.(fb .- fa) ./ max.(abs.(fa), eps())) < 1e-10
    # and it must be a real force, not a zeroed buffer that trivially "agrees"
    @test maximum(abs, fa) > 0
end
