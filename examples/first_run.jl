# A short, inexpensive calculation with the same physical model as the slides.
# Run from the project root: julia --project=. examples/first_run.jl

using LaserPlasma

const PC = LaserPlasma

function main()
    p = PC.PlasmaParams(
        24,      # fixed ions (and 24 electrons because Z = 1)
        1,       # ion charge Z
        1.442,   # Wigner-Seitz radius
        6.935,   # kT
        0.144,   # electron-ion softening
        0.144,   # tree softening (unused with direct forces)
        1.57,    # laser amplitude F0
        3.0,     # laser pulsation omega
        10,      # number of laser cycles
        2,       # relaxation cycles
        2,       # ramp-up cycles
        true,    # relativistic electron dynamics
    )
    outdir = joinpath(@__DIR__, "..", "runs", "first_run")
    cfg = PC.RunConfig(12.0, 8, outdir)
    PC.run_multistep(p, cfg, 0.01; sample = 2, thermostat = true)
    println("Results: ", outdir)
end

main()
