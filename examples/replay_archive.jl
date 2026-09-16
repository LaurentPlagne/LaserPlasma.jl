# Optional long calculation using the physical settings of one historical run.
# Run: julia --project=. examples/replay_archive.jl 79 25.0
# The second argument is the final time; 368.614 requests the archived horizon.

using LaserPlasma

const PC = LaserPlasma

function main(args)
    length(args) in (1, 2) || error("Usage: replay_archive.jl 79|157 [final_time]")
    name = args[1]
    name in ("79", "157") || error("The first argument must be 79 or 157")
    final_time = length(args) == 2 ? parse(Float64, args[2]) : 25.0
    final_time >= 1 || error("final_time must be at least 1 (the C truncates it)")

    high_field = name == "79"
    softening = high_field ? 0.05 : 0.144
    amplitude = high_field ? 79.0 : 1.57
    nb_step = high_field ? 32 : 8
    tolerance = high_field ? 0.001748 : 0.01
    p = PC.PlasmaParams(2000, 1, 1.442, 6.935, softening, softening,
                        amplitude, 3.0, 10, 2, 2, true)
    outdir = joinpath(@__DIR__, "..", "runs", "replay_mesure" * name)
    cfg = PC.RunConfig(final_time, nb_step, outdir)
    PC.run_multistep(p, cfg, tolerance; sample = 2, thermostat = true,
                     use_tree = true, theta = 0.5, usequad = 2,
                     checkpoint_every = 2 * nb_step, restart = true)
    println("Results: ", outdir)
end

main(ARGS)
