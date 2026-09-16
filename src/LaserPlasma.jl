"""
    LaserPlasma

Idiomatic Julia port of the C hot-plasma / laser simulation (RWTH Aachen,
1999-2000). Atomic units, Ewald summation on a table, relativistic dynamics in
momenta, and the fixed-step (`verlet`) and adaptive individual-step
(`multistep_rec`) integrators of the oracle.
"""
module LaserPlasma

using Printf
using Base.Cartesian: @ntuple

include("constants.jl")
include("vec3.jl")
include("rng.jl")
include("sobol.jl")
include("plasma.jl")
include("ewald.jl")
include("force_types.jl")
include("forces.jl")
include("tree.jl")
include("laser.jl")
include("schemes.jl")
include("integrators.jl")
include("multistep.jl")
include("diagnostics.jl")
include("simulation.jl")
include("analysis.jl")
include("special.jl")
include("dielectric.jl")

end # module
