#!/usr/bin/env julia
# Set up the pushdown validation project: dev-track the pushdown worktrees of
# EarthSciAST + EarthSciIO, add the other deps run-model.jl uses (minus Plots /
# OrdinaryDiffEqTsit5, which the prepare/_observed_field path never touches), plus
# GeometryOps + GeoInterface to activate the STRtree broad-phase fast path for the
# overlap join-gate.  Then instantiate + precompile.
import Pkg
Pkg.activate(@__DIR__)

# EA_PATH / IO_PATH default to the sibling checkouts; override via the env.
include(joinpath(@__DIR__, "paths.jl"))

for (label, var, path) in (("EarthSciAST", "EA_PATH", EA_PATH), ("EarthSciIO", "IO_PATH", IO_PATH))
    isdir(path) || error("$label checkout not found at $path — set the $var environment variable")
end
println("dev-tracking:\n  EarthSciAST: $EA_PATH\n  EarthSciIO:  $IO_PATH")

Pkg.develop([Pkg.PackageSpec(path=EA_PATH), Pkg.PackageSpec(path=IO_PATH)])
Pkg.add([
    Pkg.PackageSpec(name="Blosc"),
    Pkg.PackageSpec(name="JSON"),
    Pkg.PackageSpec(name="JSON3"),
    Pkg.PackageSpec(name="ZipFile"),
    Pkg.PackageSpec(name="NCDatasets"),
    Pkg.PackageSpec(name="GeometryOps"),
    Pkg.PackageSpec(name="GeoInterface"),
])
Pkg.instantiate()
Pkg.precompile()
println("SETUP OK")
