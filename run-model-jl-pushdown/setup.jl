#!/usr/bin/env julia
# Set up the pushdown validation project: dev-track the pushdown worktrees of
# EarthSciAST + EarthSciIO, add the other deps run-model.jl uses (minus Plots /
# OrdinaryDiffEqTsit5, which the prepare/_observed_field path never touches), plus
# GeometryOps + GeoInterface to activate the STRtree broad-phase fast path for the
# overlap join-gate.  Then instantiate + precompile.
import Pkg
Pkg.activate(@__DIR__)

const EA_PATH = "/Users/ctessum/code/earthsciml/EarthSciAST-wt-pushdown/pkg/EarthSciAST.jl"
const IO_PATH = "/Users/ctessum/code/earthsciml/EarthSciIO-wt-pushdown/julia"

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
