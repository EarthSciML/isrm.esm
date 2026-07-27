#!/usr/bin/env julia
# Set up the run-model project: dev-track the sibling EarthSciAST + EarthSciIO
# checkouts (the registered releases lack the zarr + ff10 readers this runner is
# built on), then instantiate + precompile.
#
# The Manifest.toml is deliberately gitignored — it pins machine-specific dev
# paths — so this script is how a fresh checkout becomes runnable.
#
#   EA_PATH   EarthSciAST.jl checkout   (default: ../../EarthSciAST/pkg/EarthSciAST.jl)
#   IO_PATH   EarthSciIO julia checkout (default: ../../EarthSciIO/julia)
import Pkg
Pkg.activate(@__DIR__)

const HERE      = @__DIR__
const REPO      = dirname(HERE)
const CODE_ROOT = dirname(REPO)

const EA_PATH = get(ENV, "EA_PATH", joinpath(CODE_ROOT, "EarthSciAST", "pkg", "EarthSciAST.jl"))
const IO_PATH = get(ENV, "IO_PATH", joinpath(CODE_ROOT, "EarthSciIO", "julia"))

for (label, var, path) in (("EarthSciAST", "EA_PATH", EA_PATH), ("EarthSciIO", "IO_PATH", IO_PATH))
    isdir(path) || error("$label checkout not found at $path — set the $var environment variable")
end
println("dev-tracking:\n  EarthSciAST: $EA_PATH\n  EarthSciIO:  $IO_PATH")

Pkg.develop([Pkg.PackageSpec(path=EA_PATH), Pkg.PackageSpec(path=IO_PATH)])
Pkg.add([
    Pkg.PackageSpec(name="Blosc"),          # blosc-compressed SR chunk decode
    Pkg.PackageSpec(name="JSON"),
    Pkg.PackageSpec(name="ZipFile"),        # the EGU FF10 zip
    Pkg.PackageSpec(name="NCDatasets"),
    Pkg.PackageSpec(name="OrdinaryDiffEqTsit5"),
    Pkg.PackageSpec(name="Plots"),          # plot.jl only
])
Pkg.instantiate()
Pkg.precompile()
println("SETUP OK")
