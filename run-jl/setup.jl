#!/usr/bin/env julia
# Set up the Phase-1 shim project: dev-track the EarthSciAST + EarthSciIO
# checkouts (env-overridable via EA_PATH / IO_PATH), add the runtime deps the
# prepare/observed-graph path touches (Blosc activates EarthSciIO's zarr codec
# ext; ZipFile for the TEMPORARY EGU ingest; GeometryOps + GeoInterface activate
# the STRtree broad-phase fast path for the overlap join-gate), then
# instantiate + precompile.
import Pkg
Pkg.activate(@__DIR__)

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
    Pkg.PackageSpec(name="GeometryOps"),
    Pkg.PackageSpec(name="GeoInterface"),
])
Pkg.instantiate()
Pkg.precompile()
println("SETUP OK")
