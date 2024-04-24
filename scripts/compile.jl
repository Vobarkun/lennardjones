using Pkg
Pkg.activate("LennardJones")

using PackageCompiler
PackageCompiler.create_sysimage(["LennardJones"]; sysimage_path="lennardjones.so", precompile_execution_file="start.jl", precompile_statements_file="glmakie_precompile.jl")