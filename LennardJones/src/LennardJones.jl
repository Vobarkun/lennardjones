module LennardJones

using GLMakie, LinearAlgebra, StaticArrays, Distances, OffsetArrays, Dates, ForwardDiff, StatsBase, Distributions, ImageFiltering
const SVec2 = SVector{2, Float64};

include("sim.jl")
include("app.jl")

end # module LennardJones