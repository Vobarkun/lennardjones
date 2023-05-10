Base.@kwdef mutable struct LJ
    ps::Vector{SVector{2, Float64}} = []
    vs::Vector{SVector{2, Float64}} = []
    fs::Vector{SVector{2, Float64}} = []
    ms::Vector{Float64} = []
    cs::Vector{Float64} = []
    σs::Vector{Float64} = []
    
    bonds::Vector{Tuple{Int64, Int64, Float64, Float64}} = []
    angles::Vector{Tuple{Int64, Int64, Int64, Float64, Float64}} = []
    bonded::Matrix{Bool} = []

    nbs::Vector{Tuple{Int64, Int64}} = []
    nblist::Vector{Vector{Int64}} = []
    lastnbps::Vector{SVector{2, Float64}} = []
    
    σ::Float64 = 0.05
    ε::Float64 = 1.0
    swstrength::Float64 = 0.0
    swangle::Float64 = 120.0
    coulomb::Float64 = 1.0
    bondk::Float64 = 5.0
    bondl::Float64 = 1.0
    wallr::Float64 = 1.0
    wallk::Float64 = 1e6
    g::Float64 = 0.0
    E::Float64 = 0.0
    
    T::Float64 = 0.0
    ts::Float64 = 1.0
    
    nsteps::Int = 0
end

function randb(r)
    r * SA[2rand() - 1, 2rand() - 1]
end
randb(lj::LJ) = randb(lj.wallr)

function randb(r, n)
    [randb(r) for i in 1:n]
end
randb(lj::LJ, n) = randb(lj.wallr, n)

function LJ(N::Int; kwargs...)
    wallr = get(kwargs, :wallr, 1.0)
    σ = get(kwargs, :σ, 0.05)
    ps = randb(wallr, 1)
    for i in 1:10000
        length(ps) >= N && break
        p = randb(wallr)
        if all(norm(p - q) > σ for q in ps)
            push!(ps, p)
        end
    end
    LJ(; ps = ps, vs = zero(ps), fs = zero(ps), ms = ones(N), cs = zeros(N), σs = ones(N), bonded = fill(false, N, N), kwargs...)
end

function updatebonded!(lj::LJ)
    lj.bonded = fill(false, length(lj), length(lj))
    for (i, j, k, l) in lj.bonds
        lj.bonded[i,j] = lj.bonded[j, i] = true 
    end
end

function Base.push!(lj::LJ, p; m = 1.0, c = 0.0, σ = 1.0)
    N = length(lj)
    push!(lj.ps, p); push!(lj.vs, zero(p)); push!(lj.fs, zero(p))
    push!(lj.ms, m); push!(lj.cs, c); push!(lj.σs, σ)
    updatebonded!(lj)
end

function Base.push!(lj::LJ; ps, ms = nothing, cs = nothing, σs = nothing, bonds = [], angles = [])
    N = length(lj)
    append!(lj.ps, ps)
    append!(lj.vs, zero(ps))
    append!(lj.fs, zero(ps))
    append!(lj.ms, isnothing(ms) ? ones(length(ps)) : ms)
    append!(lj.cs, isnothing(cs) ? zeros(length(ps)) : cs)
    append!(lj.σs, isnothing(σs) ? ones(length(ps)) : σs)
    append!(lj.bonds, [(i + N, j + N, k, l) for (i, j, k, l) in bonds])
    append!(lj.angles, [(i + N, j + N, k + N, kα, α) for (i, j, k, kα, α) in angles])
    updatebonded!(lj)
    lj
end

function Base.deleteat!(lj::LJ, ind)
    deleteat!(lj.ps, ind); deleteat!(lj.vs, ind); deleteat!(lj.fs, ind)
    deleteat!(lj.ms, ind); deleteat!(lj.cs, ind); deleteat!(lj.σs, ind)
    lj.bonds = [(i - (i > ind), j - (j > ind), k, l) for (i, j, k, l) in lj.bonds if i != ind && j != ind]
    lj.angles = [(i - (i > ind), j - (j > ind), k - (k > ind), kα, α) for (i, j, k, kα, α) in lj.angles if i != ind && j != ind && k != ind]
    updatebonded!(lj)
    lj.nbs .= Ref((0, 0))
    lj.lastnbps .*= 0
    fill!.(lj.nblist, 0)
    lj
end

function Base.empty!(lj::LJ)
    empty!(lj.ps); empty!(lj.vs); empty!(lj.fs)
    empty!(lj.ms); empty!(lj.cs); empty!(lj.σs)
    empty!(lj.bonds); empty!(lj.angles); empty!(lj.lastnbps); empty!(lj.nbs); empty!(lj.nblist)
    updatebonded!(lj)
    lj
end

function Base.length(lj::LJ)
    length(lj.ps)
end

function pushfree!(lj, N = 1; ps = [SA[0.0,0.0]], kwargs...)
    ps .-= Ref(mean(ps))
    for i in 1:N
        qs = placefree(lj, ps)
        if !isnothing(qs)
            push!(lj; ps = qs, kwargs...)
        end
    end
end

function placefree(lj::LJ, ps; r = NaN, ntries = 100)
    r = ifelse(isnan(r), lj.σ, r)
    for i in 1:ntries
        p = randfree(lj, r = r)
        ϕ = 2pi * rand()
        qs = Ref(SA[cos(ϕ) sin(ϕ); -sin(ϕ) cos(ϕ)]) .* ps .+ Ref(p)
        if all(all(abs.(q) .< lj.wallr) for q in qs) && all(norm(p - q) > r for p in lj.ps for q in qs)
            return qs
        end
    end
    return nothing
end

function randfree(lj::LJ; r = NaN, ntries = 10000)
    length(lj) == 0 && return randb(lj)
    r = ifelse(isnan(r), lj.σ, r)
    dmin = Inf
    pmin = SA[0.0,0.0]
    for i in 1:ntries
        p = randb(lj)
        d = minimum(norm(p - q) for q in lj.ps)
        d = min(d, minimum(lj.wallr .- abs.(p)) + lj.σ)
        d > r && return p
        if d < dmin
            dmin = d
            pmin = p
        end
    end
    return pmin
end

function ljf(σ, ε, r; b = 1.0)
    rr = r ⋅ r
    σr = (σ^2 / rr) ^ 3
    (σr > 1/250) * min(1e10, 24ε / rr * (2σr^2 - b * σr)) * r 
end

function ljp(σ, ε, r; b = 1.0)
    rr = r ⋅ r
    σr = (σ^2 / rr) ^ 3
    4ε * (σr^2 - b * σr)
end

function angle(a, b)
    acosd(clamp(a ⋅ b / (norm(a) * norm(b)), -1, 1))
end

function smoothstep(x, x1, x2, y1 = 0, y2 = 1)
    x = clamp((x - x1) / (x2 - x1), 0, 1)
    s = x^2 * (3 - 2x)
    y1 * (1 - s) + y2 * s
end

function getcutoff(lj::LJ)
    ifelse(any(c != 0 for c in lj.cs), 10lj.σ, 5lj.σ)
end

function swh(x1, y1, x2, y2, a, θ)
    r1 = sqrt(x1^2 + y1^2); r2 = sqrt(x2^2 + y2^2)
    c = (x1 * x2 + y1 * y2) / (r1 * r2)
    20 * exp(0.5a / (r1 - a) + 0.5a / (r2 - a)) * (c - cosd(θ))^2
end

function swpotential_(x1, y1, x2, y2, x3, y3, a, θ)
    swh(x2 - x1, y2 - y1, x3 - x1, y3 - y1, a, θ)# + swh(x1 - x2, y1 - y2, x3 - x2, y3 - y2, a) + swh(x2 - x3, y2 - y3, x1 - x3, y1 - y3, a)
end

function swgrad(r1, r2, r3, a, θ)
    if norm(r2 - r1) >= a || norm(r3 - r1) >= a
        return SA[zero(r1), zero(r2), zero(r3)]
    end
    f = -ForwardDiff.gradient(x -> swpotential_(x..., a, θ), vcat(r1, r2, r3))
    SA[f[SA[1,2]], f[SA[3,4]], f[SA[5,6]]]
end

function swforce(r1, r2, r3, a, θ)
    f1 = swgrad(r1, r2, r3, a, θ)
    # f2 = swgrad(r2, r3, r1, a)
    # f3 = swgrad(r3, r1, r2, a)
    # f1[1] + f2[3] + f3[2], f1[2] + f2[1] + f3[3], f1[3] + f2[2] + f3[1]
end

function swpotential(r1, r2, r3, a, θ)
    if norm(r2 - r1) >= a || norm(r3 - r1) >= a
        return 0.0
    end
    swpotential_(r1..., r2..., r3..., a, θ)
end

function forces!(lj::LJ; cutoff = Inf)
    (; ps, vs, ms, fs, cs, nbs, bonds, angles, σ, σs, ε, coulomb, wallr, wallk, bondk, bondl, g, E, nblist, swstrength, swangle, bonded) = lj

    fill!(fs, zero(eltype(fs)))

    if swstrength != 0
        for i in 1:length(lj)
            for j in nblist[i]
                j == 0 && break
                if norm(ps[i] - ps[j]) < 1.7σ
                    for k in nblist[i]
                        k == 0 && break
                        if j < k
                            f = swstrength * swforce(ps[i], ps[j], ps[k], 1.7σ, swangle)
                            fs[i] += f[1]; fs[j] += f[2]; fs[k] += f[3]
                        end
                    end
                end
            end
        end
    end

    for (i, j) in nbs
        if i > 0 && j > 0 && !bonded[i,j]
            v = ps[j] - ps[i]
            nv = norm(v)
            if nv < cutoff
                f = -ljf(σ * (σs[i] + σs[j]) / 2, ε, v)
                fs[i] += f
                fs[j] += -f
                
                f = -1.0 * coulomb * cs[i] * cs[j] * (1 / nv^3 - 1 / cutoff^3) * v
                fs[i] += f
                fs[j] += -f
            end
        end
    end

    for (i, j, k, l) in bonds
        if i <= length(ps) && j <= length(ps)
            v = ps[j] - ps[i]
            f = bondk * k * normalize(v) * (norm(v) - l * bondl)
            fs[i] += f
            fs[j] += -f
        end
    end
    for (i, j, k, kα, α) in angles
        v = ps[j] - ps[i]
        w = ps[k] - ps[i]
        β = angle(v, w)
        f = kα * (β - α)
        vn = normalize(SA[v[2], -v[1]])
        wn = normalize(SA[w[2], -w[1]])
        fj = vn * sign(dot(vn, w)) * f / max(0.01, norm(v)) * pi / 180
        fk = wn * sign(dot(wn, v)) * f / max(0.01, norm(w)) * pi / 180
        fs[i] -= fj + fk
        fs[j] += fj
        fs[k] += fk
    end
    for (i, p) in enumerate(ps)
        fs[i] += -wallk .* sign.(p) .* max.(0, abs.(p) .- wallr)
    end
    fs .+= ms .* Ref(SA[0, g])
    fs .+= coulomb .* cs .* Ref(SA[0, E])
    
    fs
end

function cellList(ps, cutoff)
    xmin, xmax = extrema(p -> floor(Int, p[1] / cutoff), ps)
    ymin, ymax = extrema(p -> floor(Int, p[2] / cutoff), ps)
    
    cells = OffsetArray([Tuple{Int, eltype(ps)}[] for i in xmin-1:xmax+1, j in ymin-1:ymax+1], xmin-1:xmax+1, ymin-1:ymax+1)
    for (i, p) in enumerate(ps)
        push!(cells[floor.(Int, p / cutoff)...], (i, p))
    end
    cells
end

function neighbors!(lj::LJ; cutoff = 0.0)
    if cutoff == 0
        cutoff = 5lj.σ
    end
    cutoff2 = cutoff^2

    lj.nbs .= Ref((0, 0))

    if length(lj) == 0
        return lj.nbs
    end
    
    cells = cellList(lj.ps, cutoff)
    
    dirs = CartesianIndex.(SA[(0, 0), (1, 0), (1, 1), (0, 1), (1, -1)])
    
    n = 0
    @inbounds for I in CartesianIndices(cells)
        isempty(cells[I]) && continue
        for D in dirs
            J = I + D
            for (i, p) in cells[I], (j, q) in cells[J]
                if I != J || i < j
                    v = p - q
                    if dot(v, v) < cutoff2
                        n += 1
                        if n <= length(lj.nbs)
                            lj.nbs[n] = (i, j)
                        else
                            push!(lj.nbs, (i, j))
                        end
                    end
                end
            end
        end
    end
    
    fill!.(lj.nblist, 0)
    if length(lj.nblist) < length(lj)
        append!(lj.nblist, [Int[] for i in length(lj.nblist)+1:length(lj)])
    end
    for (i, j) in lj.nbs
        if i > 0 && j > 0
            n = findfirst(isequal(0), lj.nblist[i])
            if isnothing(n)
                push!(lj.nblist[i], j)
            else
                lj.nblist[i][n] = j
            end
            n = findfirst(isequal(0), lj.nblist[j])
            if isnothing(n)
                push!(lj.nblist[j], i)
            else
                lj.nblist[j][n] = i
            end
        end
    end

    lj.nbs
end

function ensureNeighbors!(lj::LJ; forced = false)
    cutoff = getcutoff(lj)
    if length(lj.lastnbps) != length(lj.ps)
        lj.lastnbps = zero(lj.ps)
    end
    if forced || maximum(norm.(lj.lastnbps .- lj.ps)) > cutoff / 4
        neighbors!(lj, cutoff = cutoff)
        lj.lastnbps .= lj.ps
    end
    lj
end

function step!(lj::LJ; dt = 1e-3)
    lj.ps .*= min.(1, 2lj.wallr ./ norm.(lj.ps))
    lj.vs .*= min.(1, 1000 ./ norm.(lj.vs))

    lj.ps .+= 0.5dt .* lj.vs
    ensureNeighbors!(lj)
    forces!(lj, cutoff = getcutoff(lj) / 2)
    lj.vs .+= dt .* lj.fs ./ lj.ms
    lj.ps .+= 0.5dt .* lj.vs
    lj.nsteps += 1
end

function thermostat!(lj::LJ; dt = 1e-3)
    T = temperature(lj)
    dT = dt * lj.ts * (lj.T - T)
    lj.vs .*= sqrt((T + dT) / T)
    T - temperature(lj)
end

function minimize!(lj::LJ; nsteps = 1000, d = 1e-3)
    for i in 1:nsteps
        lj.ps .*= min.(1, 2lj.wallr ./ norm.(lj.ps))
        i % 10 == 1 && neighbors!(lj)
        fs = forces!(lj)
        fs .*= 1 ./ max.(1, norm.(fs)) .* d
        lj.ps .+= fs
    end
end

function temperature(lj::LJ)
    1/2length(lj) * sum(lj.ms .* norm.(lj.vs).^2, init = 0.0)
end

function potential(lj::LJ; cutoff = Inf)
    (; ps, vs, ms, fs, cs, nbs, bonds, angles, σ, σs, ε, coulomb, wallr, wallk, bondk, bondl, g, E, nblist, swstrength, swangle, bonded) = lj
    
    E = 0.0
    for (i, j) in nbs
        if i > 0 && !bonded[i,j]
            v = ps[j] - ps[i]
            nv = norm(v)
            if nv < cutoff
                E += ljp(σ * (σs[i] + σs[j]) / 2, ε, v) - ljp(σ * (σs[i] + σs[j]) / 2, ε, cutoff)
                
                E += coulomb * cs[i] * cs[j] * (1 / norm(v) - 1 / cutoff - 0.5norm(v)^2)
            end
        end
    end
    if swstrength != 0
        for i in 1:length(lj)
            for j in nblist[i]
                j == 0 && break
                if norm(ps[i] - ps[j]) < 1.7σ
                    for k in nblist[i]
                        k == 0 && break
                        if j < k
                            E += swstrength * swpotential(ps[i], ps[j], ps[k], 1.7σ, swangle)
                        end
                    end
                end
            end
        end
    end
    for (i, j, k, l) in bonds
        if i <= length(ps) && j <= length(ps)
            v = ps[j] - ps[i]
            E += 0.5bondk * k * (norm(v) - l * bondl)^2
        end
    end
    E / length(lj)
end

function potentialPerParticle(lj::LJ)
    (; ps, vs, ms, fs, cs, nbs, bonds, angles, σ, σs, ε, coulomb, wallr, wallk, bondk, bondl, g, E, nblist, swstrength, swangle, bonded) = lj

    cutoff = getcutoff(lj)
    
    Es = zeros(length(lj))
    for (i, j) in nbs
        if i > 0 && j > 0 && !bonded[i,j]
            v = ps[j] - ps[i]
            nv = norm(v)
            if nv < cutoff
                E = ljp(σ * (σs[i] + σs[j]) / 2, ε, v) - ljp(σ * (σs[i] + σs[j]) / 2, ε, cutoff)
                E += coulomb * cs[i] * cs[j] * (1 / norm(v) - 1 / cutoff - 0.5nv^2)
                Es[i] += E; Es[j] += E
            end
        end
    end
    if swstrength != 0
        for i in 1:length(lj)
            for j in nblist[i]
                j == 0 && break
                if norm(ps[i] - ps[j]) < 1.7σ
                    for k in nblist[i]
                        k == 0 && break
                        if j < k
                            E = swstrength * swpotential(ps[i], ps[j], ps[k], 1.7σ, swangle)
                            Es[i] += E; Es[j] += E; Es[k] += E
                        end
                    end
                end
            end
        end
    end
    for (i, j, k, l) in bonds
        if i <= length(ps) && j <= length(ps)
            v = ps[j] - ps[i]
            E = 0.5bondk * k * (norm(v) - l * bondl)^2
            Es[i] += E; Es[j] += E
        end
    end
    Es
end

# bij = zeros(length(lj))
# for (i, j) in nbs
#     if i != 0 && j != 0
#         nv = norm(lj.ps[i] - lj.ps[j])
#         b = smoothstep(nv, σ, 1.5σ, 1, 0)
#         bij[i] += b; bij[j] += b
#     end
# end
# bij .= smoothstep.(bij, 2, 3, 1, 0)