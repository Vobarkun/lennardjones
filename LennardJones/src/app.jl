function remove_interactions!(ax)
    deregister_interaction!(ax, :rectanglezoom)
    deregister_interaction!(ax, :dragpan)
    deregister_interaction!(ax, :scrollzoom)
    deregister_interaction!(ax, :limitreset)
end

function maptounitrange(X)
    min, max = extrema(X)
    (X .- min) ./ (max - min)
end

function main()
    lj = LJ(500; σ = 0.05, ε = 1, T = 1e-1, coulomb = 1.0, ts = 1.0)
    minimize!(lj, nsteps = 100, d = 1e-3)
    lj.vs .= 0.005 .* randn.(SVec2)
    step!(lj, dt = 0)
    lk = ReentrantLock()

    fig = Figure(resolution = (2700,1800), figure_padding = 20)
    main_axis = Axis(fig[1,1], aspect = DataAspect())
    xlims!(main_axis, -1.03, 1.03); ylims!(main_axis, -1.03, 1.03); hidedecorations!(main_axis)
    remove_interactions!(main_axis)

    gl = fig[1,2] = GridLayout()
    gl2 = GridLayout(gl[4,1], alignmode = Outside())
    colsize!(fig.layout, 1, Relative(0.666)); colsize!(fig.layout, 2, Relative(0.333))

    menugrid = gl[2,1] = GridLayout(tellwidth=false, halign = :left)
    startbutton = Button(menugrid[1,1], label = "Start", width = 120)
    freezebutton = Button(menugrid[1,2], label = "Freeze", width = 80)
    on(n -> (@lock lk lj.vs .*= 0), freezebutton.clicks)
    particlemenu = Menu(menugrid[1,3], options = ["default", "heavy", "positive", "negative", "polar pair", "neutral pair"], width = 150)
    colormenu = Menu(menugrid[1,4], options = ["potential", "velocity", "charge", "nothing"], width = 150)


    maxN = 1000
    node = Observable(lj)
    color = Observable(zeros(maxN))
    pnode = Observable(fill(SA[0.0,0.0], maxN))
    markersize = Observable(fill(lj.σ, 1000))
    on(node) do lj
        N = length(lj)
        pnode[][1:N] .= lj.ps
        pnode[][N+1:end] .= Ref(SA[NaN, NaN])
        pnode[] = pnode[]

        if colormenu.selection[] == "charge"
            color[][1:N] .= maptounitrange(lj.cs)
        elseif colormenu.selection[] == "velocity"
            color[][1:N] .= maptounitrange(norm.(lj.vs))
        elseif colormenu.selection[] == "potential"
            color[][1:N] .= sqrt.(maptounitrange(potentialPerParticle(lj) ./ lj.σs))
        elseif colormenu.selection[] == "nothing"
            color[][1:N] .= 0
        end
        color[] = color[]

        if any(markersize[][i] != lj.σ * lj.σs[i] for i in 1:N)
            markersize[][1:N] .= lj.σ .* lj.σs
            markersize[] = markersize[]
        end
    end
    node[] = ensureNeighbors!(deepcopy(lj))

    linesegments!(main_axis, lift(lj -> lj.ps[[b[i] for b in lj.bonds for i in 1:2]], node), strokewidth = lift(lj -> 0.1lj.σ, node), color = :grey)
    main_plot = scatter!(main_axis, pnode, markersize = markersize, strokewidth = 1, color = color, markerspace = :data, colorrange = (0, 1))



    temperatures = Observable(fill(NaN, 10000))
    temperature_axis = Axis(gl2[1,1], ylabel = "temperature")
    lines!(temperature_axis, -9999:1:0, temperatures)
   
    potentials = Observable(fill(NaN, 10000))
    potential_axis = Axis(gl2[2,1], ylabel = "potential energy")
    lines!(potential_axis, -9999:1:0, potentials)

    nneighbors = lift(node) do lj
        nbs = [(i,j) for (i,j) in lj.nbs if i != 0 && norm(lj.ps[i] - lj.ps[j]) < 1.5lj.σ]
        nbs = [first.(nbs); last.(nbs)]
        count.(isequal.(1:length(lj)), Ref(nbs))
    end
    neighbor_axis = Axis(gl2[3,1])
    hist!(neighbor_axis, nneighbors, bins = -0.5:1:11, normalization = :probability)
    hlines!(neighbor_axis, 0.4, color = :transparent)

    distances = lift(lj -> [norm(lj.ps[i] - lj.ps[j]) / lj.σ for (i, j) in lj.nbs if i != 0], node)
    distance_axis = Axis(gl2[4,1])
    hist!(distance_axis, distances, bins = lift(lj -> 0:0.025:5, node), normalization = :probability)
    hlines!(distance_axis, 0.15, color = :transparent)

    on(node) do lj
        temperatures[] = [temperatures[][2:end]; temperature(lj)]
        potentials[] = [potentials[][2:end]; potential(lj, cutoff = 2.5lj.σ)]
        autolimits!(temperature_axis); autolimits!(potential_axis); autolimits!(neighbor_axis); autolimits!(distance_axis)
    end    
    
    remove_interactions!(temperature_axis); remove_interactions!(potential_axis); remove_interactions!(neighbor_axis); remove_interactions!(distance_axis); 
    on(events(potential_axis).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press && (is_mouseinside(potential_axis) || is_mouseinside(temperature_axis))
            potentials[] .= NaN; potentials[] = potentials[]; autolimits!(potential_axis)
            temperatures[] .= NaN; temperatures[] = temperatures[]; autolimits!(temperature_axis)
        end
    end

    sps = Observable(0.0); fps = Observable(0); lastnsteps = Observable(0); nframes = Observable(0); lastframe = Observable(time_ns())
    Label(gl2[5,1], halign = :left, text = lift(x -> "time steps per second: $(round(Int, x))", sps), tellwidth = false)
    Label(gl2[5,1], halign = :right, text = lift(x -> "frames per second: $x", fps), tellwidth = false)
    Label(gl2[5,1], halign = :center, text = lift(x -> "number of particles: $(length(x))", node), tellwidth = false)
    on(node) do lj
        nframes[] += 1
        if time_ns() - lastframe[] > 1e9
            sps[] = (lj.nsteps - lastnsteps[]) / (time_ns() - lastframe[]) * 1e9
            lastnsteps[] = lj.nsteps; fps[] = nframes[]; nframes[] = 0; lastframe[] = time_ns()
        end
    end


    
    presets = gl[3,1] = GridLayout(tellwidth=false, halign = :left)
    begin
        ions = Button(presets[1,1], label = "Ions", width = 80)
        on(ions.clicks) do _
            @lock lk begin
                N = length(lj)
                empty!(lj)
                append!(lj, ps = randb(lj.wallr, N), cs = rand([-1.0, 1.0], N))
                minimize!(lj, nsteps = 1000, d = 1e-3)
                lj.vs .= 0.005 .* randn.(SVec2)
            end
            node[] = ensureNeighbors!(deepcopy(lj))
        end
        pairs = Button(presets[1,2], label = "bonded pairs", width = 140)
        on(pairs.clicks) do _
            @lock lk begin
                N = length(lj) ÷ 2
                empty!(lj)
                for i in 1:N
                    ps = 0.9randb(lj.wallr, 2)
                    ps[2] = ps[1] .+ 0.01randn.()
                    push!(lj, ps = ps, bonds = [(1, 2, 10000.0, 0.05)])
                end
                minimize!(lj, nsteps = 1000, d = 1e-3)
                lj.vs .= 0.005 .* randn.(SVec2)
            end
            node[] = ensureNeighbors!(deepcopy(lj))
        end
        pairs = Button(presets[1,3], label = "polar pairs", width = 130)
        on(pairs.clicks) do _
            @lock lk begin
                N = length(lj) ÷ 2
                empty!(lj)
                for i in 1:N
                    ps = 0.9randb(lj.wallr, 2)
                    ps[2] = ps[1] .+ 0.01randn.()
                    push!(lj, ps = ps, cs = [-1.0, 1.0], bonds = [(1, 2, 10000.0, 0.05)])
                end
                minimize!(lj, nsteps = 1000, d = 1e-3)
                lj.vs .= 0.005 .* randn.(SVec2)
            end
            node[] = ensureNeighbors!(deepcopy(lj))
        end
        chains = Button(presets[1,4], label = "chains", width = 120)
        on(chains.clicks) do _
            @lock lk begin
                N = length(lj) ÷ 5
                empty!(lj)
                for i in 1:N
                    ps = 0.9randb(lj.wallr, 5)
                    for i in 2:5
                        ps[i] = ps[1] .+ 0.01randn.()
                    end
                    append!(lj, ps = ps, cs = zeros(5), 
                        bonds = [(1, 2, 10000.0, 0.05), (2, 3, 10000.0, 0.05), (3, 4, 10000.0, 0.05), (4, 5, 10000.0, 0.05)], 
                        angles = [(2, 1, 3, 10, 180), (3, 2, 4, 10, 180), (4, 3, 5, 10, 180)]
                    )
                end
                minimize!(lj, nsteps = 1000, d = 1e-3)
                lj.vs .= 0.005 .* randn.(SVec2)
            end
            node[] = ensureNeighbors!(deepcopy(lj))
        end
        rings = Button(presets[1,5], label = "rings", width = 120)
        on(rings.clicks) do _
            @lock lk begin
                N = length(lj) ÷ 6
                empty!(lj)
                for i in 1:N
                    ps = 0.9randb(lj.wallr, 6)
                    for i in 2:5
                        ps[i] = ps[1] .+ 0.01randn.()
                    end
                    append!(lj, ps = ps, cs = zeros(6), 
                        bonds = [(1, 2, 10000.0, 0.05), (2, 3, 10000.0, 0.05), (3, 4, 10000.0, 0.05), (4, 5, 10000.0, 0.05), (5, 6, 10000.0, 0.05), (6, 1, 10000.0, 0.05)],
                        angles = [(2, 1, 3, 10, 180), (3, 2, 4, 10, 180), (4, 3, 5, 10, 180), (5, 4, 6, 10, 180), (6, 5, 1, 10, 180), (1, 6, 2, 10, 180)]
                    )
                end
                minimize!(lj, nsteps = 1000, d = 1e-3)
                lj.vs .= 0.005 .* randn.(SVec2)
            end
            node[] = ensureNeighbors!(deepcopy(lj))
        end
    end


    function setnumber(val)
        @lock lk begin
            while length(lj) > val
                deleteat!(lj, length(lj))
            end
            for i in 1:10000
                p = randb(lj.wallr)
                if all(norm(p - q) > lj.σ for q in lj.ps)
                    push!(lj, ps = [p])
                end
                length(lj) >= val && break
            end
        end
        node[] = ensureNeighbors!(deepcopy(lj))
    end

    dt = Observable(1e-4)
    mousestrength = Observable(500)

    settings = GridLayout(gl[1,1])
    Label(settings[1, 1], halign = :left, fontsize = 20, text = "Particle settings", tellwidth = false)
    sliderconf = (
        (label = "number", range = 1:1:1000, startvalue = length(lj)) => setnumber,
        (label = "size", range = 0.01:0.001:0.5, startvalue = lj.σ) => (val -> (lj.σ = val)),
    )
    on.(last.(sliderconf), getfield.(SliderGrid(settings[2,1], first.(sliderconf)...).sliders, :value))

    Label(settings[3, 1], halign = :left, fontsize = 20, text = "Thermostat settings", tellwidth = false)
    sliderconf = (
        (label = "temperature", range = 0:0.001:1, startvalue = lj.T) => val -> (lj.T = val),
        (label = "coupling", range = (0:0.01:10).^3, startvalue = lj.ts) => (val -> (lj.ts = val)),
    )
    on.(last.(sliderconf), getfield.(SliderGrid(settings[4,1], first.(sliderconf)...).sliders, :value))

    Label(settings[5, 1], halign = :left, fontsize = 20, text = "Internal forces", tellwidth = false)
    sliderconf = (
        (label = "Coulomb", range = (0:0.001:10).^2, startvalue = lj.coulomb) => (val -> (lj.coulomb = val)),
        (label = "Lennard-Jones", range = 0:0.01:10, startvalue = lj.ε) => (val -> (lj.ε = val)),
        (label = "Stillinger-Weber", range = 0:0.01:10, startvalue = lj.swstrength) => (val -> (lj.swstrength = val)),
        (label = "covalent angle", range = 0:180, startvalue = lj.swangle) => (val -> (lj.swangle = val)),
        (label = "bond strength", range = (0:0.01:10), startvalue = lj.bondk) => (val -> (lj.bondk = val)),
        (label = "bond length", range = (0:0.01:2), startvalue = lj.bondl) => (val -> (lj.bondl = val)),
    )
    on.(last.(sliderconf), getfield.(SliderGrid(settings[6,1], first.(sliderconf)...).sliders, :value))

    Label(settings[7, 1], halign = :left, fontsize = 20, text = "External forces", tellwidth = false)
    sliderconf = (
        (label = "gravity", range = 0:-0.01:-10, startvalue = lj.g) => (val -> (lj.g = val)),
        (label = "voltage", range = (0:0.1:100).^3, startvalue = lj.g) => (val -> (lj.E = val)),
        (label = "interactive", range = 1:1:1000, startvalue = 500) => (val -> (mousestrength[] = val)),
    )
    on.(last.(sliderconf), getfield.(SliderGrid(settings[8,1], first.(sliderconf)...).sliders, :value))

    Label(settings[9, 1], halign = :left, fontsize = 20, text = "Simulation Controls", tellwidth = false)
    sliderconf = (
        (label = "time step", range = 0:1e-5:1e-3, startvalue = 1e-4) => (val -> (dt[] = val)),
    )
    on.(last.(sliderconf), getfield.(SliderGrid(settings[10,1], first.(sliderconf)...).sliders, :value))



    interactionIndex = Observable(1)
    on(events(main_axis).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            interactionIndex[] = argmin(norm.(lj.ps .- Ref(mouseposition(main_axis))))
        end
    end
    scatter!(main_axis, lift(lj -> lj.ps[min(length(lj.ps), interactionIndex[]):min(length(lj.ps), interactionIndex[])], node), markersize = lift(lj -> 1.05lj.σ, node), 
        markerspace = :data, color = lift(lj -> ifelse(ispressed(main_axis, Mouse.left) && is_mouseinside(main_axis), :red, :transparent), node))

    
    running = Ref(false)
    interactions::Dict{Symbol, Bool} = Dict(:pullsingle => false, :spawn => false, :pullall => false, :delete => false, :slow => false)
    mousepos::SVec2 = SA[0.0,0.0]
    
    function runfunc()
        lastaction = time_ns()
        while true
            try
                t = time_ns()

                @lock lk begin 
                    step!(lj, dt = dt[])
                    thermostat!(lj, dt = dt[])

                    if interactions[:pullsingle] && interactionIndex[] <= length(lj)
                        lj.vs[interactionIndex[]] += mousestrength[] * dt[] * (mousepos - lj.ps[interactionIndex[]]) / lj.ms[interactionIndex[]]
                    end
                    if interactions[:pullall]
                        lj.vs .+= mousestrength[] / 10 * dt[] .* (Ref(mousepos) .- lj.ps) ./ (0.05 .+ norm.(Ref(mousepos) .- lj.ps)).^2 ./ lj.ms
                    end
                    if interactions[:slow]
                        lj.vs .*= ifelse.(norm.(lj.ps .- Ref(mousepos)) .< 0.05, 0.99, 1.0)
                    end

                    if time_ns() - lastaction > 1e8
                        if interactions[:spawn]
                            p = mousepos
                            if length(lj) < maxN && minimum(norm.(Ref(p) .- lj.ps)) > lj.σ
                                if particlemenu.selection[] == "default"
                                    push!(lj, ps = [p])
                                elseif particlemenu.selection[] == "positive"
                                    push!(lj, ps = [p], cs = [1], σs = [2])
                                elseif particlemenu.selection[] == "negative"
                                    push!(lj, ps = [p], cs = [-1], σs = [2])
                                elseif particlemenu.selection[] == "heavy" && minimum(norm.(Ref(p) .- lj.ps)) > sqrt(10) * lj.σ
                                    push!(lj, ps = [p], ms = [100], σs = [sqrt(10)])
                                elseif particlemenu.selection[] == "polar pair" && length(lj) < maxN - 1
                                    d = normalize(randn(SVec2)) * lj.σ / 2
                                    push!(lj, ps = [p + d, p - d], cs = [-0.3, 0.3], ms = [2, 2], bonds = [(1, 2, 10000.0, 0.05)])
                                elseif particlemenu.selection[] == "neutral pair" && length(lj) < maxN - 1
                                    d = normalize(randn(SVec2)) * lj.σ / 2
                                    push!(lj, ps = [p + d, p - d], cs = [0, 0], ms = [2, 2], bonds = [(1, 2, 10000.0, 0.05)])
                                end
                            end
                        end
                        if interactions[:delete] && length(lj) > 1
                            deleteat!(lj, argmin(norm.(lj.ps .- Ref(mousepos))))
                            interactions[:delete] = false
                        end
                        lastaction = time_ns()
                    end
                end
                
                u = time_ns()
                while time_ns() - u < 1e4 end
                while time_ns() - t < 1e5 end

                if !running[]
                    return
                end
            catch e
                stacktrace(catch_backtrace())
                running[] = false
                return
            end
        end
    end
        
    function renderfunc()
        lastplaced = now()
        for i in 1:100000000
            try
                node[] = ensureNeighbors!(deepcopy(lj), forced = true)

                interactions[:pullsingle] = is_mouseinside(main_axis) && ispressed(main_axis, Mouse.left)
                interactions[:spawn] = is_mouseinside(main_axis) && ispressed(main_axis, Mouse.right)
                interactions[:pullall] = is_mouseinside(main_axis) && ispressed(main_axis, Mouse.middle)
                interactions[:delete] = is_mouseinside(main_axis) && ispressed(main_axis, Keyboard.x)
                interactions[:slow] = is_mouseinside(main_axis) && ispressed(main_axis, Keyboard.s)
                mousepos = SVec2(mouseposition(main_axis))

                if !running[] || !events(fig).window_open[]
                    running[] = false
                    startbutton.label = "Start"
                    return
                end
                sleep(0.001)
                yield()
            catch e
                stacktrace(catch_backtrace())
                running[] = false
                return
            end
        end
    end

    on(startbutton.clicks) do n
        running[] = !running[]
        startbutton.label = ifelse(running[], "Stop", "Start")
        if running[]
            Threads.@spawn runfunc()
            schedule(Task(renderfunc))
        end
    end

    display(fig)
    fig, main_plot
end