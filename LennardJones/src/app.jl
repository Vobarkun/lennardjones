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

function updateevery(innode, dt)
    outnode = Observable(innode[])
    t = Observable(time())
    on(innode) do val
        if time() - t[] > dt
            t[] = time()
            outnode[] = val
        end
    end
    outnode
end

function main()
    lj = LJ(500)
    step!(lj, dt = 0.0)
    lk = ReentrantLock()

    fig = Figure(resolution = (3200,1800), figure_padding = 20)
    main_axis::Axis = Axis(fig[1,1], aspect = DataAspect())
    xlims!(main_axis, -1.03, 1.03); ylims!(main_axis, -1.03, 1.03); hidedecorations!(main_axis)
    remove_interactions!(main_axis)

    rightarea = fig[1,2] = GridLayout()
    plotgrid = GridLayout(rightarea[3,1], alignmode = Outside())
    colsize!(fig.layout, 1, Aspect(1, 1.0))#; colsize!(fig.layout, 2, Relative(7/16))

    menugrid = rightarea[2,1] = GridLayout(tellwidth=false, halign = :left)
    
    startbutton::Button = Button(menugrid[1,1], label = "Start", width = 120)
    freezebutton::Button = Button(menugrid[1,2], label = "Freeze", width = 80)
    on(n -> (@lock lk lj.vs .*= 0), freezebutton.clicks)
    
    Label(menugrid[1,3], " Right click:", fontsize = 20, justification = :right)
    particlemenu::Menu = Menu(menugrid[1,4], options = ["default", "heavy", "very heavy", "positive", "negative", "polar pair", "neutral pair"], width = 120)
    colgap!(menugrid, 3, 10)
    
    Label(menugrid[1,5], " Color:", fontsize = 20, justification = :right)
    colormenu::Menu = Menu(menugrid[1,6], options = ["potential", "velocity", "charge", "nothing"], width = 120)
    colgap!(menugrid, 5, 10)


    maxN = 1000
    node = Observable(lj)
    color = Observable(zeros(maxN))
    colorrange = Observable((0.0, 1.0))
    pnode = Observable(fill(SA[0.0,0.0], maxN))
    markersize = Observable(fill(lj.σ, 1000))
    on(node) do lj
        N = length(lj)
        pnode[][1:N] .= lj.ps
        pnode[][N+1:end] .= Ref(SA[NaN, NaN])
        pnode[] = pnode[]

        if colormenu.selection[] == "charge"
            color[][1:N] .= lj.cs
        elseif colormenu.selection[] == "velocity"
            color[][1:N] .= 0.5 .* color[][1:N] .+ 0.5 .* norm.(lj.vs)
        elseif colormenu.selection[] == "potential"
            cs = potentialPerParticle(lj) ./ lj.σs
            color[][1:N] .= sqrt.(cs .- minimum(cs))
        elseif colormenu.selection[] == "nothing"
            color[][1:N] .= 0
        end
        color[] = color[]
        colorrange[] = (0.9colorrange[][1] + 0.1minimum(color[]), 0.9colorrange[][2] + 0.1maximum(color[]))

        if any(markersize[][i] != lj.σ * lj.σs[i] for i in 1:N)
            markersize[][1:N] .= lj.σ .* lj.σs
            markersize[] = markersize[]
        end
    end
    node[] = ensureNeighbors!(deepcopy(lj))
    colorrange[] = extrema(color[])

    linesegments!(main_axis, lift(lj -> lj.ps[[b[i] for b in lj.bonds for i in 1:2]], node), strokewidth = lift(lj -> 0.1lj.σ, node), color = :grey)
    main_plot::Scatter{Tuple{Vector{Point{2, Float32}}}} = scatter!(main_axis, pnode, markersize = markersize, strokewidth = 1, color = color, markerspace = :data, colormap = :isoluminant_cm_70_c39_n256, colorrange = colorrange)


    # settings
    begin 
        function setnumber(val)
            @lock lk begin
                while length(lj) > val
                    deleteat!(lj, length(lj))
                end
                while length(lj) < val
                    push!(lj, randfree(lj))
                end
                node[] = ensureNeighbors!(deepcopy(lj))
            end
        end

        dt = Observable(1e-4)
        mousestrength = Observable(500)
        numberofatoms = Observable(length(lj))

        settings = GridLayout(rightarea[1,1])
        Label(settings[1, 1], halign = :left, fontsize = 20, text = "Particle settings", tellwidth = false)
        sliderconf = (
            (label = "number", range = 1:1:1000, startvalue = length(lj)) => (val -> (numberofatoms[] = val; setnumber(val))),
            (label = "size", range = 0.01:0.001:0.5, startvalue = lj.σ) => (val -> (lj.σ = val)),
        )
        on.(last.(sliderconf), getfield.(SliderGrid(settings[2,1], first.(sliderconf)...).sliders, :value))

        Label(settings[3, 1], halign = :left, fontsize = 20, text = "Thermostat settings", tellwidth = false)
        sliderconf = (
            (label = "temperature", range = 0:0.001:1, startvalue = lj.T) => val -> (lj.T = val),
            (label = "strength", range = (0:0.01:10).^3, startvalue = lj.ts) => (val -> (lj.ts = val)),
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
        mousestrengthslider::Slider = contents(settings[8,1])[1].sliders[3]

        Label(settings[9, 1], halign = :left, fontsize = 20, text = "Simulation Controls", tellwidth = false)
        sliderconf = (
            (label = "time step", range = range(0, sqrt(1e-3), length = 1000).^2, startvalue = 1e-4) => (val -> (dt[] = val)),
        )
        on.(last.(sliderconf), getfield.(SliderGrid(settings[10,1], first.(sliderconf)...).sliders, :value))
    end


    #presets
    begin
        options = [
            ("Default", N -> pushfree!(lj, N)),
            ("Ions", N -> for i in 1:N pushfree!(lj, cs = [-sign(sum(lj.cs) - 1e-10)]) end),
            ("diatomic", N -> pushfree!(lj, N ÷ 2, 
                ps = [SA[lj.bondl * 0.05, 0.0], SA[0.0, 0.0]], 
                bonds = [(1, 2, 10000.0, 0.05)]
            )),
            ("polar", N -> pushfree!(lj, N ÷ 2, 
                ps = [SA[lj.bondl * 0.05, 0.0], SA[0.0, 0.0]], 
                cs = [1.0, -1.0], bonds = [(1, 2, 10000.0, 0.05)]
            )),
            ("chains", N -> pushfree!(lj, N ÷ 5, 
                ps = [SA[i * lj.bondl * 0.05, 0.0] for i in 1:5], 
                bonds = [(1, 2, 10000.0, 0.05), (2, 3, 10000.0, 0.05), (3, 4, 10000.0, 0.05), (4, 5, 10000.0, 0.05)], 
                angles = [(2, 1, 3, 10, 180), (3, 2, 4, 10, 180), (4, 3, 5, 10, 180)]
            )),
            ("rings", N -> pushfree!(lj, N ÷ 6, 
                ps = [lj.bondl * 0.05 * SA[cos(ϕ), sin(ϕ)] for ϕ in 0:2pi/6:6], 
                bonds = [(1, 2, 10000.0, 0.05), (2, 3, 10000.0, 0.05), (3, 4, 10000.0, 0.05), 
                         (4, 5, 10000.0, 0.05), (5, 6, 10000.0, 0.05), (6, 1, 10000.0, 0.05)],
                angles = [(2, 1, 3, 10, 180), (3, 2, 4, 10, 180), (4, 3, 5, 10, 180), 
                          (5, 4, 6, 10, 180), (6, 5, 1, 10, 180), (1, 6, 2, 10, 180)]
            )),
        ]
        Label(menugrid[1, 7], "  Presets:", fontsize = 20, justification = :right)
        on(Menu(menugrid[1, 8], options = options, width = 120).selection) do func
            lock(lk) do 
                N = numberofatoms[]; empty!(lj)
                func(N)
                node[] = ensureNeighbors!(deepcopy(lj))
            end
        end
        colgap!(menugrid, 7, 10)

        sps = Observable(0.0); fps = Observable(0); lastnsteps = Observable(0); nframes = Observable(0); lastframe = Observable(time_ns())
        Label(menugrid[1,9], halign = :right, text = lift((s, f, lj) -> "tps: $(round(Int, s))  fps: $f  N: $(length(lj))", sps, fps, node), tellwidth = false)
        on(node) do lj
            nframes[] += 1
            if time_ns() - lastframe[] > 1e9
                sps[] = (lj.nsteps - lastnsteps[]) / (time_ns() - lastframe[]) * 1e9
                lastnsteps[] = lj.nsteps; fps[] = nframes[]; nframes[] = 0; lastframe[] = time_ns()
            end
        end
    end


    # plots
    begin
        histsteps = 1000
        Label(plotgrid[1, 1], "Temperature", rotation = pi/2, tellheight = false, fontsize = 20)
        temperatures = Observable(fill(NaN, histsteps))
        temperature_axis::Axis = Axis(plotgrid[1,2], ylabelsize = 20)
        hlines!(temperature_axis, [0, 0.1], color = :transparent)
        lines!(temperature_axis, -histsteps+1:1:0, temperatures)
    
        Label(plotgrid[2, 1], "Potential energy", rotation = pi/2, tellheight = false, fontsize = 20)
        potentials = Observable(fill(NaN, histsteps))
        potential_axis::Axis = Axis(plotgrid[2,2], ylabelsize = 20)
        hlines!(potential_axis, 0, color = :transparent)
        lines!(potential_axis, -histsteps+1:1:0, potentials)

        Label(plotgrid[3, 1], "Number of neighbors", rotation = pi/2, tellheight = false, fontsize = 20)
        nneighbors = lift(node) do lj
            nbs = [(i,j) for (i,j) in lj.nbs if i != 0 && norm(lj.ps[i] - lj.ps[j]) < 1.5lj.σ]
            nbs = [first.(nbs); last.(nbs)]
            count.(isequal.(1:length(lj)), Ref(nbs))
        end
        neighbor_axis::Axis = Axis(plotgrid[3,2], xticks = 0:10, yticks = 0:0.2:1, ylabelsize = 20)
        hist!(neighbor_axis, nneighbors, bins = -0.5:1:11, normalization = :probability)
        hlines!(neighbor_axis, 0.4, color = :transparent)

        Label(plotgrid[4, 1], "Distance distribution", rotation = pi/2, tellheight = false, fontsize = 20)
        distances = lift(lj -> [norm(lj.ps[i] - lj.ps[j]) / lj.σ for (i, j) in lj.nbs if i != 0], node)
        distance_axis::Axis = Axis(plotgrid[4,2], xticks = MultiplesTicks(6, 1, "σ"), ylabelsize = 20)
        hist!(distance_axis, distances, bins = lift(lj -> 0:0.05:5, node), normalization = :probability)
        hlines!(distance_axis, 0.15, color = :transparent)

        on(node) do lj
            T = temperature(lj); T != temperatures[][end] && (temperatures[] = [temperatures[][2:end]; temperature(lj)])
            V = potential(lj); V != potentials[][end] && (potentials[] = [potentials[][2:end]; V])
            autolimits!(temperature_axis); autolimits!(potential_axis); autolimits!(neighbor_axis); autolimits!(distance_axis)
        end    
        
        remove_interactions!(temperature_axis); remove_interactions!(potential_axis); 
        remove_interactions!(neighbor_axis); remove_interactions!(distance_axis); 
        on(events(potential_axis).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press && (is_mouseinside(potential_axis) || is_mouseinside(temperature_axis))
                potentials[] .= NaN; potentials[] = potentials[]; autolimits!(potential_axis)
                temperatures[] .= NaN; temperatures[] = temperatures[]; autolimits!(temperature_axis)
            end
        end
    end



    
    running = Observable(false)
    keymap = (
        :pullsingle => Mouse.left, :spawn => Mouse.right, :pullall => Mouse.middle, :pushall => Keyboard.c,
        :delete => Keyboard.x, :cool => Keyboard.s, :heat => Keyboard.d, :stir => Keyboard.a,
    )
    interactions = Dict(first.(keymap) .=> false)
    mousepos = Observable(SA[0.0,0.0])

    interactionIndex = Observable(1)
    on(events(main_axis).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            interactionIndex[] = argmin(norm.(lj.ps .- Ref(mouseposition(main_axis))))
        end
    end
    scatter!(main_axis, 
        lift(lj -> (i = clamp(interactionIndex[], 1, length(lj)); lj.ps[i:i]), node), 
        markersize = lift(lj -> 1.05lj.σ, node), markerspace = :data, 
        color = lift(lj -> ifelse(interactions[:pullsingle], :red, :transparent), node)
    )

    scrollnode = Observable(0.0)
    on(updateevery(scrollnode, 0.05)) do val
        set_close_to!(mousestrengthslider, round(mousestrengthslider.value[] + 10val, digits = -1))
        scrollnode.val = 0.0
    end 
    register_interaction!(main_axis, :scroll) do event::ScrollEvent, axis
        scrollnode[] += event.y
    end
    
    function runfunc()
        while true
            try
                t = time_ns()

                @lock lk begin 
                    step!(lj, dt = dt[])
                    thermostat!(lj, dt = dt[])

                    if interactions[:pullsingle] && interactionIndex[] <= length(lj)
                        lj.vs[interactionIndex[]] += mousestrength[] * dt[] * (mousepos[] - lj.ps[interactionIndex[]]) / lj.ms[interactionIndex[]]
                    end
                    if interactions[:pullall]
                        lj.vs .+= mousestrength[] / 10 * dt[] .* (Ref(mousepos[]) .- lj.ps) ./ (0.05 .+ norm.(Ref(mousepos[]) .- lj.ps)).^2 ./ lj.ms
                    end
                    if interactions[:pushall]
                        lj.vs .+= mousestrength[] * dt[] .* exp.(-20 .* norm.(lj.ps .- Ref(mousepos[])).^2) .* (lj.ps .- Ref(mousepos[]))
                    end
                    if interactions[:cool]
                        lj.vs .*= 1 .- (0.2mousestrength[] * dt[]) .* exp.(-10 .* norm.(lj.ps .- Ref(mousepos[])).^2)
                    end
                    if interactions[:heat]
                        lj.vs .+= 2mousestrength[] * dt[] .* exp.(-100 .* norm.(lj.ps .- Ref(mousepos[])).^2) .* randn.(SVec2)
                    end
                    if interactions[:stir]
                        lj.vs .+= mousestrength[] / 5 * dt[] .* exp.(-20 .* norm.(lj.ps .- Ref(mousepos[])).^2) .* Ref(SA[0.7 1; -1 0.7]) .* (Ref(mousepos[]) .- lj.ps)
                    end
                end
                
                u = time_ns()
                while time_ns() - u < 1e4 end
                while time_ns() - t < 1e5 end

                if !running[]
                    return
                end
            catch e
                running[] = false
                rethrow(e)
            end
        end
    end
        
    lastaction = Observable(time())

    screen = GLMakie.Screen(renderloop = GLMakie.renderloop, framerate = 120)
    on(screen.render_tick) do x
        node[] = ensureNeighbors!(deepcopy(lj))

        for (name, key) in keymap
            if is_mouseinside(main_axis) || interactions[name]
                interactions[name] = ispressed(main_axis, key)
            end
        end
        mousepos[] = SVec2(mouseposition(main_axis))
        
        if time() - lastaction[] > 0.05
            if interactions[:spawn]
                lock(lk) do 
                    p = mousepos[]
                    if length(lj) < maxN && minimum(norm.(Ref(p) .- lj.ps)) > lj.σ
                        if particlemenu.selection[] == "default"
                            push!(lj, ps = [p])
                        elseif particlemenu.selection[] == "positive"
                            push!(lj, ps = [p], cs = [1], σs = [2])
                        elseif particlemenu.selection[] == "negative"
                            push!(lj, ps = [p], cs = [-1], σs = [2])
                        elseif particlemenu.selection[] == "heavy" && minimum(norm.(Ref(p) .- lj.ps)) > sqrt(10) * lj.σ
                            push!(lj, ps = [p], ms = [100], σs = [sqrt(10)])
                        elseif particlemenu.selection[] == "very heavy" && minimum(norm.(Ref(p) .- lj.ps)) > sqrt(10) * lj.σ
                            push!(lj, ps = [p], ms = [10000], σs = [sqrt(10)])
                        elseif particlemenu.selection[] == "polar pair" && length(lj) < maxN - 1
                            d = normalize(randn(SVec2)) * lj.σ / 2
                            push!(lj, ps = [p + d, p - d], cs = [-0.3, 0.3], ms = [2, 2], bonds = [(1, 2, 10000.0, 0.05)])
                        elseif particlemenu.selection[] == "neutral pair" && length(lj) < maxN - 1
                            d = normalize(randn(SVec2)) * lj.σ / 2
                            push!(lj, ps = [p + d, p - d], cs = [0, 0], ms = [2, 2], bonds = [(1, 2, 10000.0, 0.05)])
                        end
                    end
                end
                lastaction[] = time()
            end
            if interactions[:delete] && length(lj) > 1
                lock(lk) do 
                    deleteat!(lj, argmin(norm.(lj.ps .- Ref(mousepos[]))))
                end
                lastaction[] = time()
            end
        end

        if !running[] || !events(fig).window_open[]
            running[] = false
            startbutton.label = "Start"
            return
        end
    end
    on(screen.window_open) do val
        if !val
            running[] = false
        end
    end

    on(startbutton.clicks) do n
        running[] = !running[]
        startbutton.label = ifelse(running[], "Stop", "Start")
        if running[]
            Threads.@spawn runfunc()
        end
    end

    display(screen, fig)

    fig, main_plot, node, screen
end