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

function logrange(a, b; length)
    exp.(range(log(a), log(b), length = length))
end

function paddedextrema(X; rpad = 0.1, apad = 1e-3)
    a = minimum(x -> ifelse(isnan(x), Inf, x), X)
    b = maximum(x -> ifelse(isnan(x), -Inf, x), X)
    if !isfinite(a) || !isfinite(b)
        return (zero(a), one(b))
    end
    a - rpad * (b - a) - apad, b + rpad * (b - a) + apad
end

function darktheme!()
    Makie.COLOR_ACCENT[] = RGBf(((79, 122, 214) ./ 255 .* 1.5)...)
    Makie.COLOR_ACCENT_DIMMED[] = RGBf(((174, 192, 230) ./ 255 ./ 2)...);
    theme = Theme(
        backgroundcolor = :gray10,
        textcolor = :white,
        palette = (
            color = [Makie.COLOR_ACCENT[]],
        ),
        Hist = (cycle = :color, ),
        Lines = (cycle = :color, ),
        Axis = (
            backgroundcolor = :grey15,
            xgridcolor = (:white, 0.09),
            ygridcolor = (:white, 0.09),
            leftspinevisible = false,
            rightspinevisible = false,
            bottomspinevisible = false,
            topspinevisible = false,
            xminorticksvisible = false,
            yminorticksvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            xlabelpadding = 3,
            ylabelpadding = 3
        ),
        Slider = (
            color_inactive = :grey20,
            text_color = :grey80,
        ),
        Button = (
            buttoncolor = :grey20,
        ),
        Menu = (
            cell_color_inactive_even = :grey20,
            cell_color_inactive_odd = :grey20,
            cell_color_active = RGBf(0.4, 0.4, 0.4),
            color = :grey20,
            backgroundcolor = :grey20,
            textcolor = :white,
            selection_cell_color_inactive = :grey20
        )
    )
    set_theme!(theme);
end

function lighttheme!()
    Makie.COLOR_ACCENT[] = RGBf(((79, 122, 214) ./ 255)...)
    Makie.COLOR_ACCENT_DIMMED[] = RGBf(((174, 192, 230) ./ 255)...);
    set_theme!()
end

# pentagonal crystal: lj 1.8, sw 0.7
# nitinol: coulomb 0.4, lj 1.3, σ 0.35

function main(; decorated = true)
    lj = LJ(500)
    step!(lj, dt = 0.0)
    node = Observable(lj)

    lk = ReentrantLock()

    fig = Figure(resolution = (3200,1800), figure_padding = 20)

    rightarea = fig[1,2] = GridLayout()
    colsize!(fig.layout, 1, Aspect(1, 1.0))


    # settings
    begin 
        function setnumber(val)
            lock(lk) do 
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
        Label(settings[1, 1], halign = :left, fontsize = 25, text = "Particle settings", tellwidth = false)
        sliderconf1 = (
            (label = "number", range = 1:1:1000, startvalue = length(lj)) => (val -> (numberofatoms[] = val; setnumber(val))),
            (label = "size", range = 0.01:0.0001:0.5, startvalue = lj.σ) => (val -> (lj.σ = val)),
        )
        on.(last.(sliderconf1), getfield.(SliderGrid(settings[2,1], first.(sliderconf1)...).sliders, :value))
        sizeslider::Slider = contents(settings[2,1])[1].sliders[2]

        Label(settings[3, 1], halign = :left, fontsize = 25, text = "Thermostat settings", tellwidth = false)
        sliderconf2 = (
            (label = "temperature", range = (0:0.001:2).^log2(10), startvalue = lj.T) => val -> (lj.T = val),
            (label = "strength", range = (0:0.01:10).^3, startvalue = lj.ts) => (val -> (lj.ts = val)),
        )
        on.(last.(sliderconf2), getfield.(SliderGrid(settings[4,1], first.(sliderconf2)...).sliders, :value))

        Label(settings[5, 1], halign = :left, fontsize = 25, text = "Internal forces", tellwidth = false)
        sliderconf3 = (
            (label = "Coulomb", range = (0:0.001:10).^2, startvalue = lj.coulomb) => (val -> (lj.coulomb = val)),
            (label = "Lennard-Jones", range = 0:0.01:10, startvalue = lj.ε) => (val -> (lj.ε = val)),
            (label = "Stillinger-Weber", range = 0:0.01:10, startvalue = lj.swstrength) => (val -> (lj.swstrength = val)),
            (label = "covalent angle", range = 0:180, startvalue = lj.swangle) => (val -> (lj.swangle = val)),
            (label = "bond strength", range = (0:0.01:10), startvalue = lj.bondk) => (val -> (lj.bondk = val)),
            (label = "bond length", range = (0:0.01:2), startvalue = lj.bondl) => (val -> (lj.bondl = val)),
        )
        on.(last.(sliderconf3), getfield.(SliderGrid(settings[6,1], first.(sliderconf3)...).sliders, :value))

        Label(settings[7, 1], halign = :left, fontsize = 25, text = "External forces", tellwidth = false)
        sliderconf4 = (
            (label = "gravity", range = 0:-0.01:-10, startvalue = lj.g) => (val -> (lj.g = val)),
            (label = "voltage", range = (0:0.1:100).^3, startvalue = lj.g) => (val -> (lj.E = val)),
            (label = "magnetic", range = (0:0.1:100), startvalue = lj.g) => (val -> (lj.B = val)),
            (label = "interactive", range = 1:1:1000, startvalue = 500) => (val -> (mousestrength[] = val)),
        )
        on.(last.(sliderconf4), getfield.(SliderGrid(settings[8,1], first.(sliderconf4)...).sliders, :value))

        Label(settings[9, 1], halign = :left, fontsize = 25, text = "Simulation Controls", tellwidth = false)
        sliderconf5 = (
            (label = "time step", range = logrange(1e-6, 1e-3, length = 1000), startvalue = 1e-4) => (val -> (dt[] = val)),
        )
        on.(last.(sliderconf5), getfield.(SliderGrid(settings[10,1], first.(sliderconf5)...).sliders, :value))
    end


    # buttons
    begin
        menugrid = rightarea[2,1] = GridLayout(tellwidth = false, halign = :left)
    
        startbutton::Button = Button(menugrid[1,1], label = "Start", width = 120, fontsize = 20)
        on(n -> lock(() -> lj.vs .*= 0, lk), Button(menugrid[1,2], label = "Freeze", width = 120, fontsize = 20).clicks)

        Label(menugrid[1,4], "   place:", fontsize = 20, justification = :right)
        particlemenu::Menu = Menu(menugrid[1,5], options = ["default", "heavy", "very heavy", "positive", "negative", "polar pair", "neutral pair"], width = 120, fontsize = 20)
        colgap!(menugrid, 4, 10)
        
        Label(menugrid[1,6], " Color:", fontsize = 20, justification = :right)
        colormenu::Menu = Menu(menugrid[1,7], options = ["potential", "velocity", "virial", "charge", "nothing"], width = 120, fontsize = 20)
        colgap!(menugrid, 6, 10)

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
        Label(menugrid[1, 8], "  Presets:", fontsize = 20, justification = :right)
        presetmenu::Menu = Menu(menugrid[1, 9], options = options, width = 120, fontsize = 20)
        on(presetmenu.selection) do func
            lock(lk) do 
                N = numberofatoms[]; empty!(lj)
                func(N)
                node[] = ensureNeighbors!(deepcopy(lj))
            end
        end
        colgap!(menugrid, 8, 10)

        on(Button(menugrid[1,3], label = "Reset", width = 120, fontsize = 20).clicks) do _
            notify(presetmenu.selection)
        end

        sps = Observable(0.0); fps = Observable(0); lastnsteps = Observable(0); nframes = Observable(0); lastframe = Observable(time_ns())
        Label(menugrid[1,10], halign = :right, text = lift((s, f, lj) -> "tps: $(round(Int, s))  fps: $f  N: $(length(lj))", sps, fps, node), tellwidth = false)
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
        plotgrid = GridLayout(rightarea[3,1], alignmode = Outside())

        mcolor = cgrad(:isoluminant_cm_70_c39_n256)[1.0]

        histsteps = 1001
        Label(plotgrid[1, 1], "Temperature", rotation = pi/2, tellheight = false, fontsize = 20)
        temperatures = Observable(fill(NaN, histsteps))
        temperature_axis::Axis = Axis(plotgrid[1,2], ylabelsize = 20)
        xlims!(temperature_axis, -1.02histsteps, 0.02histsteps)
        lines!(temperature_axis, -histsteps+1:1:0, temperatures)

        Label(plotgrid[2, 1], "Potential", rotation = pi/2, tellheight = false, fontsize = 20)
        potentials = Observable(fill(NaN, histsteps))
        potential_axis::Axis = Axis(plotgrid[2,2], ylabelsize = 20)
        xlims!(potential_axis, -1.02histsteps, 0.02histsteps)
        lines!(potential_axis, -histsteps+1:1:0, potentials)

        Label(plotgrid[3, 1], "Pressure", rotation = pi/2, tellheight = false, fontsize = 20)
        pressures = Observable(fill(NaN, histsteps))
        pressure_axis::Axis = Axis(plotgrid[3,2], ylabelsize = 20)
        xlims!(pressure_axis, -1.02histsteps, 0.02histsteps)
        lines!(pressure_axis, -histsteps+1:1:0, pressures)

        # lines!(temperature_axis, -histsteps+1:1:0, lift(x -> imfilter(x, KernelFactors.gaussian(10)), temperatures), linewidth = 2, color = mcolor)
        # lines!(potential_axis, -histsteps+1:1:0, lift(x -> imfilter(x, KernelFactors.gaussian(10)), potentials), linewidth = 2, color = mcolor)
        # lines!(pressure_axis, -histsteps+1:1:0, lift(x -> imfilter(x, KernelFactors.gaussian(10)), pressures), linewidth = 2, color = mcolor)

        Label(plotgrid[4, 1], "Number of neighbors", rotation = pi/2, tellheight = false, fontsize = 20)
        nneighbors = lift(node) do lj
            nbs = [(i,j) for (i,j) in lj.nbs if i != 0 && j != 0 && norm(lj.ps[i] - lj.ps[j]) < 1.5lj.σ]
            nbs = [first.(nbs); last.(nbs)]
            count.(isequal.(1:length(lj)), Ref(nbs))
        end
        neighbor_axis::Axis = Axis(plotgrid[4,2], xticks = 0:10, yticks = 0:0.2:1, ylabelsize = 20)
        hist!(neighbor_axis, nneighbors, bins = -0.5:1:11, normalization = :probability)
        hlines!(neighbor_axis, 0.4, color = :transparent)

        Label(plotgrid[5, 1], "Distance distribution", rotation = pi/2, tellheight = false, fontsize = 20)
        distances = lift(lj -> [norm(lj.ps[i] - lj.ps[j]) / lj.σ for (i, j) in lj.nbs if i != 0 && j != 0], node)
        distance_axis::Axis = Axis(plotgrid[5,2], xticks = MultiplesTicks(6, 1, "σ"), ylabelsize = 20)
        hist!(distance_axis, distances, bins = 0:0.05:5, normalization = :probability)
        hlines!(distance_axis, 0.15, color = :transparent)

        colgap!(plotgrid, 1, 10)

        function updatePlots!(lj)
            T = temperature(lj)
            if T != temperatures[][end]
                temperatures[][1:end-1] .= temperatures[][2:end]
                temperatures[][end] = T
                notify(temperatures)
            end
            V = potential(lj)
            if V != potentials[][end] 
                potentials[][1:end-1] .= potentials[][2:end]
                potentials[][end] = V
                notify(potentials)
            end
            P = T / 4 + virial(lj) / 8 / length(lj)
            if P != pressures[][end]
                pressures[][1:end-1] .= pressures[][2:end]
                pressures[][end] = P#isnan(pressures[][end]) ? P : 0.99pressures[][end] + 0.01P
                notify(pressures)
            end
            
            ylims!(temperature_axis, paddedextrema(temperatures[], rpad = 0.12, apad = 0.0012))
            ylims!(potential_axis, paddedextrema(potentials[], rpad = 0.12, apad = 0.0012))
            ylims!(pressure_axis, paddedextrema(pressures[], rpad = 0.12, apad = 0.0012))
            autolimits!(neighbor_axis); autolimits!(distance_axis)
        end    
        on(updatePlots!, node)
        
        remove_interactions!(temperature_axis); remove_interactions!(potential_axis); 
        remove_interactions!(neighbor_axis); remove_interactions!(distance_axis); 

        function resethistory(event)
            potentials[] .= NaN; 
            temperatures[] .= NaN; 
            pressures[] .= NaN; 
        end        
        onmouseleftdown(resethistory, temperature_axis.mouseeventhandle)
        onmouseleftdown(resethistory, potential_axis.mouseeventhandle)
        onmouseleftdown(resethistory, pressure_axis.mouseeventhandle)
    end

    
    # main plot
    begin
        main_axis::Axis = Axis(fig[1,1], aspect = DataAspect())
        xlims!(main_axis, -1.03, 1.03); ylims!(main_axis, -1.03, 1.03); hidedecorations!(main_axis)
        remove_interactions!(main_axis)

        maxN = 1000
        color = Observable(zeros(maxN))
        colorrange = Observable((0.0, 1.0))
        pnode = Observable(fill(SA[0.0,0.0], maxN))
        markersize = Observable(fill(lj.σ, 1000))

        function plotfunc(lj)
            N = length(lj)
            pnode[][1:N] .= lj.ps
            pnode[][N+1:end] .= Ref(SA[NaN, NaN])
            notify(pnode)

            if colormenu.selection[] == "charge"
                color[][1:N] .= lj.cs
            elseif colormenu.selection[] == "velocity"
                color[][1:N] .= 0.5 .* color[][1:N] .+ 0.5 .* norm.(lj.vs)
            elseif colormenu.selection[] == "potential"
                cs = potentialPerParticle(lj) ./ lj.σs
                color[][1:N] .= sqrt.(cs .- minimum(cs))
            elseif colormenu.selection[] == "virial"
                color[][1:N] .= virialPerParticle(lj)# ./ max.(1, neighborsPerParticle(lj, 1.2lj.σ))
            elseif colormenu.selection[] == "nothing"
                color[][1:N] .= 0
            end
            color[] = color[]
            colorrange[] = (0.8colorrange[][1] + 0.2minimum(color[][1:N]), 0.8colorrange[][2] + 0.2maximum(color[][1:N]))

            if any(markersize[][i] != lj.σ * lj.σs[i] for i in 1:N)
                markersize[][1:N] .= lj.σ .* lj.σs
                notify(markersize)
            end
        end
        on(plotfunc, node) 
        node[] = ensureNeighbors!(deepcopy(lj))
        colorrange[] = extrema(color[])

        linesegments!(main_axis, lift(lj -> lj.ps[[b[i] for b in lj.bonds for i in 1:2]], node), strokewidth = lift(lj -> 0.1lj.σ, node), color = :grey)
        main_plot::Scatter{Tuple{Vector{Point{2, Float32}}}} = scatter!(main_axis, pnode, markersize = markersize, strokewidth = 1, color = color, markerspace = :data, colormap = :isoluminant_cm_70_c39_n256, colorrange = colorrange)
    end


    # interactions
    begin
        keymap = (
            :pullsingle => Mouse.left, :pullall => Mouse.middle, :pushall => Keyboard.c,
            :spawn => Mouse.right, :delete => Keyboard.x, 
            :stirleft => Keyboard.q, :stirright => Keyboard.e,
            :pushleft => Keyboard.a, :pushright => Keyboard.d,
            :pushup => Keyboard.w, :pushdown => Keyboard.s,
            :cool => Keyboard.g, :coolall => Keyboard.f, 
            :heat => Keyboard.t,  :heatall => Keyboard.r, 
        )
        interactions = Dict(first.(keymap) .=> false)
        function handleInteractions!(lj, interactions, index, mousepos, strength, dt)
            if interactions[:pullsingle] && index <= length(lj)
                lj.vs[index] += index * dt * (mousepos - lj.ps[index]) / lj.ms[index]
            end
            if interactions[:pullall]
                lj.vs .+= 0.1strength * dt .* (Ref(mousepos) .- lj.ps) ./ (0.05 .+ norm.(Ref(mousepos) .- lj.ps)).^2 ./ lj.ms
                lj.vs .*= 1 .- (0.1strength * dt) .* exp.(-10 .* norm.(lj.ps .- Ref(mousepos)).^2)
            end
            if interactions[:pushall]
                lj.vs .+= 30strength * dt .* exp.(-200 .* norm.(lj.ps .- Ref(mousepos)).^2) .* (lj.ps .- Ref(mousepos)) ./ lj.ms
            end
            if interactions[:cool]
                lj.vs .*= max.(0, 1 .- (0.2strength * dt) .* exp.(-100 .* norm.(lj.ps .- Ref(mousepos)).^2))
            end
            if interactions[:coolall]
                lj.vs .*= max.(0, 1 .- (0.05strength * dt))
            end
            if interactions[:heat]
                lj.vs .+= 2strength * dt .* exp.(-100 .* norm.(lj.ps .- Ref(mousepos)).^2) .* randn.(SVec2) ./ lj.ms
            end
            if interactions[:heatall]
                lj.vs .+= 1strength * dt .* randn.(SVec2) ./ lj.ms
            end
            if interactions[:stirleft]
                lj.vs .+= strength / 5 * dt .* exp.(-20 .* norm.(lj.ps .- Ref(mousepos)).^2) .* Ref(SA[0.7 1; -1 0.7]) .* (Ref(mousepos) .- lj.ps) ./ lj.ms
            end
            if interactions[:stirright]
                lj.vs .+= strength / 5 * dt .* exp.(-20 .* norm.(lj.ps .- Ref(mousepos)).^2) .* Ref(SA[0.7 -1; 1 0.7]) .* (Ref(mousepos) .- lj.ps) ./ lj.ms
            end
            for (i, dir) in zip((:pushleft, :pushright, :pushup, :pushdown), (SA[-1, 0], SA[1, 0], SA[0, 1], SA[0, -1]))
                if interactions[i]
                    lj.vs .+= 0.05strength .* dt .* Ref(dir) ./ lj.ms
                end
            end
        end
        function updateInteractions!(interactions, keymap, axis)
            for (name, key) in keymap
                if key isa Keyboard.Button || is_mouseinside(axis) || interactions[name]
                    interactions[name] = ispressed(axis, key)
                end
            end
        end
        
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
            set_close_to!(sizeslider, sizeslider.value[] + val)
            scrollnode.val = 0.0
        end 
        register_interaction!(main_axis, :scroll) do event::ScrollEvent, axis
            lj.σ += 0.0001event.y
            scrollnode[] += 0.0001event.y
        end
    end


    # threads
    running = Observable(false)
    function runfunc()
        while true
            try
                t = time_ns()

                @lock lk begin 
                    step!(lj, dt = dt[])
                    vrescale!(lj, dt = dt[])
                    handleInteractions!(lj, interactions, interactionIndex[], mousepos[], mousestrength[], dt[])
                end
                
                u = time_ns()
                while time_ns() - u < 2e4 end
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

    function renderfunc(x)
        node[] = ensureNeighbors!(lock(() -> deepcopy(lj), lk), forced = false)

        updateInteractions!(interactions, keymap, main_axis)
        mousepos[] = SVec2(mouseposition(main_axis))
        
        if time() - lastaction[] > 0.05
            if interactions[:spawn] && is_mouseinside(main_axis)
                p = mousepos[]
                lock(lk) do 
                    if length(lj) < maxN && minimum(norm.(Ref(p) .- lj.ps)) > lj.σ
                        s = particlemenu.selection[]
                        if s == "default"
                            push!(lj, ps = [p])
                        elseif s == "positive"
                            push!(lj, ps = [p], cs = [1], σs = [2])
                        elseif s == "negative"
                            push!(lj, ps = [p], cs = [-1], σs = [2])
                        elseif s == "heavy" && minimum(norm.(Ref(p) .- lj.ps)) > sqrt(10) * lj.σ
                            push!(lj, ps = [p], ms = [100], σs = [sqrt(10)])
                        elseif s == "very heavy" && minimum(norm.(Ref(p) .- lj.ps)) > sqrt(10) * lj.σ
                            push!(lj, ps = [p], ms = [10000], σs = [sqrt(10)])
                        elseif s == "polar pair" && length(lj) < maxN - 1
                            d = normalize(randn(SVec2)) * lj.σ / 2
                            push!(lj, ps = [p + d, p - d], cs = [-0.3, 0.3], ms = [2, 2], bonds = [(1, 2, 10000.0, 0.05)])
                        elseif s == "neutral pair" && length(lj) < maxN - 1
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
    
    screen::GLMakie.Screen{GLMakie.GLFW.Window} = GLMakie.Screen(renderloop = GLMakie.renderloop, framerate = 120, decorated = decorated)
    on(renderfunc, screen.render_tick)

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

    fig, main_plot, node, screen, runfunc, renderfunc
end