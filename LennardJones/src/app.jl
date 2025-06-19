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

function ssqrt(x)
    sign(x) * sqrt(abs(x))
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
        Legend = (backgroundcolor = :grey20, ),
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
            font_size = 10
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

function primary_resolution()
    monitor = GLMakie.GLFW.GetPrimaryMonitor()
    videomode = GLMakie.MonitorProperties(monitor).videomode
    return (videomode.width, videomode.height)
end

const keymap = (
    :pullsingle => Mouse.right, :pullall => Mouse.middle, :pushall => Keyboard.c,
    :spawn => Mouse.left, :delete => Keyboard.x, 
    :stirleft => Keyboard.q, :stirright => Keyboard.e,
    :pushleft => Keyboard.a, :pushright => Keyboard.d,
    :pushup => Keyboard.w, :pushdown => Keyboard.s,
    :cool => Keyboard.g, :coolall => Keyboard.f, 
    :heat => Keyboard.t,  :heatall => Keyboard.r, 
)

struct Event
    type::Symbol
    position::SVector{2, Float64}
    data::Float64
    
    function Event(type, position = SA[NaN, NaN], data = NaN)
        new(type, position, data)
    end
end

mutable struct InteractiveLJ
    lj::LJ
    outlj::LJ

    interactionIndex::Int
    mouseposition::SVector{2, Float64}
    mousestrength::Float64
    dt::Float64
    interactions::Dict{Symbol, Bool}
    lastaction::Float64
    placetype::Symbol
    preset::Symbol
    
    running::Bool
    terminate::Bool
    events::Channel{Event}
end

function InteractiveLJ(N)
    lj = LJ(N)
    events = Channel{Event}(Inf)
    interactions = Dict(first.(keymap) .=> false)

    InteractiveLJ(lj, deepcopy(lj), 1, SA[0.0, 0.0], 500.0, 1e-4, interactions, time(), :default, :default, false, false, events)
end

function handleInteractions!(sim::InteractiveLJ)
    lj = sim.lj
    dt = sim.dt
    strength = sim.mousestrength
    mousepos = sim.mouseposition
    index = sim.interactionIndex
    interactions = sim.interactions

    if sim.running
        if interactions[:pullsingle] && index <= length(lj)
            lj.vs[index] += strength * dt * (mousepos - lj.ps[index])# / lj.ms[index]
        end
        if interactions[:pullall]
            lj.vs .+= 0.04strength * dt .* (Ref(mousepos) .- lj.ps) ./ (0.05 .+ norm.(Ref(mousepos) .- lj.ps)).^2 ./ lj.ms
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
    if time() - sim.lastaction > 0.05
        if interactions[:spawn]
            p = mousepos
            maxN = 1000
            if length(lj) < maxN && minimum(norm.(Ref(p) .- lj.ps)) > lj.σ
                s = sim.placetype
                if s == :default
                    push!(lj, ps = [p])
                elseif s == :positive
                    push!(lj, ps = [p], cs = [1])
                elseif s == :negative
                    push!(lj, ps = [p], cs = [-1])
                elseif s == :heavy && minimum(norm.(Ref(p) .- lj.ps)) > sqrt(10) * lj.σ
                    push!(lj, ps = [p], ms = [100], σs = [sqrt(10)])
                elseif s == :veryheavy && minimum(norm.(Ref(p) .- lj.ps)) > sqrt(10) * lj.σ
                    push!(lj, ps = [p], ms = [10000], σs = [sqrt(10)])
                elseif s == :polarpair && length(lj) < maxN - 1
                    d = normalize(randn(SVec2)) * lj.σ / 2
                    push!(lj, ps = [p + d, p - d], cs = [-0.3, 0.3], ms = [2, 2], bonds = [(1, 2, 10000.0, 0.05)])
                elseif s == :neutralpair && length(lj) < maxN - 1
                    d = normalize(randn(SVec2)) * lj.σ / 2
                    push!(lj, ps = [p + d, p - d], cs = [0, 0], ms = [2, 2], bonds = [(1, 2, 10000.0, 0.05)])
                end
            end
            sim.lastaction = time()
        end
        if interactions[:delete] && length(lj) > 1
            deleteat!(lj, argmin(norm.(lj.ps .- (mousepos,))))
            sim.lastaction = time()
        end
    end
    ensureNeighbors!(sim.lj)
end

function reset!(sim)
    lj = sim.lj
    N = length(lj)
    if sim.preset == :default
        empty!(lj)
        pushfree!(lj, N)
    elseif sim.preset == :randomneutral
        empty!(lj)
        for i in 1:N 
            m = (rand() + 1) / 2
            pushfree!(lj, 1, ms = [m], σs = [sqrt(m)])
        end
    elseif sim.preset == :randomcharge
        empty!(lj)
        for i in 1:N 
            m = (rand() + 1) / 2
            pushfree!(lj, 1, ms = [m], σs = [sqrt(m)], cs = [2rand()-1])
        end
    elseif sim.preset == :ions
        empty!(lj)
        for i in 1:N 
            pushfree!(lj, cs = [-sign(sum(lj.cs) - 1e-10)]) 
        end
    elseif sim.preset == :diatomic
        empty!(lj)
        pushfree!(lj, N ÷ 2, 
            ps = [SA[lj.bondl * 0.05, 0.0], SA[0.0, 0.0]], 
            bonds = [(1, 2, 10000.0, 0.05)]
        )
    elseif sim.preset == :polar
        empty!(lj)
        pushfree!(lj, N ÷ 2, 
            ps = [SA[lj.bondl * 0.05, 0.0], SA[0.0, 0.0]], 
            cs = [1.0, -1.0], bonds = [(1, 2, 10000.0, 0.05)]
        )
    elseif sim.preset == :chains
        empty!(lj)
        pushfree!(lj, N ÷ 5, 
            ps = [SA[i * lj.bondl * 0.025, 0.0] for i in 1:5], 
            bonds = [(1, 2, 10000.0, 0.025), (2, 3, 10000.0, 0.025), (3, 4, 10000.0, 0.025), (4, 5, 10000.0, 0.025)], 
            angles = [(2, 1, 3, 10, 180), (3, 2, 4, 10, 180), (4, 3, 5, 10, 180)]
        )
    elseif sim.preset == :rings
        empty!(lj)
        pushfree!(lj, N ÷ 6, 
            ps = [lj.bondl * 0.025 * SA[cos(ϕ), sin(ϕ)] for ϕ in 0:2pi/6:6], 
            bonds = [(1, 2, 10000.0, 0.025), (2, 3, 10000.0, 0.025), (3, 4, 10000.0, 0.025), 
                        (4, 5, 10000.0, 0.025), (5, 6, 10000.0, 0.025), (6, 1, 10000.0, 0.025)],
            angles = [(2, 1, 3, 10, 180), (3, 2, 4, 10, 180), (4, 3, 5, 10, 180), 
                        (5, 4, 6, 10, 180), (6, 5, 1, 10, 180), (1, 6, 2, 10, 180)]
        )
    end
    ensureNeighbors!(lj)
end

function start(sim::InteractiveLJ)
    function runfunc(sim)
        lj = sim.lj
        t = time()
        n = 0
        while true
            if sim.terminate
                return
            end
            for i in 1:100
                while isready(sim.events)
                    event = take!(sim.events)

                    if event.type == :changenumber
                        while length(lj) > event.data
                            deleteat!(lj, length(lj))
                        end
                        pushfree!(lj, round(Int, event.data - length(lj)))
                    elseif event.type == :reset
                        reset!(sim)
                    elseif event.type == :setindex
                        sim.interactionIndex = argmin(norm.(sim.lj.ps .- (sim.mouseposition,)))
                    elseif event.type == :freeze
                        sim.lj.vs .*= 0
                    end
                end
                handleInteractions!(sim)
                if sim.running
                    step!(lj, dt = sim.dt)
                    vrescale!(lj, dt = sim.dt)
                end
                n += 1
                if n > 1500 * (time() - t) - 0.25
                    break
                end
            end
            n = 0
            t = time()
            sim.outlj = deepcopy(ensureNeighbors!(lj))
            sleep(1e-3)
        end
    end
    Threads.@spawn runfunc(sim)
end

function stop(sim::InteractiveLJ)
    sim.terminate = true
end

function locale(lang)
    if lang == :de
        (; 
            particles = "Teilchen",
            number = "Anzahl",
            size = "Größe",
            thermostat = "Thermostat",
            temperature = "Zieltemperatur",
            tstrength = "Kopplungsstärke",
            forces = "Kräfte",
            coulomb = "Coulomb",
            lennardjones = "Lennard Jones",
            bondstrength = "Bindungsstärke",
            bondlength = "Bindungslänge",
            gravity = "Schwerkraft",
            voltage = "Elektrisches Feld",
            interactive = "Interaktiv",
            timestep = "Zeitschritt",
            start = "Start",
            stop = "Stop",
            freeze = "Einfrieren",
            placedefault = "Platzieren: Standard",
            placeheavy = "Platzieren: Schwer",
            placeveryheavy = "Platzieren: Sehr schwer",
            placepostive = "Platzieren: Positiv",
            placenegative = "Platzieren: Negativ",
            placepolarpair = "Platzieren: Polares Paar",
            placeneutralpair = "Platzieren: Neutrales Paar",
            colorpotential = "Farbe: Potential",
            colorvelocity = "Farbe: Geschwindigkeit",
            colorvirial = "Farbe: Druck",
            colorcharge = "Farbe: Ladung",
            colornothing = "Farbe: Keine",
            presetdefault = "Preset: Standard",
            presetrandomneutral = "Preset: Zufällige Größe",
            presetrandomcharge = "Preset: Zufällige Ladung",
            presetions = "Preset: Ionen",
            presetdiatomic = "Preset: Zweiatomig",
            presetpolar = "Preset: Polar",
            presetchains = "Preset: Ketten",
            presetrings = "Preset: Ringe",
            reset = "Reset",
            kinetic = "kinetische Energie",
            potential = "potentielle Energie",
            pressure = "Druck",
            controlstext = rich(
                rich("Linksklick: Neues Teilchen hinzufügen, Mittelklick: alle Teilchen bewegen, Rechtsklick: einzelnes Teilchen bewegen\n", fontsize = 20),
                rich("WASD: Kraft nach oben/links/unten/rechts, Q/E: links/rechts drehen, R/F: alle erhitzen/kühlen, T/G: lokal erhitzen/kühlen, X: löschen, C: Teilchen", fontsize = 15)
            ),
        )
    else
        (; 
            particles = "Particle",
            number = "number",
            size = "size",
            thermostat = "Thermostat",
            temperature = "target temperature",
            tstrength = "coupling strength",
            forces = "Forces",
            coulomb = "Coulomb",
            lennardjones = "Lennard Jones",
            bondstrength = "bond strength",
            bondlength = "bond length",
            gravity = "gravity",
            voltage = "voltage",
            interactive = "interactive",
            timestep = "time step",
            start = "Start",
            stop = "Stop",
            freeze = "Freeze",
            placedefault = "Place: default",
            placeheavy = "Place: heavy", 
            placeveryheavy = "Place: very heavy",
            placepostive = "Place: positive",
            placenegative = "Place: negative",
            placepolarpair = "Place: polar pair",
            placeneutralpair = "Place: neutral pair",
            colorpotential = "Color: potential",
            colorvelocity = "Color: velocity",
            colorvirial = "Color: virial",
            colorcharge = "Color: charge",
            colornothing = "Color: nothing",
            presetdefault = "Preset: default",
            presetrandomneutral = "Preset: random neutral",
            presetrandomcharge = "Preset: random charge",
            presetions = "Preset: ions",
            presetdiatomic = "Preset: diatomic",
            presetpolar = "Preset: polar", 
            presetchains = "Preset: chains", 
            presetrings = "Preset: rings",
            reset = "Reset",
            kinetic = "kinetic energy",
            potential = "potential energy",
            pressure = "pressure",
            controlstext = rich(
                rich("Left click: add new particle, Middle click: move all particles, Right click: move single particle\n", fontsize = 20),
                rich("WASD: force up/left/down/right, Q/E: rotate left/right, R/F: heat/cool all, T/G: heat/cool locally, X: delete particle, C: repel particles", fontsize = 15)
            ),
        )
    end
end

function makesettings(gridpos, sim, loc)
    settings = GridLayout(gridpos)
    sliderconf = (
        (label = loc.number, range = 1:1:1000, startvalue = length(sim.lj)) => (val -> put!(sim.events, Event(:changenumber, SA[NaN, NaN], 1.0val))),
        (label = loc.size, range = 0.01:0.0001:0.2, startvalue = sim.lj.σ) => (val -> (sim.lj.σ = val)),
        (label = loc.temperature, range = (0:0.001:2).^log2(10), startvalue = sim.lj.T) => val -> (sim.lj.T = val),
        (label = loc.tstrength, range = (0:0.01:10).^3, startvalue = sim.lj.ts) => (val -> (sim.lj.ts = val)),
        (label = loc.coulomb, range = (0:0.001:10).^2, startvalue = sim.lj.coulomb) => (val -> (sim.lj.coulomb = val)),
        (label = loc.lennardjones, range = 0:0.01:10, startvalue = sim.lj.ε) => (val -> (sim.lj.ε = val)),
        # (label = "Stillinger-Weber", range = 0:0.01:10, startvalue = sim.lj.swstrength) => (val -> (sim.lj.swstrength = val)),
        # (label = "covalent angle", range = 0:180, startvalue = sim.lj.swangle) => (val -> (sim.lj.swangle = val)),
        (label = loc.bondstrength, range = (0:0.01:10), startvalue = sim.lj.bondk) => (val -> (sim.lj.bondk = val)),
        (label = loc.bondlength, range = (0:0.01:2), startvalue = sim.lj.bondl) => (val -> (sim.lj.bondl = val)),
        (label = loc.gravity, range = 0:-0.01:-10, startvalue = sim.lj.g) => (val -> (sim.lj.g = val)),
        (label = loc.voltage, range = (0:0.1:100).^3, startvalue = sim.lj.g) => (val -> (sim.lj.E = val)),
        # (label = "magnetic", range = (0:0.1:100), startvalue = sim.lj.g) => (val -> (sim.lj.B = val)),
        (label = loc.interactive, range = 1:1:1000, startvalue = sim.mousestrength) => (val -> (sim.mousestrength = val)),
        (label = loc.timestep, range = logrange(1e-6, 1e-3, length = 1000), startvalue = sim.dt) => (val -> (sim.dt = val)),
    )
    on.(last.(sliderconf), getfield.(SliderGrid(settings[1,1], first.(sliderconf)...).sliders, :value))
    rowgap!(contents(settings[1,1])[1].layout, 7)
    for gap in [2, 4, 6, 8, 10, 11]
        rowgap!(contents(settings[1,1])[1].layout, gap, 21)
    end
end

function main(; lang = :en, history = true)
    sim = InteractiveLJ(500)
    node = Observable(sim.outlj)
    
    loc = locale(lang)

    fig = Figure(size = (1200, 700), fontsize = 20)

    makesettings(fig[1,2][1,1], sim, loc)

    main_axis::Axis = Axis(fig[1,1], aspect = DataAspect())
    xlims!(main_axis, -1.03, 1.03); ylims!(main_axis, -1.03, 1.03); hidedecorations!(main_axis)
    remove_interactions!(main_axis)
    colsize!(fig.layout, 1, Aspect(1, 1.0))

    color = Observable(zeros(1000))
    colorrange = Observable((0.0, 1.0))
    positions = Observable(fill(SA[0.0,0.0], 1000))
    markersize = Observable(fill(node[].σ, 1000))
    charge = Observable(fill(0, 1000))
    colorby = Observable(:potential)
    
    function plotfunc(lj)
        N = length(lj)
        positions[][1:N] .= lj.ps
        positions[][N+1:end] .= Ref(SA[NaN, NaN])

        charge[][1:N] .= round.(Int, sign.(lj.cs))
        charge[][N+1:end] .= 0
        notify(charge)

        if colorby[] == :charge
            color[][1:N] .= lj.cs
        elseif colorby[] == :velocity
            color[][1:N] .= 0.5 .* color[][1:N] .+ 0.5 .* norm.(lj.vs)
        elseif colorby[] == :potential
            cs = potentialPerParticle(lj) ./ lj.σs
            color[][1:N] .= sqrt.(cs .- minimum(cs))
        elseif colorby[] == :virial
            color[][1:N] .= ssqrt.(virialPerParticle(lj))
        elseif colorby[] == :nothing
            color[][1:N] .= 0
        end
        colorrange[] = (0.8colorrange[][1] + 0.2minimum(color[][1:N]), 0.8colorrange[][2] + 0.2maximum(color[][1:N]))
        notify(color)

        if any(markersize[][i] != lj.σ * lj.σs[i] for i in 1:N)
            markersize[][1:N] .= lj.σ .* lj.σs
            notify(markersize)
        end
        notify(positions)
    end
    on(plotfunc, node) 
    node[] = sim.outlj
    colorrange[] = extrema(color[])

    linesegments!(main_axis, 
        lift(lj -> lj.ps[[b[i] for b in lj.bonds for i in 1:2]], node), 
        linewidth = lift(lj -> 50lj.σ, node), color = cgrad(:isoluminant_cm_70_c39_n256)[0.5]
    )
    scatter!(main_axis, positions, 
        markersize = markersize, strokewidth = 0, color = color, 
        markerspace = :data, colormap = :isoluminant_cm_70_c39_n256, 
        colorrange = colorrange
    )
    scatter!(main_axis, positions, 
        markersize = lift((ms, cs) -> ms .* ifelse.(cs .== 0, 0, ifelse.(cs .> 0, 1.1, 1.5)), markersize, charge), color = :grey10, 
        markerspace = :data, marker = lift(cs -> getindex.((('-', 'n', '+'),), cs .+ 2), charge)
    )
    scatter!(main_axis, 
        lift(lj -> (i = clamp(sim.interactionIndex, 1, length(lj)); lj.ps[i:i]), node), 
        markersize = lift(lj -> 1.05lj.σ, node), markerspace = :data, 
        color = lift(lj -> ifelse(sim.interactions[:pullsingle], :red, :transparent), node)
    )

    menufontsize = 25
    menugrid = GridLayout(fig[1,2][2,1], tellwidth = false, tellheight = true, halign = :left)
    startbutton = Button(menugrid[1,1:2], label = loc.start, width = Relative(1), fontsize = menufontsize)
    on(Button(menugrid[1,3:4], label = loc.freeze, width = Relative(1), fontsize = menufontsize).clicks) do _
        put!(sim.events, Event(:freeze))
    end

    particlemenu = Menu(menugrid[2,1:2], width = Relative(1), fontsize = menufontsize,
        options = [
            (loc.placedefault, :default),
            (loc.placeheavy, :heavy),
            (loc.placeveryheavy, :veryheavy),
            (loc.placepostive, :positive),
            (loc.placenegative, :negative),
            (loc.placepolarpair, :polarpair),
            (loc.placeneutralpair, :neutralpair)
        ]
    )
    on(x -> sim.placetype = x, particlemenu.selection)

    colormenu = Menu(menugrid[2,3:4], width = Relative(1), fontsize = menufontsize, 
        options = [
            (loc.colorpotential, :potential), 
            (loc.colorvelocity, :velocity),
            (loc.colorvirial, :virial), 
            (loc.colorcharge, :charge), 
            (loc.colornothing, :nothing)
        ]
    )
    on(x -> colorby[] = x, colormenu.selection)

    presetmenu::Menu = Menu(menugrid[2, 5:6], width = Relative(1), fontsize = menufontsize,
        options = [
            (loc.presetdefault, :default),
            (loc.presetrandomneutral, :randomneutral),
            (loc.presetrandomcharge, :randomcharge),
            (loc.presetions, :ions),
            (loc.presetdiatomic, :diatomic),
            (loc.presetpolar, :polar),
            (loc.presetchains, :chains),
            (loc.presetrings, :rings)
        ]
    )
    on(presetmenu.selection) do sel
        sim.preset = sel
        put!(sim.events, Event(:reset))
    end

    on(Button(menugrid[1,5:6], label = loc.reset, width = Relative(1), fontsize = menufontsize).clicks) do _
        put!(sim.events, Event(:reset))
    end


    sps = Observable(0.0); fps = Observable(0); lastnsteps = Observable(0); nframes = Observable(0); lastframe = Observable(time_ns())
    Label(menugrid[3,5:6], halign = :right, 
        text = updateevery(lift((s, f, lj) -> "tps: $(round(Int, s))\nfps: $f\nN: $(length(lj))", sps, fps, node), 0.2), 
        tellwidth = false, justification = :right, color = :grey60, fontsize = 12
    )
    on(node) do lj
        nframes[] += 1
        if time_ns() - lastframe[] > 1e9
            sps[] = (lj.nsteps - lastnsteps[]) / (time_ns() - lastframe[]) * 1e9
            lastnsteps[] = lj.nsteps; fps[] = nframes[]; nframes[] = 0; lastframe[] = time_ns()
        end
    end

    Label(menugrid[3,1:5], loc.controlstext, justification = :left, halign = :left, color = :grey60, padding = (0,0,0,0))
    rowgap!(menugrid, 1, 15)
        
    if history
        plotgrid = GridLayout(fig[1,2][3,1], alignmode = Outside())
        rowgap!(plotgrid, 5)
        mcolor = cgrad(:isoluminant_cm_70_c39_n256)[1.0]

        histsteps = 1001
        axis = Axis(plotgrid[1,1], ylabelsize = 20, xticklabelsvisible = false)
        xlims!(axis, -1.02histsteps, 0.02histsteps)
        
        temperatures = Observable(fill(NaN, histsteps))
        lines!(axis, -histsteps+1:1:0, temperatures, color = Makie.wong_colors()[1], label = loc.kinetic)
        potentials = Observable(fill(NaN, histsteps))
        lines!(axis, -histsteps+1:1:0, potentials, color = Makie.wong_colors()[2], label = loc.potential)
        pressures = Observable(fill(NaN, histsteps))
        lines!(axis, -histsteps+1:1:0, pressures, color = Makie.wong_colors()[3], label = loc.pressure)
        Legend(plotgrid[2,1], axis, position = :lt, orientation = :horizontal, framevisible = false, tellwidth = false, fontsize = 20)

        function updatePlots(lj)
            T = temperature(lj)
            if T != temperatures[][end]
                temperatures[][1:end-1] .= temperatures[][2:end]
                temperatures[][end] = T
                notify(temperatures)
            end
            V = potential(lj) + 3
            if V != potentials[][end] 
                potentials[][1:end-1] .= potentials[][2:end]
                potentials[][end] = V
                notify(potentials)
            end
            P = T / 4 + virial(lj) / 8 / length(lj)
            if P != pressures[][end]
                pressures[][1:end-1] .= pressures[][2:end]
                pressures[][end] = P
                notify(pressures)
            end
            
            reset_limits!(axis)
        end    
        on(updatePlots, node)
        
        remove_interactions!(axis)
        onmouseleftdown(axis.mouseeventhandle) do _
            potentials[] .= NaN; 
            temperatures[] .= NaN; 
            pressures[] .= NaN; 
        end
    else
        plotgrid = GridLayout(fig[1,2][3,1], alignmode = Outside())
    end


    on(events(fig).mouseposition, priority = 1) do event
        sim.mouseposition = SVector{2, Float64}(mouseposition(main_axis))
    end

    on(events(fig).mousebutton, priority = 10) do event
       if event.button == Mouse.right && event.action == Mouse.press && is_mouseinside(main_axis)
            put!(sim.events, Event(:setindex))
        end
        map(keymap) do (name, key)
            if key isa Mouse.Button && (is_mouseinside(main_axis) || sim.interactions[name])
                sim.interactions[name] = ispressed(fig, key)
            end
        end
        return Consume(is_mouseinside(main_axis))
    end

    on(events(fig).keyboardbutton, priority = 1) do event
        map(keymap) do (name, key)
            if key isa Keyboard.Button
                sim.interactions[name] = ispressed(main_axis, key)
            end
        end
    end


    function renderfunc(x)
        node[] = sim.outlj
    end

    on(startbutton.clicks) do n
        sim.running = !sim.running
        startbutton.label = sim.running ? loc.stop : loc.start
    end
    screen::GLMakie.Screen{GLMakie.GLFW.Window} = GLMakie.Screen(renderloop = GLMakie.renderloop, framerate = 60.0, vsync = false, render_on_demand = true)
    on(renderfunc, events(fig).tick)

    on(screen.window_open) do val
        if !val
            stop(sim)
        end
    end

    display(screen, fig)
    GLMakie.GLFW.SetWindowPos(screen.glscreen, 0, 0)
    start(sim)
    
    screen
end

