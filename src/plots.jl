using Compose

function plot(layout::RQCLayout{VT, PT, SCGate}, vs = collect(vertices(layout))) where {VT, PT}
    vs = sort!(intersect(vs, collect(vertices(layout))))
    qubit_map = gen_qubit_map(layout, vs)
    gate_locs = gen_gate_locs(layout, vs)
    ht = length(vs)
    wd = maximum([maximum(gate_locs[v]) for v in vs]) + 1
    set_default_graphic_size((wd+1)*cm, (ht+1)*cm)
    ctx_gates = Compose.compose(Compose.context())
    ctx_lines = Compose.compose(Compose.context())
    for v in vs
        ctx_lines = Compose.compose(ctx_lines, 
            (Compose.context(), line([(0, qubit_map[v]), (wd, qubit_map[v])]), stroke("black"), linewidth(2))
        )
        for d in 1:depth(layout)
            x1 = gate_locs[v][d]
            y1 = qubit_map[v]
            if gates(layout, v)[d] === Noise
                ctx_gates = Compose.compose(ctx_gates, plot_noise(x1, y1))
            elseif gates(layout, v)[d] === Single
                ctx_gates = Compose.compose(ctx_gates, plot_single(x1, y1))
            elseif gates(layout, v)[d] === FSim
                u = partner(layout, v, depth_to_cycle(d))
                x2 = gate_locs[u][d]
                y2 = qubit_map[u]
                y1 < y2 || continue
                ctx_gates = Compose.compose(ctx_gates, plot_fsim(x1, y1, x2, y2))
            end
        end
    end
    ctx = Compose.compose(
        Compose.context(units = UnitBox(-0.5, -0.5, wd+1, ht+1)), 
        ctx_gates, 
        ctx_lines)
    return ctx
end

function gen_qubit_map(layout::RQCLayout{VT, PT, SCGate}, vs = collect(vertices(layout))) where {VT, PT}
    vs = sort!(intersect(vs, collect(vertices(layout))))
    qubit_map = Dict(vs[i] => i for i = 1:length(vs))
    return qubit_map
end

function gen_gate_locs(layout::RQCLayout{VT, PT, SCGate}, vs = collect(vertices(layout))) where {VT, PT}
    qubit_map = gen_qubit_map(layout, vs)
    vs = sort!(intersect(vs, collect(vertices(layout))), by = v -> qubit_map[v])
    gate_locs = Dict(v => Float64[] for v in vs)
    frontier = Dict(v => 1 for v in vs)
    for d = 1:depth(layout)
        for v in vs
            if gates(layout, v)[d] in (Single, Id, Noise)
                push!(gate_locs[v], frontier[v])
                frontier[v] += 1
            else
                u = partner(layout, v, depth_to_cycle(d))
                qubit_map[u] > qubit_map[v] || continue
                loc = max(frontier[u], frontier[v])
                push!(gate_locs[u], loc)
                push!(gate_locs[v], loc)
                for w in vs
                    if qubit_map[v] <= qubit_map[w] <= qubit_map[u]
                        frontier[w] = loc + 1
                    end
                end
            end
        end
    end
    return gate_locs
end

plot_noise(x, y) = (Compose.context(), xgon(x, y, 0.25, 4), 
    fill("red"), stroke("red"))
plot_single(x, y) = (Compose.context(), rectangle(x-0.2, y-0.2, 0.4, 0.4), 
    fill("black"), stroke("black"))
plot_fsim(x1, y1, x2, y2) = (Compose.context(),
    (Compose.context(), circle(x1, y1, 0.2), fill("black"), stroke("black")),
    (Compose.context(), circle(x2, y2, 0.2), fill("black"), stroke("black")),
    (Compose.context(), line([(x1, y1), (x2, y2)]), stroke("black"), linewidth(2))
)