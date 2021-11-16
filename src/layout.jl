using LinearAlgebra
using Graphs
using OMEinsum

@enum SCGate begin 
    Single
    FSim
    Id
    Noise
end
struct RQCLayout{VT, PT, GT}
    g::SimpleGraph{VT}
    edge_patterns::Dict{Edge{VT}, PT}
    gates::Dict{VT, Vector{GT}}
    pattern_loop::Vector{PT}
end

function RQCLayout{VT, PT, GT}(nbits::VT, edge_pattern::Dict{Edge{VT}, PT}, 
        pattern_loop::Vector{PT}) where {VT, PT, GT}
    g = SimpleGraph{VT}(nbits)
    for (e, et) in edge_pattern
        add_edge!(g, e)
    end
    gates = Dict(v => GT[] for v in vertices(g))
    return RQCLayout{VT, PT, GT}(g, edge_pattern, gates, pattern_loop)
end

Base.deepcopy(layout::RQCLayout{VT, PT, GT}) where {VT, PT, GT} = 
    RQCLayout{VT, PT, GT}(deepcopy(layout.g), deepcopy(layout.edge_patterns), 
        deepcopy(layout.gates), deepcopy(layout.pattern_loop))
Graphs.vertices(layout::RQCLayout) = vertices(layout.g)
Graphs.has_vertex(layout::RQCLayout, v) = has_vertex(layout.g, v)
Graphs.neighbors(layout::RQCLayout, v) = neighbors(layout.g, v)
function Graphs.connected_components(g::RQCLayout{VT, PT, SCGate}) where {VT, PT}
    vs_remain = Set(collect(vertices(g)))
    comps = []
    while !isempty(vs_remain)
        v0 = pop!(vs_remain)
        comp_v0 = [v0]
        comp_remain = [v0]
        while !isempty(comp_remain)
            v = pop!(comp_remain)
            gates_v = gates(g, v)
            for i = 1:length(gates_v)
                gates_v[i] === FSim || continue
                u = partner(g, v, depth_to_cycle(i))
                !(u in comp_v0) && push!(comp_remain, u)
                !(u in comp_v0) && push!(comp_v0, u)
                delete!(vs_remain, u)
            end
        end
        push!(comps, sort!(comp_v0))
    end
    return sort!(comps; by = x -> first(x))
end

gates(layout::RQCLayout, v) = layout.gates[v]

function google_layout_53(nbits::VT = 53, ncycles::Integer = 20) where {VT}
    edge_patterns = Dict{Edge{VT}, Symbol}()
    pattern_A = ((31,32),(21,22),(8,11),(24,29),(7,18),(1,5),(26,40),(15,25),
        (6,16),(44,53),(42,48),(46,51),(4,14),(13,27),(28,39),(2,3),(9,17),
        (23,30),(10,12),(19,20),(33,34),(41,47),(43,50),(45,49))
    pattern_B = ((32,37),(21,24),(18,26),(25,44),(22,35),(7,8),(5,15),(16,42),
        (1,4),(2,6),(12,51), (14,36),(3,13),(9,10),(20,41),(27,38),(17,28),
        (19,23),(34,43))
    pattern_C = ((32,52),(24,31),(26,29),(40,44),(22,37),(7,21),(15,18),(25,42),
        (11,35),(1,8),(5,6),(16,51),(3,4),(2,10),(12,41),(27,36),(13,17),(9,19),
        (20,43),(38,39),(28,30),(23,33),(34,45))
    pattern_D = ((21,32),(18,24),(25,26),(44,48),(8,22),(5,7),(15,16),(42,46),
        (4,11),(1,2),(6,12),(47,51),(13,14),(3,9),(10,20),(41,50),(27,28),
        (17,23),(19,34),(43,49))
    pattern_loop = [:A, :B, :C, :D, :C, :D, :A, :B]
    
    for e in pattern_A
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :A
    end
    for e in pattern_B
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :B
    end
    for e in pattern_C
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :C
    end
    for e in pattern_D
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :D
    end
    
    layout = RQCLayout{VT, Symbol, SCGate}(nbits, edge_patterns, pattern_loop)
    if ncycles >= 1
        for i = 1:ncycles
            for v in vertices(layout)
                push!(gates(layout, v), Single)
                if has_partner(layout, v, i)
                    push!(gates(layout, v), FSim)
                else
                    push!(gates(layout, v), Id)
                end
            end
        end
    end
    for v in vertices(layout)
        push!(gates(layout, v), Single)
    end
    return layout
end

cycle_pattern(layout::RQCLayout, ncycles) = getindex(layout.pattern_loop, rem(ncycles-1, length(layout.pattern_loop))+1)
edge_pattern(layout::RQCLayout, u, v) = layout.edge_patterns[Edge(min(u, v), max(u, v))]
function partner(layout::RQCLayout, v, ncycles)
    pt = cycle_pattern(layout, ncycles)
    nbs = neighbors(layout, v)
    idx = findfirst(u -> edge_pattern(layout, u, v) === pt, nbs)
    idx !== nothing && return nbs[idx]
    return nothing
end
has_partner(g, v, ncycles) = (partner(g, v, ncycles) !== nothing)
depth_to_cycle(d) = d÷2
