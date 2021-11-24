using LinearAlgebra
using Multigraphs, Graphs
using OMEinsum

@enum SCGate begin 
    Single
    FSim
    Id
    Noise
end
struct RQCLayout{VT, PT, GT}
    g::Multigraph{VT}
    edge_patterns::Dict{Edge{VT}, PT}
    gates::Dict{VT, Vector{GT}}
    pattern_loop::Vector{PT}
end

function RQCLayout{VT, PT, GT}(vs::Vector{VT}, edge_pattern::Dict{Edge{VT}, PT}, 
        pattern_loop::Vector{PT}) where {VT, PT, GT}
    g = Multigraph{VT}(Dict(v => VT[] for v in vs), maximum(vs))
    for (e, et) in edge_pattern
        add_edge!(g, e)
    end
    gates = Dict(v => GT[] for v in vertices(g))
    return RQCLayout{VT, PT, GT}(g, edge_pattern, gates, pattern_loop)
end
RQCLayout{VT, PT, GT}(nbits::VT, edge_pattern::Dict{Edge{VT}, PT}, 
        pattern_loop::Vector{PT}) where {VT, PT, GT} = 
    RQCLayout{VT, PT, GT}(VT[i for i = 1:nbits], edge_pattern, pattern_loop)

Base.copy(layout::RQCLayout) = deepcopy(layout)
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

function depth(g::RQCLayout) 
    v0 = first(keys(g.gates))
    return length(gates(g, v0))
end

function flux(g::RQCLayout{VT, PT, SCGate}, vs = collect(vertices(g))) where {VT, PT}
    M = [gates(g, v)[i] !== Noise for v in vs, i in 1:depth(g)]
    return minimum(sum(M, dims = 1))
end