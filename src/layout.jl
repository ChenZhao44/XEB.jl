using LinearAlgebra
using Graphs, MetaGraphs
using OMEinsum

const pattern_A = [(31,32),(21,22),(8,11),(24,29),(7,18),(1,5),(26,40),(15,25),
    (6,16),(44,53),(42,48),(46,51),(4,14),(13,27),(28,39),(2,3),(9,17),
    (23,30),(10,12),(19,20),(33,34),(41,47),(43,50),(45,49)]
const pattern_B = [(32,37),(21,24),(18,26),(25,44),(22,35),(7,8),(5,15),(16,42),
    (1,4),(2,6),(12,51), (14,36),(3,13),(9,10),(20,41),(27,38),(17,28),
    (19,23),(34,43)]    
const pattern_C = [(32,52),(24,31),(26,29),(40,44),(22,37),(7,21),(15,18),(25,42),
    (11,35),(1,8),(5,6),(16,51),(3,4),(2,10),(12,41),(27,36),(13,17),(9,19),
    (20,43),(38,39),(28,30),(23,33),(34,45)]
const pattern_D = [(21,32),(18,24),(25,26),(44,48),(8,22),(5,7),(15,16),(42,46),
    (4,11),(1,2),(6,12),(47,51),(13,14),(3,9),(10,20),(41,50),(27,28),
    (17,23),(19,34),(43,49)]

function google_layout_53(nbits = 53, ncycle = 20)
    g = MetaGraph(nbits)
    for (es, et) in ((pattern_A, :A), (pattern_B, :B), (pattern_C, :C), (pattern_D, :D))
        for e in es
            add_edge!(g, e) && set_prop!(g, e[1], e[2], :pattern, et)
        end
    end
    for v in vertices(g)
        set_prop!(g, v, :gates, [])
    end
    if ncycle >= 1
        for i = 1:ncycle-1
            for v in vertices(g)
                push!(get_prop(g, v, :gates), :single)
                if has_partner(g, v, i)
                    push!(get_prop(g, v, :gates), :fsim)
                else
                    push!(get_prop(g, v, :gates), :id)
                end
            end
        end
        for v in vertices(g)
            push!(get_prop(g, v, :gates), :single)
        end
    end
    return g
end

pattern_type(ncycle) = getindex((:A, :B, :C, :D, :C, :D, :A, :B), rem(ncycle-1,8)+1)
function partner(g, v, ncycle)
    pt = pattern_type(ncycle)
    nbs = neighbors(g, v)
    idx = findfirst(u -> get_prop(g, v, u, :pattern) === pt, nbs)
    idx !== nothing && return nbs[idx]
    return nothing
end
has_partner(g, v, ncycle) = (partner(g, v, ncycle) !== nothing)
gates(g, v) = get_prop(g, v, :gates)
