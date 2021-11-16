function generate_cut(layout::RQCLayout, s, cycle_start, cycle_end)
    ncycles = length(gates(layout, 1)) ÷ 2
    vs = sort!(collect(vertices(layout)))
    s1 = filter!(v -> has_vertex(layout, v), s)
    s2 = setdiff(vs, s)
    cycle_end > ncycles && (cycle_end = ncycles)
    cuts = []
    for c in cycle_start:cycle_end
        for v1 in s1
            v2 = partner(layout, v1, c)
            if v2 in s2
                push!(cuts, (v1, 2*(c-1)+1), (v2, 2*(c-1)+1))
            end
        end
    end
    return cuts
end
