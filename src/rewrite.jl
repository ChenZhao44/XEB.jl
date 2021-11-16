function add_noise!(layout::RQCLayout{VT, PT, Symbol}, v, d) where {VT, PT}
    if gates(layout, v)[d] !== :fsim
        gates(layout, v)[d] = :noise
        return true
    end
    return false
end
function update_noise!(layout::RQCLayout{VT, PT, Symbol}) where {VT, PT}
    to_update = true
    while to_update
        to_update = false
        for v in vertices(layout)
            gates_v = gates(layout, v)
            for i in 1:length(gates_v)
                gates_v[i] === :noise || continue
                for j in (i-1, i+1)
                    1 <= j <= length(gates_v) || continue
                    if gates_v[j] === :fsim
                        u = partner(layout, v, depth_to_cycle(j))
                        gates_u = gates(layout, u)
                        if :noise in (gates_u[i], gates_u[2j-i], gates_v[2j-i])
                            gates_v[j] = :noise
                            gates_u[j] = :noise
                            to_update = true
                        end
                    elseif gates_v[j] in (:id, :single)
                        gates_v[j] = :noise
                        to_update = true
                    end
                end
            end
        end
    end
    return layout
end
function update_id!(layout::RQCLayout{VT, PT, Symbol}) where {VT, PT}
    for v in vertices(layout)
        gates_v = gates(layout, v)
        for i = 1:(length(gates_v)-1)
            if gates_v[i] === :noise && gates_v[i+1] === :noise
                i > 1 && (gates_v[i] = :id)
            end
        end
    end
    for v in vertices(layout)
        gates_v = gates(layout, v)
        if gates_v[1] === :noise && gates_v[2] === :noise
            gates_v[2] = :id
        end
    end
    return layout
end
simplify!(layout::RQCLayout{VT, PT, Symbol}) where {VT, PT} = update_id!(update_noise!(layout))
function simplify!(layout::RQCLayout{VT, PT, Symbol}, cuts) where {VT, PT}
    for c in cuts
        add_noise!(layout, c[1], c[2])
    end
    return simplify!(layout)
end

function my_cut_53(cycle_start, cycle_end, nbits = 53, ncycles = 20)
    layout = google_layout_53(nbits, ncycles)
    cycle_end > ncycles && (cycle_end = ncycles)
    vs_cut = [1, 2, 4, 6, 11, 12, 47, 51]
    filter!(x -> has_vertex(layout, x), vs_cut)
    cuts = []
    for c in cycle_start:cycle_end
        for v1 in vs_cut
            v2 = partner(layout, v1, c)
            if v2 in vs_cut && v1 < v2
                push!(cuts, (v1, 2*(c-1)+1), (v2, 2*(c-1)+1))
            end
        end
    end
    return cuts
end
