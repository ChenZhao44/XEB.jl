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
