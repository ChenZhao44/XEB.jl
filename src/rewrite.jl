function add_noise!(layout::RQCLayout{VT, PT, SCGate}, v, d) where {VT, PT}
    if gates(layout, v)[d] !== FSim
        gates(layout, v)[d] = Noise
        return true
    end
    return false
end

function update_noise!(layout::RQCLayout{VT, PT, SCGate}) where {VT, PT}
    to_update = true
    while to_update
        to_update = false
        for v in vertices(layout)
            gates_v = gates(layout, v)
            for i in 1:length(gates_v)
                gates_v[i] === Noise || continue
                for j in (i-1, i+1)
                    1 <= j <= length(gates_v) || continue
                    if gates_v[j] === FSim
                        u = partner(layout, v, depth_to_cycle(j))
                        gates_u = gates(layout, u)
                        if Noise in (gates_u[i], gates_u[2j-i], gates_v[2j-i])
                            gates_v[j] = Noise
                            gates_u[j] = Noise
                            to_update = true
                        end
                    elseif gates_v[j] in (Id, Single)
                        gates_v[j] = Noise
                        to_update = true
                    end
                end
            end
        end
    end
    return layout
end

function update_id!(layout::RQCLayout{VT, PT, SCGate}) where {VT, PT}
    for v in vertices(layout)
        gates_v = gates(layout, v)
        for i = 1:(length(gates_v)-1)
            if gates_v[i] === Noise && gates_v[i+1] === Noise
                i > 1 && (gates_v[i] = Id)
            end
        end
    end
    for v in vertices(layout)
        gates_v = gates(layout, v)
        if gates_v[1] === Noise && gates_v[2] === Noise
            gates_v[2] = Id
        end
    end
    return layout
end

simplify!(layout::RQCLayout{VT, PT, SCGate}) where {VT, PT} = update_id!(update_noise!(layout))
function simplify!(layout::RQCLayout{VT, PT, SCGate}, cuts) where {VT, PT}
    for c in cuts
        add_noise!(layout, c[1], c[2])
    end
    return simplify!(layout)
end
