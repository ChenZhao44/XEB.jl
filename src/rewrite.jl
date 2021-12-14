function add_noise!(layout::RQCLayout{VT, PT, SCGate}, v, d) where {VT, PT}
    if gates(layout, v)[d] !== FSim
        gates(layout, v)[d] = Noise
        return true
    end
    return false
end

function update_noise!(layout::RQCLayout{VT, PT, SCGate}; approx = false) where {VT, PT}
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
                        if Noise in (approx ? (gates_u[i], gates_v[2j-i], gates_u[2j-i]) : (gates_u[i], gates_v[2j-i]))
                            gates_v[j] = Noise
                            gates_u[j] = Noise
                            # u == 53 && println("FSim @ $(depth_to_cycle(j)) from $v")
                            # v == 53 && println("FSim @ $(depth_to_cycle(j))")
                            to_update = true
                        end
                    elseif gates_v[j] in (Id, Single)
                        gates_v[j] = Noise
                        # v == 53 && println("Single @ $(depth_to_cycle(j))")
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

simplify!(layout::RQCLayout{VT, PT, SCGate}; approx = false) where {VT, PT} = update_noise!(layout; approx = approx)
function simplify!(layout::RQCLayout{VT, PT, SCGate}, cuts; approx = false) where {VT, PT}
    for c in cuts
        add_noise!(layout, c[1], c[2])
    end
    return simplify!(layout; approx = approx)
end

function add_noise_greedy!(layout::RQCLayout, vs = collect(vertices(layout)), 
    flux_min = 15, out_min = 25)
    update_noise!(layout)
    to_update = true
    while to_update
        to_update = false
        for d in 1:depth(layout), v in vs
            if gates(layout, v)[d] in (Single, Id)
                lo_new = copy(layout)
                gates(lo_new, v)[d] = Noise
                update_noise!(lo_new)
                if flux(lo_new, vs) >= flux_min && 
                        sum(gates(lo_new, v)[end] !== Noise for v in vs) >= out_min
                    add_noise!(layout, v, d)
                    update_noise!(layout)
                    to_update = true
                end
            end
        end
    end
    return layout
end