function add_noise!(g, v, d) 
    if gates(g, v)[d] !== :fsim
        gates(g, v)[d] = :noise
        return true
    end
    return false
end
function update_noise!(g)
    to_update = true
    while to_update
        to_update = false
        for v in vertices(g)
            gates_v = gates(g, v)
            for i in 1:length(gates_v)
                gates_v[i] === :noise || continue
                for j in (i-1, i+1)
                    1 <= j <= length(gates_v) || continue
                    if gates_v[j] === :fsim
                        ncycle = j÷2
                        u = partner(g, v, ncycle)
                        gates_u = gates(g, u)
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
    return g
end
function update_id!(g)
    for v in vertices(g)
        gates_v = gates(g, v)
        for i = 1:(length(gates_v)-1)
            if gates_v[i] === :noise && gates_v[i+1] === :noise
                i > 1 && (gates_v[i] = :id)
            end
        end
    end
    for v in vertices(g)
        gates_v = gates(g, v)
        if gates_v[1] === :noise && gates_v[2] === :noise
            gates_v[2] = :id
        end
    end
    return g
end
simplify!(g) = update_id!(update_noise!(g))
function simplify!(g, cuts)
    for c in cuts
        add_noise!(g, c[1], c[2])
    end
    return simplify!(g)
end

function my_cut_53(d_start, d_end, nbits = 53, depth = 20)
    g = google_layout_53(nbits, depth)
    vs_cut = [1, 2, 4, 6, 11, 12, 47, 51]
    filter!(x -> x <= nv(g), vs_cut)
    cuts = []
    for c in d_start:d_end
        for v1 in vs_cut
            v2 = partner(g, v1, c)
            if v2 in vs_cut
                push!(cuts, (v1, 2*(c-1)+1))
                push!(cuts, (v2, 2*(c-1)+1))
            end
        end
    end
    return cuts
end
