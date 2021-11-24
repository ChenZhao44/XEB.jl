using IterTools

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

function gen_cut_info(layout::RQCLayout, vs_cut)
    vs = collect(vertices(layout))
    cuts_info = Dict{Vector{Int}, Tuple{Int, Int}}()
    ncycle = depth_to_cycle(depth(layout))
    for S in subsets(vs_cut)
        L = setdiff(vs, S)
        layout = google_layout_53(53, ncycle)
        cuts = XEB.generate_cut(layout, L, 1, ncycle)
        XEB.simplify!(layout, cuts);
        
        n_open = length(filter(v -> gates(layout, v)[end] !== XEB.Noise, L))
        cuts_info[S] = (flux(layout, L), n_open)
    end
    return cuts_info
end
