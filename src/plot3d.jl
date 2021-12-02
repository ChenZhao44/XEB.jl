using GLMakie
export plot3d, plot3d_space_bottlenecks!

function plot3d(l::RQCLayout, vs = vertices(l))
    vs = sort!(intersect(vs, vertices(l)))
    D = depth(l)
    scene = meshscatter([l.qubit_locs[v][1] for v in vs], 
        [l.qubit_locs[v][2] for v in vs], 
        zeros(length(vs)) .- 0.01 * rand(length(vs)),
        color = :blue,
        markersize = 0.15
    )
    for v in vs
        loc = l.qubit_locs[v]
        text!("$v", position = (loc[1], loc[2], -0.5))
        d1 = 0
        d2 = 1
        while d1 < D
            if d2 <= D && gates(l, v)[d2] !== Noise
                d2 += 1
            else
                lines!([loc[1], loc[1]], [loc[2], loc[2]], [0.5d1, 0.5d2])
                d1 = d2 + 1
                while d1 < D && gates(l, v)[d1+1] === Noise
                    d1 += 1
                end
                d2 = d1 + 1
                if d1 === D && gates(l, v)[d1] === Noise
                    lines!([loc[1], loc[1]], [loc[2], loc[2]], [0.5d1, 0.5d2])
                end
            end
        end
    end
    
    for v1 in vs, d = 1:D
        if gates(l, v1)[d] === FSim
            v2 = partner(l, v1, depth_to_cycle(d))
            v2 in vs || continue
            v1 > v2 && continue
            loc1 = l.qubit_locs[v1]
            loc2 = l.qubit_locs[v2]
            lines!([loc1[1], loc2[1]], [loc1[2], loc2[2]], [d*0.5, d*0.5]; color = :red)
        end
    end
    
    v_xs = Float64[]
    v_ys = Float64[]
    v_zs = Float64[]
    v_sizes = Float64[]
    v_colors = Symbol[]
    for d = 1:D
        for v in vs
            push!(v_xs, l.qubit_locs[v][1])
            push!(v_ys, l.qubit_locs[v][2])
            push!(v_zs, d * 0.5)
            if gates(l, v)[d] === FSim
                push!(v_sizes, 0.1)
                push!(v_colors, :red)
            elseif gates(l, v)[d] === Noise
                push!(v_sizes, 0.1)
                if (d == 1 || d == D) || gates(l, v)[d-1] !== Noise || gates(l, v)[d+1] !== Noise
                    push!(v_colors, :white)
                else
                    push!(v_colors, :transparent)
                end
            elseif gates(l, v)[d] === Single
                push!(v_sizes, 0.05)
                push!(v_colors, :black)
            else
                push!(v_sizes, 0.05)
                push!(v_colors, :transparent)
            end
        end
    end
    meshscatter!(v_xs, v_ys, v_zs; color = v_colors, markersize = v_sizes)
    meshscatter!([l.qubit_locs[v][1] for v in vs], 
        [l.qubit_locs[v][2] for v in vs], 
        ones(length(vs)) * (D+1) * 0.5,
        color = :blue,
        markersize = 0.15
    )
    return scene
end

function plot3d_space_bottlenecks!(layout, ec, ids_dic, ids_loc)
    max_ids = max_space_ids(ec, ids_dic)
    scene = current_figure()
    for i in max_ids
        loc = layout.qubit_locs[ids_loc[i][1]]
        scene = lines!([loc[1], loc[1]], [loc[2], loc[2]], 
            [(ids_loc[i][2])/2, (ids_loc[i][2]+1)/2],
            color = :green,
            linewidth = 10
        )
    end
    return scene
end
