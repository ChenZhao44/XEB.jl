function tensor_init_xeb(; enable_cuda = false)
    tn = [1.0, 0, 0]
    enable_cuda && return CuArray(tn)
    return tn
end

function tensor_single_xeb(; enable_cuda = false)
    tn = 1/2 * [1.0 1 1; 1 1 -1; 1 -1 1]
    enable_cuda && return CuArray(tn)
    return tn
end

function tensor_measure_xeb(; enable_cuda = false)
    tn = [2.0; 0; 0]
    enable_cuda && return CuArray(tn)
    return tn
end

function tensor_id_xeb(; enable_cuda = false)
    tn = [1.0 0 0; 0 1 0; 0 0 1]
    enable_cuda && return CuArray(tn)
    return tn    
end
function tensor_noise_xeb(; enable_cuda = false)
    tn = [sqrt(1/2), 0, sqrt(1/2)]
    enable_cuda && return CuArray(tn)
    return tn
end

function tensor_fsim_xeb(; enable_cuda = false)
    tn = zeros(3, 3, 3, 3)
    for i = 1:3, j = 1:3
        tn[i, j, j, i] = 1.0
    end
    tn[2, 3, 3, 2] = -sqrt(3)/2
    tn[3, 2, 2, 3] = -sqrt(3)/2
    enable_cuda && return CuArray(tn)
    return tn
end

"""
    to_ein_code_xeb(layout, [s; haar = false, enable_cuda = false])

Generate the tensor network description for `layout` on a subset `s`. If the is a FSim gate out side `s`, 
it will be replaced by an identity gate.
"""
function to_ein_code_xeb(layout::RQCLayout{VT, PT, SCGate}, s = collect(vertices(layout)); haar = false, enable_cuda = false) where {VT, PT}
    layout = deepcopy(layout)
    # update_id!(layout)
    vs = sort!(intersect(s, collect(vertices(layout))))
    nbits = length(vs)
    ids_size = Dict{Int, Int}(i => 3 for i in 1:nbits)
    open_ids_dict = Dict(vs[i] => i for i = 1:nbits)
    tensor_ids = [[id] for id in 1:nbits]
    max_id = nbits
    tensors = Any[]
    for v in vs
        if gates(layout, v)[1] === Noise
            push!(tensors, sqrt(1/2)*tensor_noise_xeb(enable_cuda = enable_cuda))
        else
            push!(tensors, tensor_init_xeb(enable_cuda = enable_cuda))
        end
    end

    D = length(gates(layout, 1))
    for d = 1:D
        for v in vs
            if gates(layout, v)[d] === Single
                ts = tensor_single_xeb(; enable_cuda = enable_cuda)
                max_id += 1
                new_id_v = max_id
                push!(tensors, ts)
                push!(tensor_ids, [open_ids_dict[v], new_id_v])
                ids_size[new_id_v] = 3
                open_ids_dict[v] = new_id_v
            elseif gates(layout, v)[d] === Id
                ts = tensor_id_xeb(; enable_cuda = enable_cuda)
                max_id += 1
                new_id_v = max_id
                push!(tensors, ts)
                push!(tensor_ids, [open_ids_dict[v], new_id_v])
                ids_size[new_id_v] = 3
                open_ids_dict[v] = new_id_v
            elseif gates(layout, v)[d] === Noise
                max_id += 1
                new_id_v = max_id

                push!(tensors, tensor_noise_xeb(enable_cuda = enable_cuda), tensor_noise_xeb(enable_cuda = enable_cuda))
                push!(tensor_ids, [open_ids_dict[v]], [new_id_v])
                ids_size[new_id_v] = 3
                open_ids_dict[v] = new_id_v
            elseif gates(layout, v)[d] === FSim
                c = depth_to_cycle(d)
                u = partner(layout, v, c)
                if u < v && u in vs
                    max_id += 2
                    new_id_u = max_id - 1
                    new_id_v = max_id
                    push!(tensors, tensor_fsim_xeb(enable_cuda = enable_cuda))
                    push!(tensor_ids, [open_ids_dict[u], open_ids_dict[v], new_id_u, new_id_v])
                    ids_size[new_id_u] = 3
                    ids_size[new_id_v] = 3
                    open_ids_dict[u] = new_id_u
                    open_ids_dict[v] = new_id_v
                end
            end
        end
    end
    for v in vs
        push!(tensors, tensor_measure_xeb(enable_cuda = enable_cuda))
        push!(tensor_ids, [open_ids_dict[v]])
        open_ids_dict[v] = 0
    end
    
    open_ids = Int[]
    for v in vs
        open_ids_dict[v] > 0 && push!(open_ids, open_ids_dict[v])
    end
    # return tensor_ids, open_ids
    return EinCode(tensor_ids, open_ids), tensors, ids_size
end
