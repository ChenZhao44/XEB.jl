using RandomMatrices
using CUDA

function tensor_init_pauli(; enable_cuda = false)
    tn = Float64[1, 0, 0, 1] / sqrt(2)
    enable_cuda && return CuArray(tn)
    return tn
end
function tensor_measure_pauli(; enable_cuda = false)
    tn = (1/sqrt(2))*Float64[1 1; 0 0; 0 0; 1 -1]
    enable_cuda && return CuArray(tn)
    return tn
end
function tensor_single_pauli(g_prev = :null; haar = false, enable_cuda = false)
    I = [1 0; 0 1]
    X = [0 1; 1 0]
    Y = [0 -im; im 0]
    Z = [1 0; 0 -1]
    paulis = (I, X, Y, Z)
    if haar
        M = rand(Haar(2), 2)
        g_new = :null
    else    # Haar random matrix
        g_prev === :null && (g_new = rand((:sx, :sy, :sw)))
        g_prev === :sx && (g_new = rand((:sy, :sw)))
        g_prev === :sy && (g_new = rand((:sx, :sw)))
        g_prev === :sw && (g_new = rand((:sx, :sy)))

        g_new === :sx && (M = sqrt(X))
        g_new === :sy && (M = sqrt(Y))
        g_new === :sw && (M = sqrt(sqrt(1/2)*(X+Y)))
        M = [1 0; 0 exp(rand()*2π*im)]*M*[1 0; 0 exp(rand()*2π*im)]
    end
    ts = zeros(Float64, 4, 4)
    for i = 1:4, j = 1:4
        ts[i, j] = real(tr(M*paulis[i]*M'*paulis[j]/2))
    end
    enable_cuda && return CuArray(ts), g_new
    return ts, g_new
end
function tensor_noise_pauli(; enable_cuda = false)
    tn = zeros(Float64, 4)
    tn[1] = 1
    enable_cuda && return CuArray(tn)
    return tn
end
function tensor_fsim_pauli(; enable_cuda = false)
    fsim_gate = zeros(ComplexF64, 4, 4)
    fsim_gate[1, 1] = 1
    fsim_gate[2, 3] = -im
    fsim_gate[3, 2] = -im
    fsim_gate[4, 4] = exp(-pi/6*im)
    tf = zeros(Float64, 4, 4, 4, 4)
    I = [1 0; 0 1]
    X = [0 1; 1 0]
    Y = [0 -im; im 0]
    Z = [1 0; 0 -1]
    paulis = (I, X, Y, Z)
    for a = 1:4, b = 1:4, c = 1:4, d = 1:4
        tf[a, b, c, d] = real(tr(fsim_gate * kron(paulis[a], paulis[b]) * fsim_gate' * kron(paulis[c], paulis[d])) / 4)
    end
    enable_cuda && return CuArray(tf)
    return tf
end

"""
    to_ein_code_pauli(layout, [s; haar = false, enable_cuda = false])

Generate the tensor network description for `layout` on a subset `s`. If the is a FSim gate out side `s`, 
it will be replaced by an identity gate.
"""
function to_ein_code_pauli(layout::RQCLayout{VT, PT, SCGate}, s = collect(vertices(layout)); haar = false, enable_cuda = false) where {VT, PT}
    layout = deepcopy(layout)
    update_id!(layout)
    vs = sort!(intersect(s, collect(vertices(layout))))
    nbits = length(vs)
    ids_size = Dict{Int, Int}(i => 4 for i in 1:nbits)
    open_ids_dict = Dict(vs[i] => i for i = 1:nbits)
    tensor_ids = [[id] for id in 1:nbits]
    max_id = nbits
    tensors = Any[]
    for v in vs
        if gates(layout, v)[1] === Noise
            push!(tensors, sqrt(1/2)*tensor_noise_pauli(enable_cuda = enable_cuda))
            gates(layout, v)[1] = Id
        else
            push!(tensors, tensor_init_pauli(enable_cuda = enable_cuda))
        end
    end
    last_single = Dict(vs[i] => :null for i = 1:nbits)
    D = length(gates(layout, 1))
    for d = 1:D
        for v in vs
            if gates(layout, v)[d] === Single
                ts, last_single[v] = tensor_single_pauli(last_single[v]; haar = haar, enable_cuda = enable_cuda)
                max_id += 1
                new_id_v = max_id
                push!(tensors, ts)
                push!(tensor_ids, [open_ids_dict[v], new_id_v])
                ids_size[new_id_v] = 4
                open_ids_dict[v] = new_id_v
            elseif gates(layout, v)[d] === Noise && d < D
                max_id += 1
                new_id_v = max_id
                push!(tensors, tensor_noise_pauli(enable_cuda = enable_cuda), tensor_noise_pauli(enable_cuda = enable_cuda))
                push!(tensor_ids, [open_ids_dict[v]], [new_id_v])
                ids_size[new_id_v] = 4
                open_ids_dict[v] = new_id_v
            elseif gates(layout, v)[d] === FSim
                c = depth_to_cycle(d)
                u = partner(layout, v, c)
                if u < v && u in vs
                    max_id += 2
                    new_id_u = max_id - 1
                    new_id_v = max_id
                    push!(tensors, tensor_fsim_pauli(enable_cuda = enable_cuda))
                    push!(tensor_ids, [open_ids_dict[u], open_ids_dict[v], new_id_u, new_id_v])
                    ids_size[new_id_u] = 4
                    ids_size[new_id_v] = 4
                    open_ids_dict[u] = new_id_u
                    open_ids_dict[v] = new_id_v
                end
            end
        end
    end
    for v in vs
        if gates(layout, v)[end] === Noise
            push!(tensors, sqrt(2)*tensor_noise_pauli(enable_cuda = enable_cuda))
            push!(tensor_ids, [open_ids_dict[v]])
            open_ids_dict[v] = 0
        else
            max_id += 1
            new_id_v = max_id
            push!(tensors, tensor_measure_pauli(enable_cuda = enable_cuda))
            push!(tensor_ids, [open_ids_dict[v], new_id_v])
            ids_size[new_id_v] = 2
            open_ids_dict[v] = new_id_v
        end
    end
    
    open_ids = Int[]
    for v in vs
        open_ids_dict[v] > 0 && push!(open_ids, open_ids_dict[v])
    end
    # return tensor_ids, open_ids
    return EinCode(tensor_ids, open_ids), tensors, ids_size
end
