tensor_init() = ComplexF64[1, 0, 0, 1] / sqrt(2)
tensor_measure() = (1/sqrt(2))*ComplexF64[1 1; 0 0; 0 0; 1 -1]
function tensor_single(g = :null)
    g === :null && (new_g = rand((:sx, :sy, :sw)))
    g === :sx && (new_g = rand((:sy, :sw)))
    g === :sy && (new_g = rand((:sx, :sw)))
    g === :sw && (new_g = rand((:sx, :sy)))

    I = [1 0; 0 1]
    X = [0 1; 1 0]
    Y = [0 -im; im 0]
    Z = [1 0; 0 -1]
    paulis = (I, X, Y, Z)
    new_g === :sx && (M = sqrt(X))
    new_g === :sy && (M = sqrt(Y))
    new_g === :sw && (M = sqrt(0.5*(X+Y)))
    M = [1 0; 0 exp(rand()*2π*im)]*M*[1 0; 0 exp(rand()*2π*im)]
    ts = zeros(ComplexF64, 4, 4)
    
    for i = 1:4, j = 1:4
        ts[i, j] = tr(M*paulis[i]*M'*paulis[j]/2)
    end
    return ts, new_g
end
function tensor_noise()
    tn = zeros(ComplexF64, 4)
    tn[1] = 1
    return tn
end
function tensor_fsim()
    fsim_gate = zeros(ComplexF64, 4, 4)
    fsim_gate[1, 1] = 1
    fsim_gate[2, 3] = -im
    fsim_gate[3, 2] = -im
    fsim_gate[4, 4] = exp(-pi/6*im)
    tf = zeros(ComplexF64, 4, 4, 4, 4)
    I = [1 0; 0 1]
    X = [0 1; 1 0]
    Y = [0 -im; im 0]
    Z = [1 0; 0 -1]
    paulis = (I, X, Y, Z)
    for a = 1:4, b = 1:4, c = 1:4, d = 1:4
        tf[a, b, c, d] = tr(fsim_gate * kron(paulis[a], paulis[b]) * fsim_gate' * kron(paulis[c], paulis[d])) / 4
    end
    return tf
end

function to_ein_code(g)
    g = deepcopy(g)
    simplify!(g)
    nbits = nv(g)
    ids_size = Dict{Int, Int}(i => 4 for i in 1:nbits)
    open_ids = [i for i = 1:nbits]
    tensor_ids = [[id] for id in open_ids]
    max_id = nbits
    tensors = Any[]
    for v in 1:nv(g)
        if gates(g, v)[1] === :noise
            push!(tensors, tensor_noise())
            gates(g, v)[1] = :id
        else
            push!(tensors, tensor_init())
        end
    end
    last_single = [:null for _ = 1:nbits]
    D = length(gates(g, 1))
    for d = 1:D
        for v in vertices(g)
            if gates(g, v)[d] === :single
                ts, last_single[v] = tensor_single(last_single[v])
                push!(tensors, ts)
                max_id += 1
                new_id_v = max_id
                push!(tensor_ids, [open_ids[v], new_id_v])
                ids_size[new_id_v] = 4
                open_ids[v] = new_id_v
            elseif gates(g, v)[d] === :noise && d < D
                push!(tensors, tensor_noise(), tensor_noise())
                max_id += 1
                new_id_v = max_id
                push!(tensor_ids, [open_ids[v]], [new_id_v])
                ids_size[new_id_v] = 4
                open_ids[v] = new_id_v
            elseif gates(g, v)[d] === :fsim
                c = d÷2
                u = partner(g, v, c)
                if u < v
                    push!(tensors, tensor_fsim())
                    max_id += 2
                    new_id_u = max_id - 1
                    new_id_v = max_id
                    push!(tensor_ids, [open_ids[u], open_ids[v], new_id_u, new_id_v])
                    ids_size[new_id_u] = 4
                    ids_size[new_id_v] = 4
                    open_ids[u] = new_id_u
                    open_ids[v] = new_id_v
                end
            end
        end
    end
    for v in nv(g):-1:1
        if gates(g, v)[end] === :noise
            push!(tensors, tensor_noise())
            push!(tensor_ids, [open_ids[v]])
            deleteat!(open_ids, v)
        else
            push!(tensors, tensor_measure())
            max_id += 1
            new_id_v = max_id
            push!(tensor_ids, [open_ids[v], new_id_v])
            ids_size[new_id_v] = 2
            open_ids[v] = new_id_v
        end
    end
    # return tensor_ids, open_ids
    return EinCode(tensor_ids, open_ids), tensors, ids_size
end