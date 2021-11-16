using RandomMatrices

tensor_init() = Float64[1, 0, 0, 1] / sqrt(2)
tensor_measure() = (1/sqrt(2))*Float64[1 1; 0 0; 0 0; 1 -1]
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
    new_g === :sw && (M = sqrt(sqrt(1/2)*(X+Y)))
    M = [1 0; 0 exp(rand()*2π*im)]*M*[1 0; 0 exp(rand()*2π*im)]
    M = rand(Haar(2), 2)    # Haar random matrix
    ts = zeros(Float64, 4, 4)
    
    for i = 1:4, j = 1:4
        ts[i, j] = real(tr(M*paulis[i]*M'*paulis[j]/2))
    end
    return ts, new_g
end
function tensor_noise()
    tn = zeros(Float64, 4)
    tn[1] = 1
    return tn
end
function tensor_fsim()
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
    return tf
end

function to_ein_code(g, s = collect(vertices(g)))
    g = deepcopy(g)
    simplify!(g)
    vs = sort!(intersect(s, collect(vertices(g))))
    nbits = length(vs)
    ids_size = Dict{Int, Int}(i => 4 for i in 1:nbits)
    open_ids_dict = Dict(vs[i] => i for i = 1:nbits)
    tensor_ids = [[id] for id in 1:nbits]
    max_id = nbits
    tensors = Any[]
    for v in vs
        if gates(g, v)[1] === :noise
            push!(tensors, sqrt(1/2)*tensor_noise())
            gates(g, v)[1] = :id
        else
            push!(tensors, tensor_init())
        end
    end
    last_single = Dict(vs[i] => :null for i = 1:nbits)
    D = length(gates(g, 1))
    for d = 1:D
        for v in vs
            if gates(g, v)[d] === :single
                ts, last_single[v] = tensor_single(last_single[v])
                max_id += 1
                new_id_v = max_id
                push!(tensors, ts)
                push!(tensor_ids, [open_ids_dict[v], new_id_v])
                ids_size[new_id_v] = 4
                open_ids_dict[v] = new_id_v
            elseif gates(g, v)[d] === :noise && d < D
                max_id += 1
                new_id_v = max_id
                push!(tensors, tensor_noise(), tensor_noise())
                push!(tensor_ids, [open_ids_dict[v]], [new_id_v])
                ids_size[new_id_v] = 4
                open_ids_dict[v] = new_id_v
            elseif gates(g, v)[d] === :fsim
                c = d÷2
                u = partner(g, v, c)
                if u < v && u in vs
                    max_id += 2
                    new_id_u = max_id - 1
                    new_id_v = max_id
                    push!(tensors, tensor_fsim())
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
        if gates(g, v)[end] === :noise
            push!(tensors, sqrt(2)*tensor_noise())
            push!(tensor_ids, [open_ids_dict[v]])
            open_ids_dict[v] = 0
        else
            max_id += 1
            new_id_v = max_id
            push!(tensors, tensor_measure())
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

function conn_comps(g)
    vs_remain = Set(collect(vertices(g)))
    comps = []
    while !isempty(vs_remain)
        v0 = pop!(vs_remain)
        comp_v0 = [v0]
        comp_remain = [v0]
        while !isempty(comp_remain)
            v = pop!(comp_remain)
            gates_v = gates(g, v)
            for i = 1:length(gates_v)
                gates_v[i] === :fsim || continue
                u = partner(g, v, i÷2)
                !(u in comp_v0) && push!(comp_remain, u)
                !(u in comp_v0) && push!(comp_v0, u)
                delete!(vs_remain, u)
            end
        end
        push!(comps, comp_v0)
    end
    return comps
end