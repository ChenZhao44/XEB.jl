using Yao, YaoExtensions

@const_gate sqrtX = ComplexF64[0.5+0.5im 0.5-0.5im; 0.5-0.5im 0.5+0.5im]
@const_gate sqrtY = ComplexF64[0.5+0.5im -0.5-0.5im; 0.5+0.5im 0.5+0.5im]
@const_gate sqrtW = ComplexF64[0.5+0.5im -sqrt(1/2)*im; sqrt(1/2) 0.5+0.5im]

function to_yao_block(layout::RQCLayout{VT, PT, SCGate}, s = collect(vertices(layout)); 
        haar = false) where {VT, PT}
    vs = sort!(intersect(collect(vertices(layout)), s))
    qubit_map = gen_qubit_map(layout, s)
    N = length(vs)
    D = depth(layout)
    blk = chain(N)
    single_prev = [I2 for i = 1:N]
    for d in 1:D
        for v in vs
            q = qubit_map[v]
            if gates(layout, v)[d] === Single
                if haar
                    M = rand(Haar(2), 2)
                    g = matblock(M)
                else
                    single_prev[q] == I2 && (g2 = rand((sqrtX, sqrtY, sqrtW)))
                    single_prev[q] == sqrtX && (g2 = rand((sqrtY, sqrtW)))
                    single_prev[q] == sqrtY && (g2 = rand((sqrtX, sqrtW)))
                    single_prev[q] == sqrtW && (g2 = rand((sqrtX, sqrtY)))
                    g = chain(Rz(rand()*2π), g2, Rz(rand()*2π))
                end
                push!(blk, put(q => g))
            elseif gates(layout, v)[d] === Noise
                g = rand((I2, X, Y, Z))
                push!(blk, put(q => g))
            elseif gates(layout, v)[d] === FSim
                u = partner(layout, v, depth_to_cycle(d))
                u in vs || continue
                q1, q2 = q, qubit_map[u]
                push!(blk, put((q1, q2) => FSimGate(π/2, π/6)))
            end
        end
    end
    return blk
end