using XEB
using OMEinsum, OMEinsumContractionOrders
using Plots

D = 16
N1, N2 = 25, 63
xebs_1d = Float64[]

for N = N1:N2
    l = XEB.brick_layout(N, D)
    L = collect(1:(N÷2))
    cuts = XEB.generate_cut(l, L, 1, D)
    XEB.simplify!(l, cuts)
    ec, ts, ids = XEB.to_ein_code_xeb(l)
    ec_opt = optimize_code(ec, ids, GreedyMethod())
    if timespace_complexity(ec_opt, ids)[2] > 28
        println("TreeSA...")
        ec_opt = optimize_code(ec_opt, ids, TreeSA())
    else
        @show timespace_complexity(ec_opt, ids)
    end
    
    if timespace_complexity(ec_opt, ids)[2] > 28
        @show timespace_complexity(ec_opt, ids)
        println("space complexity too large")
        push!(xebs_1d, NaN)
    else
        xeb_1d = ec_opt(ts...)[] - 1
        @show xeb_1d, N
        push!(xebs_1d, xeb_1d)
    end
end

Plots.plot(N1:N2, (xebs_1d), label = "D = $(D)")