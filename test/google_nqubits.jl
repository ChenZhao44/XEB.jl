using XEB
using OMEinsum, OMEinsumContractionOrders
using Plots

D = 16
width = 7
N1, N2 = 28, 63
xebs_2d = Float64[]

for N = N1:N2
    l = XEB.square_layout(N, D, width; start_with_zero = true)
    L = collect(1:2:N2)
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
        push!(xebs_2d, NaN)
    else
        xeb_2d = ec_opt(ts...)[] - 1
        @show xeb_2d, N
        push!(xebs_2d, xeb_2d)
    end
end

Plots.plot(N1:N2, sqrt.(xebs_2d), label = "Width = $(width), D = $(D)") |> display