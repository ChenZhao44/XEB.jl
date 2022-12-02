using XEB
using OMEinsum, OMEinsumContractionOrders

D = 16
N = 18
xebs_1d_end = Float64[]
for D_cut = 1:D
    l = brick_layout(N, D)
    plot(l)
    cuts = XEB.generate_cut(l, collect(1:(N÷2)), 1, ceil(Int, D_cut/2)-1)
    XEB.simplify!(l, cuts)
    cuts = XEB.generate_cut(l, collect(1:(N÷2)), D - floor(Int, D_cut/2), D)
    XEB.simplify!(l, cuts)
    ec, ts, ids = XEB.to_ein_code_xeb(l)
    ec_opt = optimize_code(ec, ids, GreedyMethod())
    @show timespace_complexity(ec_opt, ids)
    if timespace_complexity(ec_opt, ids)[2] > 28
        println("TreeSA...")
        ec_opt = optimize_code(ec, ids, TreeSA())
    end
    @show timespace_complexity(ec_opt, ids)
    if timespace_complexity(ec_opt, ids)[2] > 28
        push!(xebs_1d_end, NaN)
        println("space complexity too large")
    else
        xeb_1d = ec_opt(ts...)[] - 1
        @show D_cut, xeb_1d
        push!(xebs_1d_end, xeb_1d)
    end
end

xebs_1d_start = Float64[]
for D_cut = 1:D
    l = brick_layout(N, D)
    plot(l)
    cuts = XEB.generate_cut(l, collect(1:(N÷2)), 1, D_cut)
    XEB.simplify!(l, cuts)
    ec, ts, ids = XEB.to_ein_code_xeb(l)
    ec_opt = optimize_code(ec, ids, GreedyMethod())
    @show timespace_complexity(ec_opt, ids)
    if timespace_complexity(ec_opt, ids)[2] > 28
        println("TreeSA...")
        ec_opt = optimize_code(ec, ids, TreeSA())
    end
    @show timespace_complexity(ec_opt, ids)
    if timespace_complexity(ec_opt, ids)[2] > 28
        push!(xebs_1d_end, NaN)
        println("space complexity too large")
    else
        xeb_1d = ec_opt(ts...)[] - 1
        @show D_cut, xeb_1d
        push!(xebs_1d_start, xeb_1d)
    end
end

Plots.plot(xebs_1d_end; yscale = :log10, label = "two ends")
Plots.plot!(xebs_1d_start; yscale = :log10, label = "one end")
Plots.plot!([xebs_1d_end[end]*100 for _ = 1:length(xebs_1d_end)]; yscale = :log10, label = "threshold")
