using XEB
using OMEinsum, OMEinsumContractionOrders
using Plots

D = 16
width = 7
D1, D2 = 0, 9
N1, N2 = 15, 32
L = collect(1:2:N2)
xebs_2d = Float64[]

nd = 4
l_9 = XEB.square_layout(27, 20, 9; start_with_zero = false)
cuts_9 = XEB.generate_cut(l_9, collect(1:2:27), 1, nd)
XEB.simplify!(l_9, cuts_9)
cuts_9 = XEB.generate_cut(l_9, collect(1:2:27), 20-nd+1, 20)
XEB.simplify!(l_9, cuts_9)
l_9 |> XEB.plot3d |> display
ec_9, ts_9, ids_9 = XEB.to_ein_code_xeb(l_9)
ec_opt_9 = optimize_greedy(ec_9, ids_9; nrepeat = 30)
timespace_complexity(ec_opt_9, ids_9)
ec_opt_9 = optimize_tree(ec_opt_9, ids_9)
timespace_complexity(ec_opt_9, ids_9)
ec_opt_9(ts_9...)[] - 1

for N = N1:N2
    l = XEB.square_layout(N, D, width; start_with_zero = false)
    L0 = intersect(XEB.vertices(l), L)
    # cuts = XEB.generate_cut(l, L0, 1, D1)
    # XEB.simplify!(l, cuts)
    cuts = XEB.generate_cut(l, L0, D2, D)
    XEB.simplify!(l, cuts)
    ec, ts, ids = XEB.to_ein_code_xeb(l)
    ec_opt = optimize_greedy(ec, ids)
    if timespace_complexity(ec_opt, ids)[2] > 28
        println("TreeSA...")
        ec_opt = optimize_tree(ec_opt, ids)
    end
    @show timespace_complexity(ec_opt, ids)
    xeb_2d = ec_opt(ts...)[] - 1
    @show xeb_2d, N
    push!(xebs_2d, xeb_2d)
end

Plots.plot!(N1:N2, (xebs_2d), label = "D = $(D) (1:$(D1), $(D2):$D)") |> display