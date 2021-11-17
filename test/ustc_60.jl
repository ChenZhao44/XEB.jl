using XEB
using OMEinsum, OMEinsumContractionOrders
using Graphs
using Statistics

vs_enable = setdiff(collect(0:65), [4, 5, 11, 17, 22, 65])
L1 = [0, 1, 2, 3, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 18, 19, 20,
    21, 24, 25, 26, 30, 31, 32, 36, 37, 42, 43, 48, 54]
R1 = setdiff(vs_enable, L1)
L2 = [0, 1, 2, 3, 6, 7, 8, 9, 10, 12, 13, 14, 15, 18, 19, 20, 21, 24, 25, 
    26, 27, 30, 31, 32, 33, 36, 37, 38, 42, 43, 44, 48, 49, 54, 55, 60]
R2 = setdiff(vs_enable, L2)
L3 = [0, 1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21, 24, 25, 
    26, 27, 30, 31, 32, 33, 36, 37, 38, 42, 43, 44, 48, 49, 54, 55, 60]
R3 = setdiff(vs_enable, L3)
L4 = [0, 1, 6, 7, 8, 12, 13, 14, 18, 19, 20, 21, 24, 25, 
    26, 30, 31, 32, 36, 37, 42, 43, 48, 49, 54, 55, 60]
R4 = setdiff(vs_enable, L4)

L_ustc, R_ustc = L4, R4
D_ustc = 24
g_ustc = ustc_layout_60(60, D_ustc)
cuts_ustc = XEB.generate_cut(g_ustc, L_ustc, 1, D_ustc)
XEB.simplify!(g_ustc, cuts_ustc);
connected_components(g_ustc)
ec_L, ts_L, ids_size_L = to_ein_code(g_ustc, L_ustc)
ec_R, ts_R, ids_size_R = to_ein_code(g_ustc, R_ustc)

ec_L_greedy = optimize_code(ec_L, ids_size_L, GreedyMethod())
ec_L_opt = optimize_code(ec_L_greedy, ids_size_L, TreeSA())
timespace_complexity(ec_L_opt, ids_size_L)
ec_R_greedy = optimize_code(ec_R, ids_size_R, GreedyMethod())
ec_R_opt = optimize_code(ec_R_greedy, ids_size_R, TreeSA())
timespace_complexity(ec_R_opt, ids_size_R)
@time ec_L_opt(ts_L...);
@time ec_R_opt(ts_R...);
xebs_L = Float64[]
xebs_R = Float64[]
xebs = Float64[]

for i = 1:100
    g = google_layout_53(53, D_ustc);
    cuts = XEB.generate_cut(g, L_ustc, 1, D_ustc)
    XEB.simplify!(g, cuts);
    _, ts_L, _ = to_ein_code(g, L_ustc; haar = true)
    _, ts_R, _ = to_ein_code(g, R_ustc; haar = true)
    p_L = ec_L_opt(ts_L...);
    p_R = ec_R_opt(ts_R...);
    # @show sum(p_L), sum(p_R)
    # p_L = p_L / sum(p_L);
    # p_R = p_R / sum(p_R);
    xeb_L = 2^(length(ec_L.iy))*sum(p_L .* p_L) - 1
    push!(xebs_L, xeb_L)
    xeb_R = 2^(length(ec_R.iy))*sum(p_R .* p_R) - 1
    push!(xebs_R, xeb_R)
    xeb = xeb_L + xeb_R
    push!(xebs, xeb)
    println("$(i): XEB_L = $(mean(xebs_L)), XEB_R = $(mean(xebs_R))")
end
