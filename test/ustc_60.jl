module USTC_60

using Printf
using XEB
using OMEinsum, OMEinsumContractionOrders
using CUDA
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

L, R = L4, R4
D = 24
g = ustc_layout_60(60, D)
cuts = XEB.generate_cut(g, L, 1, D)
XEB.simplify!(g, cuts);
flux_g_L = 2.0^(flux(g, L) - length(L))
flux_g_R = 2.0^(flux(g, R) - length(R))
@show flux_g_L, flux_g_R
ec_L, ts_L, ids_size_L = to_ein_code(g, L; haar = true, enable_cuda = true)
ec_R, ts_R, ids_size_R = to_ein_code(g, R; haar = true, enable_cuda = true)
length(ec_L.iy), length(ec_R.iy)

ec_L_greedy = optimize_code(ec_L, ids_size_L, GreedyMethod())
ec_L_opt = optimize_code(ec_L_greedy, ids_size_L, TreeSA())
@show timespace_complexity(ec_L_opt, ids_size_L)
ec_R_greedy = optimize_code(ec_R, ids_size_R, GreedyMethod())
ec_R_opt = optimize_code(ec_R_greedy, ids_size_R, TreeSA())
@show timespace_complexity(ec_R_opt, ids_size_R)
CUDA.@time ec_L_opt(ts_L...);
CUDA.@time ec_R_opt(ts_R...);
xebs_L = Float64[]
xebs_R = Float64[]
maxs_L = Float64[]
maxs_R = Float64[]

for i = 1:100
    g = ustc_layout_60(60, D);
    cuts = XEB.generate_cut(g, L, 1, D)
    XEB.simplify!(g, cuts);
    _, ts_L, _ = to_ein_code(g, L; haar = false, enable_cuda = true)
    _, ts_R, _ = to_ein_code(g, R; haar = false, enable_cuda = true)
    p_L = ec_L_opt(ts_L...);
    p_R = ec_R_opt(ts_R...);
    # @show sum(p_L), sum(p_R)
    # p_L = p_L / sum(p_L);
    # p_R = p_R / sum(p_R);
    xeb_L = length(p_L)*sum(p_L .* p_L) - 1
    push!(xebs_L, xeb_L)
    xeb_R = length(p_R)*sum(p_R .* p_R) - 1
    push!(xebs_R, xeb_R)
    push!(maxs_L, length(p_L)*maximum(p_L) - 1)
    push!(maxs_R, length(p_R)*maximum(p_R) - 1)
    print("$(i): ")
    @printf "XEB_max = %.3e = %.3e + %.3e (flux = %.3e + %.3e)\n" (mean(maxs_L)+mean(maxs_R)) mean(maxs_L) mean(maxs_R) flux_g_L flux_g_R
    # println("$(i): XEB_max = $(mean(maxs_L)+mean(maxs_R)) = $(mean(maxs_L)) + $(mean(maxs_R))")
    # println("$(i): XEB_L = $(mean(xebs_L)), XEB_R = $(mean(xebs_R))")
end

end