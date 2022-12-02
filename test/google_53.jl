# module Google_53

using Printf
using XEB
using OMEinsum, OMEinsumContractionOrders
using CUDA
using Graphs
using Statistics

L = sort!([52, 37, 35, 32, 22, 11, 31, 21, 8, 24, 7, 1, 
    29, 18, 5, 26, 15, 6, 40, 25, 16, 44, 42, 51, 53, 48, 46])
# L = setdiff(1:53, [32, 23, 34, 43, 51, 42, 44, 26])
# L = setdiff(1:53, [27, 38, 39, 17, 28, 30, 19, 23, 33, 43, 34, 45, 50, 49])
R = setdiff(1:53, L)

D = 20
layout = google_layout_53(53, D)
cuts = XEB.generate_cut(layout, L, 1, D)
XEB.simplify!(layout, cuts);

ec, ts, ids, ids_loc = XEB.to_ein_code_xeb(layout, L)
ec_opt = optimize_code(ec, ids, GreedyMethod())
xeb = ec_opt(ts...)[] - 1

ec_L, ts_L, ids_size_L = to_ein_code_pauli(layout, L; enable_cuda = false)
ec_R, ts_R, ids_size_R = to_ein_code_pauli(layout, R; enable_cuda = false)

ec_L_greedy = optimize_code(ec_L, ids_size_L, GreedyMethod())
ec_L_opt = ec_L_greedy
@show timespace_complexity(ec_L_greedy, ids_size_L)
ec_R_greedy = optimize_code(ec_R, ids_size_R, GreedyMethod())
ec_R_opt = ec_R_greedy
@show timespace_complexity(ec_R_opt, ids_size_R)

xebs_L = Float64[]
xebs_R = Float64[]
maxs_L = Float64[]
maxs_R = Float64[]

for i = 1:1000
    g = layout
    _, ts_L, _ = to_ein_code_pauli(g, L; haar = false, enable_cuda = false)
    _, ts_R, _ = to_ein_code_pauli(g, R; haar = false, enable_cuda = false)
    p_L = ec_L_opt(ts_L...);
    p_R = ec_R_opt(ts_R...);
    # @show sum(p_L), sum(p_R)
    # p_L = p_L / sum(p_L);
    # p_R = p_R / sum(p_R);
    xeb_L = 2^(length(ec_L.iy))*sum(p_L .* p_L) - 1
    push!(xebs_L, xeb_L)
    xeb_R = 2^(length(ec_R.iy))*sum(p_R .* p_R) - 1
    push!(xebs_R, xeb_R)
    push!(maxs_L, (2^(length(ec_L.iy))*maximum(p_L) - 1))
    push!(maxs_R, (2^(length(ec_R.iy))*maximum(p_R) - 1))
    println("$(i): ")
    @printf "XEB = %.3e = %.3e + %.3e \n" (mean(xebs_L)+mean(xebs_R)) mean(xebs_L) mean(xebs_R)
    @printf "XEB_max = %.3e = %.3e + %.3e \n" (mean(maxs_L)+mean(maxs_R)) mean(maxs_L) mean(maxs_R)
    println("$(i): max_L = $(mean(maxs_L)), max_R = $(mean(maxs_R)), max = $(mean(maxs_R)+mean(maxs_L))")
    println("$(i): XEB_L = $(mean(xebs_L)), XEB_R = $(mean(xebs_R))")
end

mean(xebs_L)+mean(xebs_R)

std(xebs_L) / 100
std(xebs_R) / 100

# end