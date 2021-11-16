using XEB
using OMEinsum, OMEinsumContractionOrders
using Graphs
using Statistics

L = sort!([52, 37, 35, 32, 22, 11, 31, 21, 8, 24, 7, 1, 
    29, 18, 5, 26, 15, 6, 40, 25, 16, 44, 42, 51, 53, 48, 46])
R = setdiff(1:53, L)

g = google_layout_53(53, 20)
cuts = XEB.generate_cut(g, L, 1, 20)
XEB.simplify!(g, cuts);
ec_L, ts_L, ids_size_L = to_ein_code(g, L)
ec_R, ts_R, ids_size_R = to_ein_code(g, R)

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
    g = google_layout_53(53, 20);
    cuts = XEB.generate_cut(g, L, 1, 20)
    XEB.simplify!(g, cuts);
    _, ts_L, _ = to_ein_code(g, L; haar = true)
    _, ts_R, _ = to_ein_code(g, R; haar = true)
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

# q_remain = findall(!isequal(Noise), [gates(g, v)[end] for v = 1:nv(g)])
# intersect(q_remain, L)
# intersect(q_remain, R)

# XEB_Haar = 1.8*10^(-5)


# # ec_opt = optimize_code(ec_opt, ids_size, TreeSA())
# function small_xeb(nbits, ncycles)
#     g = google_layout_53(nbits, ncycles)
#     ec, ts, ids_size = to_ein_code(g)
#     ec_opt = optimize_code(ec, ids_size, GreedyMethod())
#     timespace_complexity(ec_opt, ids_size)
#     p = ec_opt(ts...)
#     return 2^nbits * sum(p .* p) - 1
# end

# gates(google_layout_53(1, 20), 1)

# mean(small_xeb(4, 10) for i = 1:1000)

# using RandomMatrices
# function test_haar(nbits, N)
#     xebs = 0.0
#     for i = 1:N
#         V = rand(Haar(2), 2^nbits)
#         st0 = zeros(2^nbits)
#         st0[1] = 1
#         st = V*st0
#         p = abs2.(st)
#         xeb = 2^nbits*p'*p - 1
#         xebs += xeb
#     end
#     return xebs/N
# end
# test_haar(4, 1000)