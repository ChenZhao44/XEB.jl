using XEB
using OMEinsum, OMEinsumContractionOrders

D = 12
width = 5
N = 30
L = collect(1:2:N)
R = setdiff(1:N, L)

l = XEB.square_layout(N, D, width; start_with_zero = true)
cuts = XEB.generate_cut(l, L, 1, D)
XEB.simplify!(l, cuts)
ec, ts, ids = XEB.to_ein_code_xeb(l)
ec_opt = optimize_code(ec, ids, GreedyMethod())
@show timespace_complexity(ec_opt, ids)
xeb = ec_opt(ts...)[] - 1

xebs_pauli = Float64[]
nsample = 1e3

ec_L, _, ids_L = to_ein_code_pauli(l, L)
ec_opt_L = optimize_code(ec_L, ids_L, GreedyMethod())
ec_R, _, ids_R = to_ein_code_pauli(l, R)
ec_opt_R = optimize_code(ec_R, ids_R, GreedyMethod())

for i = 1:nsample
    @show i
    _, ts_L, _ = to_ein_code_pauli(l, L)
    p_L = ec_opt_L(ts_L...)[:]
    xeb_L = 2^(length(ec_L.iy))*sum(p_L .* p_L) - 1

    _, ts_R, _ = to_ein_code_pauli(l, R)
    p_R = ec_opt_R(ts_R...)[:]
    xeb_R = 2^(length(ec_R.iy))*sum(p_R .* p_R) - 1
    push!(xebs_pauli, xeb_L + xeb_R)
end
