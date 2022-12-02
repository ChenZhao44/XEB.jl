using XEB
using OMEinsum, OMEinsumContractionOrders

D = 20
width = 9
N = 53
R1 = collect(2:2:N)
R2 = collect(10:2:N)
R3 = collect(16:2:N)
R4 = collect(26:2:N)
R5 = setdiff(2:2:N, [8])
R6 = setdiff(2:2:N, [4, 8])
R7 = sort([53, 35, 17, 18, 36, 
    54, 37, 19, 1, 16, 34, 52, 
    51, 33, 15, 2, 20, 38,
    39, 21, 3, 14, 32, 50, 
    49, 31, 13, 4])
xebs_2d = Dict()

for R in [R1, R2, R3, R4, R5, R6, R7]
    l = XEB.square_layout(N, D, width; start_with_zero = false)
    L = setdiff(1:N, R)
    cuts = XEB.generate_cut(l, L, 1, D)
    XEB.simplify!(l, cuts)
    ec, ts, ids = XEB.to_ein_code_xeb(l)
    ec_opt = optimize_code(ec, ids, GreedyMethod())
    if timespace_complexity(ec_opt, ids)[2] > 28
        println("TreeSA...")
        ec_opt = optimize_code(ec, ids, TreeSA())
    end
    @show timespace_complexity(ec_opt, ids)
    if timespace_complexity(ec_opt, ids)[2] > 28
        xebs_2d[R] = nothing
        println("space complexity too large")
    else
        xeb_2d = ec_opt(ts...)[] - 1
        @show xeb_2d
        xebs_2d[R] = xeb_2d
    end
end

for R in [R1, R2, R3, R4, R5, R6, R7]
    @show xebs_2d[R]
end

# XEBs[R1]: 1.5226433891513125e-8
# XEBs[R2]: 1.5078167159288114e-8
# XEBs[R3]: 1.5537246156327456e-8
# XEBs[R4]: 1.2291113238305229e-8
# XEBs[R5]: 1.5260782859627398e-8
# XEBs[R6]: 1.5249234763814457e-8
# XEBs[R7]: 1.6837453431506333e-8