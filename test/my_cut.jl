using XEB
using OMEinsumContractionOrders

g = google_layout_53(53, 20)
cuts = XEB.my_cut_53(8, 16, 53, 20)
XEB.simplify!(g, cuts)

ec, ts, ids_size = to_ein_code(g)
ec_opt = optimize_code(ec, ids_size, TreeSA())
timespace_complexity(ec_opt, ids_size)
# ec_opt(ts...)
