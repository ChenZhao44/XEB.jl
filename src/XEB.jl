module XEB

using LinearAlgebra
using Multigraphs, Graphs
using OMEinsum
using RandomMatrices
using CUDA
using Yao, Yao.EasyBuild

export brick_layout, google_layout_53, ustc_layout_60, 
add_noise!, gates, update_noise!
include("layout.jl")
include("predefines.jl")

export simplify!
include("rewrite.jl")
include("cut.jl")

export to_ein_code_pauli, to_ein_code_xeb, to_yao_block
include("ein_code_pauli.jl")
include("ein_code_xeb.jl")
include("yao_blocks.jl")

export plot, plot3d
include("plots.jl")
include("plot3d.jl")
include("utils.jl")

end
