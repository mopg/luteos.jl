# ---------------------------------------------------------------------------- #
#
#   master.jl
#
#   Abstract master element type
#   This allows for writing an n-dimensional version of solver methods
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    Master

Master abstract type:
Overarching abstract type for master element types. Currently triangle and
tetrahedron implemented.
"""
abstract type Master

end

println("quadline")
@time include("../integration/quadratureLine.jl")
println("quadtriangle")
@time include("../integration/quadratureTriangle.jl")
println("quadtet")
@time include("../integration/quadratureTet.jl")
println("basisline")
@time include("../integration/basisFuncLineLag.jl")
println("basistriangle")
@time include("../integration/basisFuncTriangleLag.jl")
println("basistet")
@time include("../integration/basisFuncTetLag.jl")

println("master2d")
@time include("master2D.jl")
println("master3d")
@time include("master3D.jl")
