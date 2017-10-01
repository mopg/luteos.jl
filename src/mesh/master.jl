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

include("../integration/quadratureLine.jl")
include("../integration/quadratureTriangle.jl")
include("../integration/quadratureTet.jl")
include("../integration/basisFuncLineLeg.jl")
include("../integration/basisFuncTriangleLeg.jl")
include("../integration/basisFuncTetLeg.jl")
include("../integration/basisFuncLineLag.jl")
include("../integration/basisFuncTriangleLag.jl")
include("../integration/basisFuncTetLag.jl")

include("master2D.jl")
include("master3D.jl")
