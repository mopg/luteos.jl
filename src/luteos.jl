# __precompile__()

module luteos

export hdgSolveElas, Material, Problem, Master2D, Master3D, Mesh2D, Mesh3D, compJacob!, writeTecplot

# include("integration/basisFuncLineLag.jl")
# include("integration/basisFuncLineLeg.jl")
# include("integration/basisFuncTriangleLag.jl")
# include("integration/basisFuncTriangleLeg.jl")
# include("integration/basisFuncTetLag.jl")
# include("integration/basisFuncTetLeg.jl")
#
# include("integration/quadratureLine.jl")
# include("integration/quadratureTriangle.jl")
# include("integration/quadratureTet.jl")

include("general/material.jl")
include("general/problem.jl")

include("mesh/mesh.jl")
include("mesh/master.jl")
include("mesh/compJacob.jl")

include("io/writeTecplot.jl")

include("solve/hdgSolveElas.jl")

end
