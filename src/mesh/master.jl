include("../integration/quadratureLine.jl")
include("../integration/quadratureTriangle.jl")
include("../integration/quadratureTet.jl")
include("../integration/basisFuncLineLeg.jl")
include("../integration/basisFuncTriangleLeg.jl")
include("../integration/basisFuncTetLeg.jl")
include("../integration/basisFuncLineLag.jl")
include("../integration/basisFuncTriangleLag.jl")
include("../integration/basisFuncTetLag.jl")

abstract type Master

end

include("master2D.jl")
include("master3D.jl")
