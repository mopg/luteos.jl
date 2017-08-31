# include("../src/integration/basisFuncTriangle.jl")
include("../src/mesh/master2D.jl")
include("../src/mesh/master3D.jl")


using Base.Test

println("Basis functions tests")
include("testBasisFunc.jl")

# println("Mesh tests")
# include("testMesh.jl")

println("Master tests")
include("testMaster2D.jl")
include("testMaster3D.jl")
