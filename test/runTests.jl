# include("../src/integration/basisFuncTriangle.jl")
include("../src/mesh/master2D.jl")
include("../src/mesh/mesh2D.jl")
include("../src/mesh/master3D.jl")


using Base.Test

println("BASIS FUNCTIONS TEST")
include("testBasisFunc.jl")

# println("MESH TESTS")
# include("testMesh.jl")

println("2D MASTER TESTS")
include("testMaster2D.jl")

println("3D MASTER TESTS")
include("testMaster3D.jl")
