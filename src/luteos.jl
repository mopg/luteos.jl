# ---------------------------------------------------------------------------- #
#
#   luteos.jl
#
#   High-order FEM Library in Julia
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

# __precompile__()

"""
    luteos

Julia library for high-order FEM problems.
`luteos` is the Latin transliteration of the ancient greek word λυτέος, which
means "one must solve".

Max Opgenoord

Fall 2017
"""
module luteos

export hdgSolveElas, Material, Problem, Master2D, Master3D, Mesh2D, Mesh3D, compJacob!, writeTecplot

include("general/material.jl")
include("general/problem.jl")
include("general/setup.jl")

include("mesh/mesh.jl")
include("mesh/master.jl")
include("mesh/compJacob.jl")

include("io/writeTecplot.jl")

include("solve/hdgSolveElas.jl")

end
