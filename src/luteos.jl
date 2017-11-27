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

__precompile__()

"""
    luteos

Julia library for high-order FEM problems.
`luteos` is the Latin transliteration of the ancient greek word λυτέος, which
means "one must solve".

Max Opgenoord

Fall 2017
"""
module luteos

import IterativeSolvers

export hdgSolveElas, Material, Problem, Master2D, Master3D, Mesh2D, Mesh3D, compJacob!, writeTecplot, writeTecplotCD

include("general/material.jl")
include("general/problem.jl")
include("general/setup.jl")

export PG1, PG2, PG3, PG4, PG5, PG6, PG7, PG8, PG9, PG10, PG11
include("integration/pgauss.jl")
export P1, P2, P3, P4
include("integration/porder.jl")

include("mesh/mesh.jl")
include("mesh/master.jl")
include("mesh/compJacob.jl")

include("io/writeTecplot.jl")
include("io/writeTecplotCD.jl")

include("solve/hdgSolveElas.jl")

export hdgSolveCD
include("solve/hdgSolveCD.jl")

# TODO: Use __init__() for pretty picture

function __init__()

    setup()
    
end

end
