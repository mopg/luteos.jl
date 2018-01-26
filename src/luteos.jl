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
`luteos` is the Latin transliteration of the ancient Greek word λυτέος, which
means "one must solve".

Max Opgenoord

Fall 2017
"""
module luteos

import IterativeSolvers
import ILU
import Dierckx # this is only used for getPressureFunc.jl, perhaps there is a better solution for that.

export CDR, Elas, Material

# problem setup
include("general/problem.jl")
include("general/setup.jl")

# integration
export Master2D, Master3D
export PG1, PG2, PG3, PG4, PG5, PG6, PG7, PG8, PG9, PG10, PG11
include("integration/pgauss.jl")
export P1, P2, P3, P4
include("integration/porder.jl")

# mesh
export Mesh2D, Mesh3D
include("mesh/mesh.jl")
include("mesh/master.jl")
include("mesh/compJacob.jl")

# links to meshers
include("meshers/mesher.jl")

# io
export writeTecplot
include("io/writeTecplotElas.jl")
include("io/writeTecplotCDR.jl")
include("io/readSU2.jl")

# solution methods
export hdgSolve
include("solve/hdgSolveElas.jl")
include("solve/hdgSolveCDR.jl")

# utilities
export pressureFunc
include("util/getPressureFunc.jl")

function __init__()

    setup()

end

end
