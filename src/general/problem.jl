# ---------------------------------------------------------------------------- #
#
#   problem.jl
#
#   Type for problem definition
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    Problem

Problem type:
Used as placeholder for information about the problem,
e.g name, source function, and boundary conditions.
"""
abstract type Problem

end

include("CDR.jl")
include("Elas.jl")
