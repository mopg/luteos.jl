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

type Problem

  name::String     # Name of problem

  source::Function # Rhs of problem

  bctype::Array{Int64}    # Relates each boundary to a boundary function type
  bcnorm::Bool            # Whether or not boundary functions are defined normal to the surface
  bcfunc::Array{Function} # Boundary function for each boundary face

end
