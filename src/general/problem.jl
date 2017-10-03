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
type Problem

  name::String     # Name of problem

  source::Function # Rhs of problem

  bctype::Array{Int64}    # Relates each boundary to a boundary function type
  bcnorm::Bool            # Whether or not boundary functions are defined normal to the surface
  bcfunc::Array{Function} # Boundary function for each boundary face

end

"""
    Problem(src::Function, bcf::Array{Function}, bctype::Array{Int64};
            name = "Answer to Life Universe and Everything", bcnorm = false)

Constructor for `Problem` type. `src` is a function that returns the source
function for the problem. `bcf` is an Array of Functions, which specifies the
boundary condition on each boundary surface.
"""
function Problem(src, bcf, bctype;
                 name = "Answer to Life Universe and Everything", bcnorm = false)

  setup()

  if length(bctype) != length(bcf)
    error("Problem:: bctype needs have the same size as bcf")
  end

  Problem(name, src, bctype, bcnorm, bcf)

end
