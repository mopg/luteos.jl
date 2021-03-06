# ---------------------------------------------------------------------------- #
#
#   Elas.jl
#
#   Struct for elasticity problems
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

include("material.jl")

"""
    Elas

Elas type:
Used as placeholder for information about structural elasticity problems,
e.g name, source function, and boundary conditions.
"""
struct Elas{S<:Function,B<:Function} <: Problem

  name::String     # Name of problem

  mat::Material    # Material used

  source::S               # Rhs of problem

  bctype::Vector{Int64}   # Relates each boundary to a boundary function type
  bcnorm::Bool            # Whether or not boundary functions are defined normal to the surface
  bcfunc::Vector{B} # Boundary function for each boundary face

end

"""
    Elas(src::Function, bcf::Array{Function}, bctype::Array{Int64},
         mat::Material; name = "Answer to Life Universe and Everything",
         bcnorm = false)

Constructor for `Problem` type. `src` is a function that returns the source
function for the problem. `bcf` is an Array of Functions, which specifies the
boundary condition on each boundary surface.
"""
function Elas(src::Function, bcf::Vector{Function}, bctype::Vector{Int64},
              mat::Material; name = "Answer to Life Universe and Everything",
              bcnorm = false)

  if length(bctype) != length(bcf)
    error("Problem:: bctype needs have the same size as bcf")
  end

  Elas(name, src, bctype, bcnorm, bcf)

end
