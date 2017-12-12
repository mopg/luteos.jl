# ---------------------------------------------------------------------------- #
#
#   CDR.jl
#
#   Struct for convection-diffusion-reaction problems
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    CDR

CDR type:
Used as placeholder for information about convection-diffusion-reaction problem,
e.g name, source function, and boundary conditions.
"""
struct CDR{S<:Function,B<:Function} <: Problem

  name::String        # Name of problem

  κ::Float64          # Diffusivity coefficient
  c::Vector{Float64}  # Convective velocity

  source::S               # Rhs of problem

  bctype::Vector{Int64}   # Relates each boundary to a boundary function type
  bcnorm::Bool            # Whether or not boundary functions are defined normal to the surface
  bcfunc::Vector{B} # Boundary function for each boundary face

end

"""
    CDR(src::Function, bcf::Array{Function}, bctype::Array{Int64};
        name = "Answer to Life Universe and Everything", bcnorm = false)

Constructor for `Problem` type. `src` is a function that returns the source
function for the problem. `bcf` is an Array of Functions, which specifies the
boundary condition on each boundary surface.
"""
function CDR(κ::Float64,
             c::Vector{Float64},
             src::Function,
             bcf::Vector{Function},
             bctype::Vector{Int64};
             name = "Answer to Life Universe and Everything", bcnorm = false)

  if length(bctype) != length(bcf)
    error("Problem:: bctype needs have the same size as bcf")
  end

  CDR(name, κ, c, src, bctype, bcnorm, bcf)

end
