# ---------------------------------------------------------------------------- #
#
#   porder.jl
#
#   Abstract type for polynomial order
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    Porder

Porder abstract type:
Overarching abstract type for polynomial order types.
"""
abstract type Porder

end

### P = 1
type P1 <: Porder
  p::Int64
end
function P1(  )
  P1(1)
end
function PGdef( P::P1 )
  return PG3()
end

### P = 2
type P2 <: Porder
  p::Int64
end
function P2(  )
  P2(2)
end
function PGdef( P::P2 )
  return PG6()
end

### P = 3
type P3 <: Porder
  p::Int64
end
function P3(  )
  P3(3)
end
function PGdef( P::P3 )
  return PG9()
end

### P = 4
type P4 <: Porder
  p::Int64
end
function P4(  )
  P4(4)
end
function PGdef( P::P4 )
  return PG12()
end

# ### P = 5
# type P5 <: Porder
#   p::Int64
# end
#
# function P5(  )
#   P5(5)
# end
#
# function PGdef( P::P5 )
#   return PG15()
# end
