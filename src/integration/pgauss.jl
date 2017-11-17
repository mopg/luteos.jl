# ---------------------------------------------------------------------------- #
#
#   pgauss.jl
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
    PGauss

PGauss abstract type:
Overarching abstract type for polynomial order types.
"""
abstract type PGauss

end

### PG = 1
type PG1 <: PGauss
  p::Int64
end
function PG1(  )
  PG1(1)
end
function comporder( P::PG1 )
  return ( PG1(), PG1(), PG1() )
end

### PG = 2
type PG2 <: PGauss
  p::Int64
end
function PG2(  )
  PG2(2)
end
function comporder( P::PG2 )
  return ( PG3(), PG1(), PG1() )
end

### PG = 3
type PG3 <: PGauss
  p::Int64
end
function PG3(  )
  PG3(3)
end
function comporder( P::PG3 )
  return ( PG3(), PG2(), PG2() )
end

### PG = 4
type PG4 <: PGauss
  p::Int64
end
function PG4(  )
  PG4(4)
end
function comporder( P::PG4 )
  return ( PG5(), PG3(), PG3() )
end

### PG = 5
type PG5 <: PGauss
  p::Int64
end
function PG5(  )
  PG5(5)
end
function comporder( P::PG5 )
  return ( PG5(), PG4(), PG3() )
end

### PG = 6
type PG6 <: PGauss
  p::Int64
end
function PG6(  )
  PG6(6)
end
function comporder( P::PG6 )
  return ( PG7(), PG5(), PG4() )
end

### PG = 7
type PG7 <: PGauss
  p::Int64
end
function PG7(  )
  PG7(7)
end
function comporder( P::PG7 )
  return ( PG7(), PG5(), PG5() )
end

### PG = 8
type PG8 <: PGauss
  p::Int64
end
function PG8(  )
  PG8(8)
end
function comporder( P::PG8 )
  return ( PG8(), PG5(), PG6() )
end

### PG = 9
type PG9 <: PGauss
  p::Int64
end
function PG9(  )
  PG9(9)
end
function comporder( P::PG9 )
  return ( PG9(), PG8(), PG7() )
end

### PG = 10
type PG10 <: PGauss
  p::Int64
end
function PG10(  )
  PG10(10)
end
function comporder( P::PG10 )
  return ( PG11(), PG8(), PG8() )
end

### PG = 11
type PG11 <: PGauss
  p::Int64
end
function PG11(  )
  PG11(11)
end
function comporder( P::PG11 )
  return ( PG11(), PG8(), PG9() )
end

### PG = 12
type PG12 <: PGauss
  p::Int64
end
function PG12(  )
  PG12(12)
end
function comporder( P::PG12 )
  return ( PG12(), PG8(), PG10() )
end

### PG = 13
type PG13 <: PGauss
  p::Int64
end
function PG13(  )
  PG13(13)
end
function comporder( P::PG13 )
  return ( PG13(), PG9(), PG11() )
end
