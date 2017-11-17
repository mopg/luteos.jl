# ---------------------------------------------------------------------------- #
#
#   quadratureTriangle.jl
#
#   Get Legendre quadrature points for Triangle
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

function quadratureTriangle( P::PG1 )

  nq = 1
  weight = Array{Float64}(nq)
  xc     = Array{Float64}(nq,1)
  yc     = Array{Float64}(nq,1)

  weight[1] = 1
  xc[1] = 1./3.
  yc[1] = 1./3.

  pts = hcat( xc, yc )

  weight = weight/2

  return pts, weight

end

function quadratureTriangle( P::PG2 )

  nq = 3

  weight = Array{Float64}(nq)
  xc     = Array{Float64}(nq,1)
  yc     = Array{Float64}(nq,1)

  weight[1] = 1./3.
  xc[1]     = 1./6.
  yc[1]     = 2./3.

  weight[2] = 1./3.
  xc[2]     = 2./3.
  yc[2]     = 1./6.

  weight[3] = 1./3.
  xc[3]     = 1./6.
  yc[3]     = 1./6.

  pts = hcat( xc, yc )

  weight = weight/2

  return pts, weight

end

function quadratureTriangle( P::PG3 )

  nq = 6
  weight = Array{Float64}(nq)
  xc     = Array{Float64}(nq,1)
  yc     = Array{Float64}(nq,1)

  w = 2./12.
  a = 0.65902762237409221517838077125540
  b = 0.23193336855303057249678456117469
  c = 0.10903900907287721232483466756991

  weight[1] = w
  xc[1]     = a
  yc[1]     = b

  weight[2] = w
  xc[2]     = b
  yc[2]     = c

  weight[3] = w
  xc[3]     = c
  yc[3]     = a

  weight[4] = w
  xc[4]     = a
  yc[4]     = c

  weight[5] = w
  xc[5]     = b
  yc[5]     = a

  weight[6] = w
  xc[6]     = c
  yc[6]     = b

  pts = hcat( xc, yc )

  weight = weight/2

  return pts, weight

end

function quadratureTriangle( P::PG4 )

  nq = 6
  weight = Array{Float64}(nq)
  xc     = Array{Float64}(nq,1)
  yc     = Array{Float64}(nq,1)


  x = 4./9. - sqrt(10)/18. + sqrt(950 - 220*sqrt(10))/90.
  w = 2.*(1./12. + 3*sqrt(95 - 22*sqrt(10))/496. - sqrt(950 - 220*sqrt(10))/7440.)
  xc[1]     = x
  yc[1]     = x
  weight[1] = w

  xc[2]     = x
  yc[2]     = 1 - 2*x
  weight[2] = w

  xc[3]     = 1 - 2*x
  yc[3]     = x
  weight[3] = w

  x         = 4./9. - sqrt(10)/18. - sqrt(950 - 220*sqrt(10))/90.
  w         = 2.*(1./12. - 3*sqrt(95 - 22*sqrt(10))/496. + sqrt(950 - 220*sqrt(10))/7440.)
  xc[4]     = x
  yc[4]     = x
  weight[4] = w

  xc[5]     = x
  yc[5]     = 1 - 2*x
  weight[5] = w

  xc[6]     = 1 - 2*x
  yc[6]     = x
  weight[6] = w

  pts = hcat( xc, yc )

  weight = weight/2

  return pts, weight

end

function quadratureTriangle( P::PG5 )

  nq = 7
  weight = Array{Float64}(nq)
  xc     = Array{Float64}(nq,1)
  yc     = Array{Float64}(nq,1)

  x         = 1./3.
  w         = 2.*9./80.
  xc[1]     = x
  yc[1]     = x
  weight[1] = w

  x         = 2./7. + sqrt(15)/21.
  w         = 2.*(31./480. + sqrt(15)/2400.)
  xc[2]     = x
  yc[2]     = x
  weight[2] = w

  xc[3]     = x
  yc[3]     = 1 - 2*x
  weight[3] = w

  xc[4]     = 1 - 2*x
  yc[4]     = x
  weight[4] = w

  x         = 2./7. - sqrt(15)/21.
  w         = 2.*(31./480. - sqrt(15)/2400.)
  xc[5]     = x
  yc[5]     = x
  weight[5] = w

  xc[6]     = x
  yc[6]     = 1 - 2*x
  weight[6] = w

  xc[7]     = 1 - 2*x
  yc[7]     = x
  weight[7] = w

  pts = hcat( xc, yc )

  weight = weight/2

  return pts, weight

end

function quadratureTriangle( P::PG8 )

  nq = 16
  weight = Array{Float64}(nq)
  xc     = Array{Float64}(nq,1)
  yc     = Array{Float64}(nq,1)

  x         = 1./3.
  w         = 2.*0.072157803838893584125545555244532
  xc[1]     = x
  yc[1]     = x
  weight[1] = w

  x         = 0.17056930775176020662229350149146
  w         = 2.*0.051608685267359125140895775146065
  xc[2]     = x
  yc[2]     = x
  weight[2] = w

  xc[3]     = x
  yc[3]     = 1 - 2*x
  weight[3] = w

  xc[4]     = 1 - 2*x
  yc[4]     = x
  weight[4] = w

  x         = 0.45929258829272315602881551449417
  w         = 2.*0.047545817133642312396948052194292
  xc[5]     = x
  yc[5]     = x
  weight[5] = w

  xc[6]     = x
  yc[6]     = 1 - 2*x
  weight[6] = w

  xc[7]     = 1 - 2*x
  yc[7]     = x
  weight[7] = w

  x         = 0.050547228317030975458423550596599
  w         = 2.*0.016229248811599040155462964170890
  xc[8]     = x
  yc[8]     = x
  weight[8] = w

  xc[9]     = x
  yc[9]     = 1 - 2*x
  weight[9] = w

  xc[10]     = 1 - 2*x
  yc[10]     = x
  weight[10] = w

  x          = 0.72849239295540428124100037917606
  y          = 0.26311282963463811342178578628464
  w          = 2.*0.0136151570872174971324223450369545
  xc[11]     = x
  yc[11]     = y
  weight[11] = w

  xc[12]     = y
  yc[12]     = x
  weight[12] = w

  xc[13]     = x
  yc[13]     = 1 - x - y
  weight[13] = w

  xc[14]     = 1 - x - y
  yc[14]     = x
  weight[14] = w

  xc[15]     = y
  yc[15]     = 1 - x - y
  weight[15] = w

  xc[16]     = 1 - x - y
  yc[16]     = y
  weight[16] = w

  pts = hcat( xc, yc )

  weight = weight/2

  return pts, weight

end
