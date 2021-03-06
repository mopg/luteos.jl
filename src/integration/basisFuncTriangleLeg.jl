# ---------------------------------------------------------------------------- #
#
#   basisFuncTriangleLeg.jl
#
#   Get Legendre basis functions for Triangle
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

function basisFuncTriangleLeg( P::P1, s::Vector{Float64}, t::Vector{Float64} )

  dim = 2
  nx  = length( s )

  # Linear

  nphi = 3
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  # value
  phi[:,1] = 1
  phi[:,2] = sqrt(6).*(s - t)
  phi[:,3] = 3.*sqrt(2).*(s + t - 2./3.)

  # derivative in x
  dphi[:,1,1] =  0
  dphi[:,2,1] =  sqrt(6)
  dphi[:,3,1] =  3.*sqrt(2)

  # derivative in y
  dphi[:,1,1] =  0
  dphi[:,2,1] = -sqrt(6)
  dphi[:,3,1] =  3.*sqrt(2)

  # transpose
  phi   = phi'

  dphi2 = fill( 0.0, size(dphi,2), size(dphi,1), 2 )
  dphi2[:,:,1] = dphi[:,:,1]'
  dphi2[:,:,2] = dphi[:,:,2]'

  return phi, dphi2

end

function basisFuncTriangleLeg( P::P2, s::Vector{Float64}, t::Vector{Float64} )

  dim = 2
  nx  = length( s )

  # Quadratric

  nphi = 6
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  # value
  phi[:,1] = 1
  phi[:,2] = sqrt(6).*(s - t)
  phi[:,3] = 3.*sqrt(2).*(s + t - 2./3.)
  phi[:,4] = 3.*(s - t).*(4 - 5.*(s + t))
  phi[:,5] = sqrt(3).*(3 - 12.*(s + t) + 10.*(s + t).*(s + t))
  phi[:,6] = sqrt(15).*((s + t).*(s + t) - 3.*(s - t).*(s - t))/2.

  # derivative in x
  dphi[:,1,1] =  0
  dphi[:,2,1] =  sqrt(6)
  dphi[:,3,1] =  3.*sqrt(2)
  dphi[:,4,1] =  6.*(2 - 5.*s)
  dphi[:,5,1] = -4.*sqrt(3).*(3 - 5.*(s + t))
  dphi[:,6,1] = -2.*sqrt(15).*(s - 2.*t)

  # derivative in y
  dphi[:,1,1] =  0
  dphi[:,2,1] = -sqrt(6)
  dphi[:,3,1] =  3.*sqrt(2)
  dphi[:,4,1] = -6.*(2 - 5.*t)
  dphi[:,5,1] = -4.*sqrt(3).*(3 - 5.*(s + t))
  dphi[:,6,1] = -2.*sqrt(15).*(t - 2.*s)

  # transpose
  phi   = phi'

  dphi2 = fill( 0.0, size(dphi,2), size(dphi,1), 2 )
  dphi2[:,:,1] = dphi[:,:,1]'
  dphi2[:,:,2] = dphi[:,:,2]'

  return phi, dphi2

end

function basisFuncTriangleLeg( P::P3, s::Vector{Float64}, t::Vector{Float64} )

  dim = 2
  nx  = length( s )

  nphi = 10
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  # value
  phi[:, 1] = 1
  phi[:, 2] = sqrt(6).*(s - t)
  phi[:, 3] = 3.*sqrt(2).*(s + t - 2./3.)
  phi[:, 4] = 3.*(s - t).*(4 - 5.*(s + t))
  phi[:, 5] = sqrt(3).*(3 - 12.*(s + t) + 10.*(s + t).*(s + t))
  phi[:, 6] = sqrt(15).*((s + t).*(s + t) - 3.*(s - t).*(s - t))/2.
  phi[:, 7] = 2.*sqrt(3).*(s - t).*(10 - 30.*(s + t) + 21.*(s + t).*(s + t))
  phi[:, 8] =   sqrt(7).*(s - t).*(5.*(s - t).*(s - t) - 3.*(s + t).*(s + t))
  phi[:, 9] = 2.*(5.*((7.*(s + t) - 12).*(s + t) + 6).*(s + t) - 4)
  phi[:,10] = sqrt(5).*((s + t).*(s + t) - 3.*(s - t).*(s - t)).*(6 - 7.*(s + t))

  # derivative in x
  dphi[:, 1,1] =  0
  dphi[:, 2,1] =  sqrt(6)
  dphi[:, 3,1] =  3.*sqrt(2)
  dphi[:, 4,1] =  6.*(2 - 5.*s)
  dphi[:, 5,1] = -4.*sqrt(3).*(3 - 5.*(s + t))
  dphi[:, 6,1] = -2.*sqrt(15).*(s - 2.*t)
  dphi[:, 7,1] =  2.*sqrt(3).*(10.*(1 - 6.*s) + 21.*(s + t).*(3.*s - t))
  dphi[:, 8,1] =  6.*sqrt(7).*(3.*(s - t).*(s - t) - 2.*s.*s)
  dphi[:, 9,1] =  30.*(2 - 8.*(s + t) + 7.*(s + t).*(s + t))
  dphi[:,10,1] =  6.*sqrt(5).*(7.*((s - t).*(s - t) - 2.*t.*t) - 4.*(s - 2.*t))

  dphi[:, 1,1] =  0
  dphi[:, 2,1] = -sqrt(6)
  dphi[:, 3,1] =  3.*sqrt(2)
  dphi[:, 4,1] = -6.*(2 - 5.*t)
  dphi[:, 5,1] = -4.*sqrt(3).*(3 - 5.*(s + t))
  dphi[:, 6,1] = -2.*sqrt(15).*(t - 2.*s)
  dphi[:, 7,1] = -2.*sqrt(3).*(10.*(1 - 6.*t) - 21.*(s + t).*(s - 3.*t))
  dphi[:, 8,1] = -6.*sqrt(7).*(3.*(s - t).*(s - t) - 2.*t.*t)
  dphi[:, 9,1] =  30.*(2 - 8.*(s + t) + 7.*(s + t).*(s + t))
  dphi[:,10,2] =  6.*sqrt(5).*(7.*((s - t).*(s - t) - 2.*s.*s) - 4.*(t - 2.*s))

  # transpose
  phi   = phi'

  dphi2 = fill( 0.0, size(dphi,2), size(dphi,1), 2 )
  dphi2[:,:,1] = dphi[:,:,1]'
  dphi2[:,:,2] = dphi[:,:,2]'

  return phi, dphi2

end

function basisFuncTriangleLeg( P::P4, s::Vector{Float64}, t::Vector{Float64} )

  dim = 2
  nx  = length( s )

  nphi = 15
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  u = s + t
  v = s - t

  # value
  phi[:, 1] = 1

  phi[:, 2] = sqrt(6).*v
  phi[:, 3] = 3.*sqrt(2).*(u - 2./3.)

  phi[:, 4] = 3.*v.*(4 - 5.*u)
  phi[:, 5] = sqrt(3).*(3 - 12.*u + 10.*u.*u)
  phi[:, 6] = sqrt(15).*(u.*u - 3.*v.*v)/2.

  phi[:, 7] = 2.*sqrt(3).*v.*(10 - 30.*u + 21.*u.*u)
  phi[:, 8] =   sqrt(7).*v.*(5.*v.*v - 3.*u.*u)
  phi[:, 9] = 2.*(5.*((7.*u - 12).*u + 6).*u - 4)
  phi[:,10] = sqrt(5).*(u.*u - 3.*v.*v).*(6 - 7.*u)

  phi[:,11] = sqrt(15).*v.*(20 - 105.*u + 168.*u.*u - 84.*u.*u.*u)
  phi[:,12] = 0.5.*sqrt(35).*v.*(5.*v.*v - 3.*u.*u).*(8 - 9.*u)
  phi[:,13] = sqrt(5).*(5 - 60.*u + 210.*u.*u - 280.*u.*u.*u + 126.*u.*u.*u.*u)
  phi[:,14] = 2.5.*(u.*u - 3.*v.*v).*(21 - 4.*u.*(14 - 9.*u))
  phi[:,15] = 0.375.*sqrt(5).*(3.*u.*u.*u.*u - 30.*u.*u.*v.*v + 35.*v.*v.*v.*v)

  # derivative in x
  dphi[:, 1,1] =  0

  dphi[:, 2,1] =  sqrt(6)
  dphi[:, 3,1] =  3.*sqrt(2)

  dphi[:, 4,1] =  3.*(4 - 5.*(u + v))
  dphi[:, 5,1] = -4.*sqrt(3).*(3 - 5.*u)
  dphi[:, 6,1] =  sqrt(15).*(u - 3.*v)

  dphi[:, 7,1] =  2.*sqrt(3).*(10 - 30.*v - 3.*u.*(10 - 7.*(u + 2.*v)))
  dphi[:, 8,1] = -3.*sqrt(7).*(u.*u + 2.*u.*v - 5.*v.*v)
  dphi[:, 9,1] =  30.*(2 - 8.*u + 7.*u.*u)
  dphi[:,10,1] = -3.*sqrt(5).*(7.*u.*u + v.*(12 - 7.*v) - 2.*u.*(2 + 7.*v))

  dphi[:,11,1] =  sqrt(15).*(5.*(4 - 21.*v) - 21.*u.*(5 - 16.*v - 4.*u.*(2 - u - 3.*v)))
  dphi[:,12,1] = -1.5.*sqrt(35).*(u.*u.*(8 - 9.*u - 27.*v) + u.*v.*(16 + 45.*v) - 5.*v.*v.*(8 - 3.*v))
  dphi[:,13,1] = -12.*sqrt(5).*(5 - 7.*u.*(5 - 2.*u.*(5 - 3.*u)))
  dphi[:,14,1] = -15.*(4.*u.*u.*(7 - 6.*u + 9.*v) - u.*(7 + 4.*v.*(14 - 9.*v)) + 7.*v.*(3 - 4.*v))
  dphi[:,15,1] =  1.5.*sqrt(5).*(3.*u.*u.*(u - 5.*v) - 5.*v.*v.*(3.*u - 7.*v))

  # derivative in y
  dphi[:, 1,2] =  0

  dphi[:, 2,2] = -sqrt(6)
  dphi[:, 3,2] =  3.*sqrt(2)

  dphi[:, 4,2] = -3.*(4 - 5.*(u - v))
  dphi[:, 5,2] = -4.*sqrt(3).*(3 - 5.*u)
  dphi[:, 6,2] =  sqrt(15).*(u + 3.*v)

  dphi[:, 7,2] = -2.*sqrt(3).*(10 + 30.*v - 3.*u.*(10 - 7.*(u - 2.*v)))
  dphi[:, 8,2] =  3.*sqrt(7).*(u.*u - 2.*u.*v - 5.*v.*v)
  dphi[:, 9,2] =  30.*(2 - 8.*u + 7.*u.*u)
  dphi[:,10,2] = -3.*sqrt(5).*(7.*u.*u - v.*(12 + 7.*v) - 2.*u.*(2 - 7.*v))

  dphi[:,11,2] = -sqrt(15).*(5.*(4 + 21.*v) - 21.*u.*(5 + 16.*v - 4.*u.*(2 - u + 3.*v)))
  dphi[:,12,2] =  1.5.*sqrt(35).*(u.*u.*(8 - 9.*u + 27.*v) - u.*v.*(16 - 45.*v) - 5.*v.*v.*(8 + 3.*v))
  dphi[:,13,2] = -12.*sqrt(5).*(5 - 7.*u.*(5 - 2.*u.*(5 - 3.*u)))
  dphi[:,14,2] = -15.*(4.*u.*u.*(7 - 6.*u - 9.*v) - u.*(7 - 4.*v.*(14 + 9.*v)) - 7.*v.*(3 + 4.*v))
  dphi[:,15,2] =  1.5.*sqrt(5).*(3.*u.*u.*(u + 5.*v) - 5.*v.*v.*(3.*u + 7.*v))

  # transpose
  phi   = phi'

  dphi2 = fill( 0.0, size(dphi,2), size(dphi,1), 2 )
  dphi2[:,:,1] = dphi[:,:,1]'
  dphi2[:,:,2] = dphi[:,:,2]'

  return phi, dphi2

end
