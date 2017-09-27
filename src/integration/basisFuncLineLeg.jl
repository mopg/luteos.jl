# ---------------------------------------------------------------------------- #
#
#   basisFuncLineLag.jl
#
#   Get Legendre basis functions for Line
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

function basisFuncLineLeg( ::Type{Val{0}}, s::Array{Float64} )

  dim = 1
  nx  = length( s )

  # Constant
  nphi = 1
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  # value
  phi[:,1] = 1.0

  # derivative in x
  dphi[:,1,1] = 0.0

  # transpose
  phi   = phi'
  dphi  = dphi'

  return phi, dphi

end

function basisFuncLineLeg( ::Type{Val{1}}, s::Array{Float64} )

  dim = 1
  nx  = length( s )

  # Linear
  nphi = 2
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  t = 1 - s

  phi[:,1] = 1
  phi[:,2] = sqrt(3).*(s - t)

  dphi[:,1,1] =  0
  dphi[:,2,1] =  2.*sqrt(3)

  # transpose
  phi   = phi'
  dphi  = dphi'

  return phi, dphi

end


function basisFuncLineLeg( ::Type{Val{2}}, s::Array{Float64} )

  dim = 1
  nx  = length( s )

  # Quadratric
  nphi = 3
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  t = 1 - s

  phi[:,1] = 1
  phi[:,2] = sqrt(3).*(s - t)
  phi[:,3] = sqrt(5).*(1 - 6.*s.*t)

  dphi[:,1,1] =  0
  dphi[:,2,1] =  2.*sqrt(3)
  dphi[:,3,1] =  6.*sqrt(5).*(s - t)

  # transpose
  phi   = phi'
  dphi  = dphi'

  return phi, dphi

end

function basisFuncLineLeg( ::Type{Val{3}}, s::Array{Float64} )

  dim = 1
  nx  = length( s )

  # Cubic
  nphi = 4
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  t = 1 - s

  phi[:,1] = 1
  phi[:,2] = sqrt(3).*(s - t)
  phi[:,3] = sqrt(5).*(1 - 6.*s.*t)
  phi[:,4] = sqrt(7).*(1 - 10.*s.*t).*(s - t)

  dphi[:,1,1] =  0
  dphi[:,2,1] =  2.*sqrt(3)
  dphi[:,3,1] =  6.*sqrt(5).*(s - t)
  dphi[:,4,1] = 12.*sqrt(7).*(1 - 5.*s.*t)

  # transpose
  phi   = phi'
  dphi  = dphi'

  return phi, dphi

end

function basisFuncLineLeg( ::Type{Val{4}}, s::Array{Float64} )

  dim = 1
  nx  = length( s )

  # Quartic
  nphi = 5
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  t = 1 - s

  phi[:,1] = 1
  phi[:,2] = sqrt(3).*(s - t)
  phi[:,3] = sqrt(5).*(1 - 6.*s.*t)
  phi[:,4] = sqrt(7).*(1 - 10.*s.*t).*(s - t)
  phi[:,5] = 3.*(1 - 10.*s.*t.*(2 - 7.*s.*t))

  dphi[:,1,1] =  0
  dphi[:,2,1] =  2.*sqrt(3)
  dphi[:,3,1] =  6.*sqrt(5).*(s - t)
  dphi[:,4,1] = 12.*sqrt(7).*(1 - 5.*s.*t)
  dphi[:,5,1] = 60.*(1 - 7.*s.*t).*(s - t)

  # transpose
  phi   = phi'
  dphi  = dphi'

  return phi, dphi

end

function basisFuncLineLeg( ::Type{Val{5}}, s::Array{Float64} )

  dim = 1
  nx  = length( s )

  # Quintic
  nphi = 6
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  t = 1 - s

  phi[:,1] = 1;
  phi[:,2] = sqrt(3).*(s - t);
  phi[:,3] = sqrt(5).*(1 - 6.*s.*t);
  phi[:,4] = sqrt(7).*(1 - 10.*s.*t).*(s - t);
  phi[:,5] = 3.*(1 - 10.*s.*t.*(2 - 7.*s.*t));
  phi[:,6] = sqrt(11).*(1 - 14.*s.*t.*(s - 2.*t).*(2.*s - t)).*(s - t);


  dphi[:,1,1] =  0
  dphi[:,2,1] =  2.*sqrt(3)
  dphi[:,3,1] =  6.*sqrt(5).*(s - t)
  dphi[:,4,1] = 12.*sqrt(7).*(1 - 5.*s.*t)
  dphi[:,5,1] = 60.*(1 - 7.*s.*t).*(s - t)
  dphi[:,6,1] = 30.*sqrt(11).*(1 - 14.*s.*t.*(1 - 3.*s.*t))

  # transpose
  phi   = phi'
  dphi  = dphi'

  return phi, dphi

end
