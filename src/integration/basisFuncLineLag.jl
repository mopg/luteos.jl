# ---------------------------------------------------------------------------- #
#
#   basisFuncLineLag.jl
#
#   Get Lagrangian basis functions for Line
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

function basisFuncLineLag( P::P1, s::Vector{Float64} )

  nx  = length( s )

  # Linear
  nphi = 2
  phi  = Array{Float64}(nphi, nx)
  dphi = Array{Float64}(nphi, nx)

  t = 1 - s

  phi[1,:] = t
  phi[2,:] = s

  dphi[1,:] =  -1.0
  dphi[2,:] =   1.0

  return phi, dphi

end


function basisFuncLineLag( P::P2, s::Vector{Float64} )

  nx  = length( s )

  # Quadratric
  nphi = 3
  phi  = Array{Float64}(nphi, nx)
  dphi = Array{Float64}(nphi, nx)

  t = 1 - s

  phi[1,:] =  s.*-3.0+(s.*s).*2.0+1.0
  phi[2,:] =  -s+(s.*s).*2.0
  phi[3,:] =  s.*4.0-(s.*s).*4.0

  dphi[1,:,1] =  s.*4.0-3.0
  dphi[2,:,1] =  s.*4.0-1.0
  dphi[3,:,1] =  s.*-8.0+4.0

  return phi, dphi

end

function basisFuncLineLag( P::P3, s::Vector{Float64} )

  nx  = length( s )

  # Cubic
  nphi = 4
  phi  = Array{Float64}(nphi, nx)
  dphi = Array{Float64}(nphi, nx)

  t = 1 - s

  phi[1,:] =  s.*(-1.1E1/2.0)+(s.*s).*9.0-(s.*s.*s).*(9.0/2.0)+1.0
  phi[2,:] =  s-(s.*s).*(9.0/2.0)+(s.*s.*s).*(9.0/2.0)
  phi[3,:] =  s.*9.0-(s.*s).*(4.5E1/2.0)+(s.*s.*s).*(2.7E1/2.0)
  phi[4,:] =  s.*(-9.0/2.0)+(s.*s).*1.8E1-(s.*s.*s).*(2.7E1/2.0)

  dphi[1,:,1] =  s.*1.8E1-(s.*s).*(2.7E1/2.0)-1.1E1/2.0
  dphi[2,:,1] =  s.*-9.0+(s.*s).*(2.7E1/2.0)+1.0
  dphi[3,:,1] =  s.*-4.5E1+(s.*s).*(8.1E1/2.0)+9.0
  dphi[4,:,1] =  s.*3.6E1-(s.*s).*(8.1E1/2.0)-9.0/2.0

  return phi, dphi

end

function basisFuncLineLag( P::P4, s::Vector{Float64} )

  nx  = length( s )

  # Quartic
  nphi = 5
  phi  = Array{Float64}(nphi, nx)
  dphi = Array{Float64}(nphi, nx)

  t = 1 - s

  phi[1,:] =  (1./3).*(s.*-2.5E1+(s.*s).*7.0E1-(s.*s.*s).*8.0E1+(s.*s.*s.*s).*3.2E1+3.0)
  phi[2,:] =  -(1./3).*(s.*3.0-(s.*s).*2.2E1+(s.*s.*s).*4.8E1-(s.*s.*s.*s).*3.2E1)
  phi[3,:] =  (1./3).*(s.*4.8E1-(s.*s).*2.08E2+(s.*s.*s).*2.88E2-(s.*s.*s.*s).*1.28E2)
  phi[4,:] =  -(1./3).*(s.*3.6E1-(s.*s).*2.28E2+(s.*s.*s).*3.84E2-(s.*s.*s.*s).*1.92E2)
  phi[5,:] =  (1./3).*(s.*1.6E1-(s.*s).*1.12E2+(s.*s.*s).*2.24E2-(s.*s.*s.*s).*1.28E2)

  dphi[1,:,1] =  (1./3).*(s.*1.4E2-(s.*s).*2.4E2+(s.*s.*s).*1.28E2-2.5E1)
  dphi[2,:,1] =  (1./3).*(s.*4.4E1-(s.*s).*1.44E2+(s.*s.*s).*1.28E2-3.0)
  dphi[3,:,1] =  -(1./3).*(s.*4.16E2-(s.*s).*8.64E2+(s.*s.*s).*5.12E2-4.8E1)
  dphi[4,:,1] =  (1./3).*(s.*4.56E2-(s.*s).*1.152E3+(s.*s.*s).*7.68E2-3.6E1)
  dphi[5,:,1] =  -(1./3).*(s.*2.24E2-(s.*s).*6.72E2+(s.*s.*s).*5.12E2-1.6E1)

  return phi, dphi

end

function basisFuncLineLag( P::P4, s::Vector{Float64} )

  nx  = length( s )

  # Quintic
  nphi = 6
  phi  = Array{Float64}(nphi, nx)
  dphi = Array{Float64}(nphi, nx)

  t = 1 - s

  phi[1,:] =  -(1./24).*(s.*2.74E2-(s.*s).*1.125E3+(s.*s.*s).*2.125E3-(s.*s.*s.*s).*1.875E3+(s.*s.*s.*s.*s).*6.25E2-2.4E1)
  phi[2,:] =  (1./24).*(s.*2.4E1-(s.*s).*2.5E2+(s.*s.*s).*8.75E2-(s.*s.*s.*s).*1.25E3+(s.*s.*s.*s.*s).*6.25E2)
  phi[3,:] =  (1./24).*(s.*6.0E2-(s.*s).*3.85E3+(s.*s.*s).*8.875E3-(s.*s.*s.*s).*8.75E3+(s.*s.*s.*s.*s).*3.125E3)
  phi[4,:] =  -(1./24).*(s.*6.0E2-(s.*s).*5.35E3+(s.*s.*s).*1.475E4-(s.*s.*s.*s).*1.625E4+(s.*s.*s.*s.*s).*6.25E3)
  phi[5,:] =  (1./24).*(s.*4.0E2-(s.*s).*3.9E3+(s.*s.*s).*1.225E4-(s.*s.*s.*s).*1.5E4+(s.*s.*s.*s.*s).*6.25E3)
  phi[6,:] =  -(1./24).*(s.*1.5E2-(s.*s).*1.525E3+(s.*s.*s).*5.125E3-(s.*s.*s.*s).*6.875E3+(s.*s.*s.*s.*s).*3.125E3)


  dphi[1,:,1] =  -(1./24).*(s.*-2.25E3+(s.*s).*6.375E3-(s.*s.*s).*7.5E3+(s.*s.*s.*s).*3.125E3+2.74E2)
  dphi[2,:,1] =  (1./24).*(s.*-5.0E2+(s.*s).*2.625E3-(s.*s.*s).*5.0E3+(s.*s.*s.*s).*3.125E3+2.4E1)
  dphi[3,:,1] =  (1./24).*(s.*-7.7E3+(s.*s).*2.6625E4-(s.*s.*s).*3.5E4+(s.*s.*s.*s).*1.5625E4+6.0E2)
  dphi[4,:,1] =  -(1./24).*(s.*-1.07E4+(s.*s).*4.425E4-(s.*s.*s).*6.5E4+(s.*s.*s.*s).*3.125E4+6.0E2)
  dphi[5,:,1] =  (1./24).*(s.*-7.8E3+(s.*s).*3.675E4-(s.*s.*s).*6.0E4+(s.*s.*s.*s).*3.125E4+4.0E2)
  dphi[6,:,1] =  -(1./24).*(s.*-3.05E3+(s.*s).*1.5375E4-(s.*s.*s).*2.75E4+(s.*s.*s.*s).*1.5625E4+1.5E2)

  return phi, dphi

end
