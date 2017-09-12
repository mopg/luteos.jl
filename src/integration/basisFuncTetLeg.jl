function basisFuncTetLeg( ::Type{Val{0}}, s::Array{Float64}, t::Array{Float64}, u::Array{Float64} )

  dim = 3
  nx  = length( s )

  # Constant

  nphi = 1
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  # value
  phi[:,1] = 1.0

  # derivative in x
  dphi[:,1,1] = 0.0
  # derivative in y
  dphi[:,1,2] = 0.0
  # derivative in z
  dphi[:,1,3] = 0.0

  # transpose
  phi   = phi'

  dphi2 = fill( 0.0, size(dphi,2), size(dphi,1), 3 )
  dphi2[:,:,1] = dphi[:,:,1]'
  dphi2[:,:,2] = dphi[:,:,2]'
  dphi2[:,:,3] = dphi[:,:,3]'

  return phi, dphi2

end

function basisFuncTetLeg( ::Type{Val{1}}, s::Array{Float64}, t::Array{Float64}, u::Array{Float64} )

  dim = 3
  nx  = length( s )

  # Linear

  nphi = 4
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  # value
  phi[:,1] = 1
  phi[:,2] = (sqrt(5/3.)).*(-1+(4).*(s))
  phi[:,3] = (sqrt(10/3.)).*(-1+s+(3).*(t))
  phi[:,4] = (sqrt(10)).*(-1+s+t+(2).*(u))

  # derivative in x
  dphi[:,1,1] = 0
  dphi[:,2,1] = (4).*(sqrt(5/3.))
  dphi[:,3,1] = sqrt(10/3.)
  dphi[:,4,1] = sqrt(10)

  # derivative in y
  dphi[:,1,2] = 0
  dphi[:,2,2] = 0
  dphi[:,3,2] = sqrt(30)
  dphi[:,4,2] = sqrt(10)

  # derivative in z
  dphi[:,1,3] = 0
  dphi[:,2,3] = 0
  dphi[:,3,3] = 0
  dphi[:,4,3] = (2).*(sqrt(10))

  # transpose
  phi   = phi'

  dphi2 = fill( 0.0, size(dphi,2), size(dphi,1), 3 )
  dphi2[:,:,1] = dphi[:,:,1]'
  dphi2[:,:,2] = dphi[:,:,2]'
  dphi2[:,:,3] = dphi[:,:,3]'

  return phi, dphi2

end


function basisFuncTetLeg( ::Type{Val{2}}, s::Array{Float64}, t::Array{Float64}, u::Array{Float64} )

  dim = 3
  nx  = length( s )

  # Quadratric
  nphi = 10
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  t2 = t .* t
  u2 = u .* u

  # value
  phi[:,1] = 1
  phi[:,2] = (sqrt(5/3.)).*(-1+(4).*(s))
  phi[:,3] = (sqrt(10/3.)).*(-1+s+(3).*(t))
  phi[:,4] = (sqrt(10)).*(-1+s+t+(2).*(u))
  phi[:,5] = (sqrt(7/3.)).*(1+(5).*((s).*(-2+(3).*(s))))
  phi[:,6] = (sqrt(14/3.)).*((-1+(6).*(s)).*(-1+s+(3).*(t)))
  phi[:,7] = (sqrt(14)).*((-1+(6).*(s)).*(-1+s+t+(2).*(u)))
  phi[:,8] = (sqrt(7)).*((-1+s).^2+(8).*((-1+s).*(t))+(10).*(t2))
  phi[:,9] = (sqrt(21)).*((-1+s+(5).*(t)).*(-1+s+t+(2).*(u)))
  phi[:,10] = (sqrt(35)).*((-1+s+t).^2+(6).*((-1+s+t).*(u))+(6).*(u2))

  # derivative in x
  dphi[:,1,1] = 0
  dphi[:,2,1] = (4).*(sqrt(5/3.))
  dphi[:,3,1] = sqrt(10/3.)
  dphi[:,4,1] = sqrt(10)
  dphi[:,5,1] = (10).*((sqrt(7/3.)).*(-1+(3).*(s)))
  dphi[:,6,1] = (sqrt(14/3.)).*(-7+(12).*(s)+(18).*(t))
  dphi[:,7,1] = (sqrt(14)).*(-7+(12).*(s)+(6).*(t)+(12).*(u))
  dphi[:,8,1] = (2).*((sqrt(7)).*(-1+s+(4).*(t)))
  dphi[:,9,1] = (2).*((sqrt(21)).*(-1+s+(3).*(t)+u))
  dphi[:,10,1] = (2).*((sqrt(35)).*(-1+s+t+(3).*(u)))

  # derivative in y
  dphi[:,1,2] = 0
  dphi[:,2,2] = 0
  dphi[:,3,2] = sqrt(30)
  dphi[:,4,2] = sqrt(10)
  dphi[:,5,2] = 0
  dphi[:,6,2] = (sqrt(42)).*(-1+(6).*(s))
  dphi[:,7,2] = (sqrt(14)).*(-1+(6).*(s))
  dphi[:,8,2] = (4).*((sqrt(7)).*(-2+(2).*(s)+(5).*(t)))
  dphi[:,9,2] = (2).*((sqrt(21)).*(-3+(3).*(s)+(5).*(t)+(5).*(u)))
  dphi[:,10,2] = (2).*((sqrt(35)).*(-1+s+t+(3).*(u)))

  # derivative in z
  dphi[:,1,3] = 0
  dphi[:,2,3] = 0
  dphi[:,3,3] = 0
  dphi[:,4,3] = (2).*(sqrt(10))
  dphi[:,5,3] = 0
  dphi[:,6,3] = 0
  dphi[:,7,3] = (2).*((sqrt(14)).*(-1+(6).*(s)))
  dphi[:,8,3] = 0
  dphi[:,9,3] = (2).*((sqrt(21)).*(-1+s+(5).*(t)))
  dphi[:,10,3] = (6).*((sqrt(35)).*(-1+s+t+(2).*(u)))

  # transpose
  phi   = phi'

  dphi2 = fill( 0.0, size(dphi,2), size(dphi,1), 3 )
  dphi2[:,:,1] = dphi[:,:,1]'
  dphi2[:,:,2] = dphi[:,:,2]'
  dphi2[:,:,3] = dphi[:,:,3]'

  return phi, dphi2

end

function basisFuncTetLeg( ::Type{Val{3}}, s::Array{Float64}, t::Array{Float64}, u::Array{Float64} )

  dim = 3
  nx  = length( s )

  # Cubic
  nphi = 20
  phi  = Array{Float64}(nx, nphi)
  dphi = Array{Float64}(nx, nphi, dim)

  s2 = s  .* s
  t2 = t  .* t
  u2 = u  .* u
  s3 = s2 .* s
  t3 = t2 .* t
  u3 = u2 .* u

  # value
  phi[:,1] = 1
  phi[:,2] = (sqrt(5/3.)).*(-1+(4).*(s))
  phi[:,3] = (sqrt(10/3.)).*(-1+s+(3).*(t))
  phi[:,4] = (sqrt(10)).*(-1+s+t+(2).*(u))
  phi[:,5] = (sqrt(7/3.)).*(1+(5).*((s).*(-2+(3).*(s))))
  phi[:,6] = (sqrt(14/3.)).*((-1+(6).*(s)).*(-1+s+(3).*(t)))
  phi[:,7] = (sqrt(14)).*((-1+(6).*(s)).*(-1+s+t+(2).*(u)))
  phi[:,8] = (sqrt(7)).*((-1+s).^2+(8).*((-1+s).*(t))+(10).*(t2))
  phi[:,9] = (sqrt(21)).*((-1+s+(5).*(t)).*(-1+s+t+(2).*(u)))
  phi[:,10] = (sqrt(35)).*((-1+s+t).^2+(6).*((-1+s+t).*(u))+(6).*(u2))
  phi[:,11] = (sqrt(3)).*(-1+(s).*(18+(7).*((s).*(-9+(8).*(s)))))
  phi[:,12] = (sqrt(6)).*((1+(14).*((s).*(-1+(2).*(s)))).*(-1+s+(3).*(t)))
  phi[:,13] = (3).*((sqrt(2)).*((1+(14).*((s).*(-1+(2).*(s)))).*(-1+s+t+(2).*(u))))
  phi[:,14] = (3).*((-1+(8).*(s)).*((-1+s).^2+(8).*((-1+s).*(t))+(10).*(t2)))
  phi[:,15] = (3).*((sqrt(3)).*((-1+(8).*(s)).*((-1+s+(5).*(t)).*(-1+s+t+(2).*(u)))))
  phi[:,16] = (3).*((sqrt(5)).*((-1+(8).*(s)).*((-1+s+t).^2+(6).*((-1+s+t).*(u))+(6).*(u2))))
  phi[:,17] = (2).*((sqrt(3)).*((-1+s).^3+(15).*(((-1+s).^2).*(t))+(45).*((-1+s).*(t2))+(35).*(t3)))
  phi[:,18] = (6).*(((-1+s).^2+(12).*((-1+s).*(t))+(21).*(t2)).*(-1+s+t+(2).*(u)))
  phi[:,19] = (2).*((sqrt(15)).*((-1+s+(7).*(t)).*((-1+s+t).^2+(6).*((-1+s+t).*(u))+(6).*(u2))))
  phi[:,20] = (2).*((sqrt(21)).*((-1+s+t).^3+(12).*(((-1+s+t).^2).*(u))+(30).*((-1+s+t).*(u2))+(20).*(u3)))

  # derivative in x
  dphi[:,1,1] = 0
  dphi[:,2,1] = (4).*(sqrt(5/3.))
  dphi[:,3,1] = sqrt(10/3.)
  dphi[:,4,1] = sqrt(10)
  dphi[:,5,1] = (10).*((sqrt(7/3.)).*(-1+(3).*(s)))
  dphi[:,6,1] = (sqrt(14/3.)).*(-7+(12).*(s)+(18).*(t))
  dphi[:,7,1] = (sqrt(14)).*(-7+(12).*(s)+(6).*(t)+(12).*(u))
  dphi[:,8,1] = (2).*((sqrt(7)).*(-1+s+(4).*(t)))
  dphi[:,9,1] = (2).*((sqrt(21)).*(-1+s+(3).*(t)+u))
  dphi[:,10,1] = (2).*((sqrt(35)).*(-1+s+t+(3).*(u)))
  dphi[:,11,1] = (6).*((sqrt(3)).*(3+(7).*((s).*(-3+(4).*(s)))))
  dphi[:,12,1] = (3).*((sqrt(6)).*(5+(-14).*(t)+(28).*((s).*(-1+s+(2).*(t)))))
  dphi[:,13,1] = (3).*((sqrt(2)).*(15+(-14).*(t)+(-28).*(u)+(28).*((s).*(-3+(3).*(s)+(2).*(t)+(4).*(u)))))
  dphi[:,14,1] = (6).*(5+(12).*(s2)+(4).*((t).*(-9+(10).*(t)))+(s).*(-17+(64).*(t)))
  dphi[:,15,1] = (6).*((sqrt(3)).*(5+(12).*(s2)+(-9).*(u)+(s).*(-17+(48).*(t)+(16).*(u))+(t).*(-27+(20).*(t)+(40).*(u))))
  dphi[:,16,1] = (6).*((sqrt(5)).*(5+(12).*(s2)+(4).*(t2)+(3).*((u).*(-9+(8).*(u)))+(3).*((t).*(-3+(8).*(u)))+(s).*(-17+(16).*(t)+(48).*(u))))
  dphi[:,17,1] = (6).*((sqrt(3)).*((-1+s).^2+(10).*((-1+s).*(t))+(15).*(t2)))
  dphi[:,18,1] = (6).*(3+(3).*(s2)+(-4).*(u)+(s).*(-6+(26).*(t)+(4).*(u))+(t).*(-26+(33).*(t)+(24).*(u)))
  dphi[:,19,1] = (6).*((sqrt(15)).*(1+s2+(5).*(t2)+(2).*((-2+u).*(u))+(s).*(-2+(6).*(t)+(4).*(u))+(2).*((t).*(-3+(8).*(u)))))
  dphi[:,20,1] = (6).*((sqrt(21)).*((-1+s+t).^2+(8).*((-1+s+t).*(u))+(10).*(u2)))

  # derivative in y
  dphi[:,1,2] = 0
  dphi[:,2,2] = 0
  dphi[:,3,2] = sqrt(30)
  dphi[:,4,2] = sqrt(10)
  dphi[:,5,2] = 0
  dphi[:,6,2] = (sqrt(42)).*(-1+(6).*(s))
  dphi[:,7,2] = (sqrt(14)).*(-1+(6).*(s))
  dphi[:,8,2] = (4).*((sqrt(7)).*(-2+(2).*(s)+(5).*(t)))
  dphi[:,9,2] = (2).*((sqrt(21)).*(-3+(3).*(s)+(5).*(t)+(5).*(u)))
  dphi[:,10,2] = (2).*((sqrt(35)).*(-1+s+t+(3).*(u)))
  dphi[:,11,2] = 0
  dphi[:,12,2] = (3).*((sqrt(6)).*(1+(14).*((s).*(-1+(2).*(s)))))
  dphi[:,13,2] = (3).*((sqrt(2)).*(1+(14).*((s).*(-1+(2).*(s)))))
  dphi[:,14,2] = (12).*((-1+(8).*(s)).*(-2+(2).*(s)+(5).*(t)))
  dphi[:,15,2] = (6).*((sqrt(3)).*((-1+(8).*(s)).*(-3+(3).*(s)+(5).*(t)+(5).*(u))))
  dphi[:,16,2] = (6).*((sqrt(5)).*((-1+(8).*(s)).*(-1+s+t+(3).*(u))))
  dphi[:,17,2] = (30).*((sqrt(3)).*((-1+s).^2+(6).*((-1+s).*(t))+(7).*(t2)))
  dphi[:,18,2] = (6).*(13+(13).*(s2)+(-24).*(u)+(s).*(-26+(66).*(t)+(24).*(u))+(3).*((t).*(-22+(21).*(t)+(28).*(u))))
  dphi[:,19,2] = (6).*((sqrt(15)).*(3+(3).*(s2)+(7).*(t2)+(2).*((u).*(-8+(7).*(u)))+(2).*((s).*(-3+(5).*(t)+(8).*(u)))+(2).*((t).*(-5+(14).*(u)))))
  dphi[:,20,2] = (6).*((sqrt(21)).*((-1+s+t).^2+(8).*((-1+s+t).*(u))+(10).*(u2)))

  # derivative in z
  dphi[:,1,3] = 0
  dphi[:,2,3] = 0
  dphi[:,3,3] = 0
  dphi[:,4,3] = (2).*(sqrt(10))
  dphi[:,5,3] = 0
  dphi[:,6,3] = 0
  dphi[:,7,3] = (2).*((sqrt(14)).*(-1+(6).*(s)))
  dphi[:,8,3] = 0
  dphi[:,9,3] = (2).*((sqrt(21)).*(-1+s+(5).*(t)))
  dphi[:,10,3] = (6).*((sqrt(35)).*(-1+s+t+(2).*(u)))
  dphi[:,11,3] = 0
  dphi[:,12,3] = 0
  dphi[:,13,3] = (6).*((sqrt(2)).*(1+(14).*((s).*(-1+(2).*(s)))))
  dphi[:,14,3] = 0
  dphi[:,15,3] = (6).*((sqrt(3)).*((-1+(8).*(s)).*(-1+s+(5).*(t))))
  dphi[:,16,3] = (18).*((sqrt(5)).*((-1+(8).*(s)).*(-1+s+t+(2).*(u))))
  dphi[:,17,3] = 0
  dphi[:,18,3] = (12).*((-1+s).^2+(12).*((-1+s).*(t))+(21).*(t2))
  dphi[:,19,3] = (12).*((sqrt(15)).*((-1+s+(7).*(t)).*(-1+s+t+(2).*(u))))
  dphi[:,20,3] = (24).*((sqrt(21)).*((-1+s+t).^2+(5).*((-1+s+t).*(u))+(5).*(u2)))

  # transpose
  phi   = phi'

  dphi2 = fill( 0.0, size(dphi,2), size(dphi,1), 3 )
  dphi2[:,:,1] = dphi[:,:,1]'
  dphi2[:,:,2] = dphi[:,:,2]'
  dphi2[:,:,3] = dphi[:,:,3]'

  return phi, dphi2

end
