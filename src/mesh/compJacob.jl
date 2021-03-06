# ---------------------------------------------------------------------------- #
#
#   compJacob.jl
#
#   Several functions to compute jacobians of an element for both 2D and 3D
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    compJacob( master::Master2D, nodes::Matrix{Float64},
               ∂ξ∂x::Matrix{Float64}, jcw::Vector{Float64},
               ∂x∂ξ::Matrix{Float64})

Returns Jacobian of a 2D cell.
"""
function compJacob!( master::Master2D, nodes::Matrix{Float64},
                     ∂ξ∂x::Matrix{Float64}, jcw::Matrix{Float64},
                     ∂x∂ξ::Matrix{Float64} )

  ∂x∂ξ[:,1] = master.∇ϕ[:,:,1]' * nodes[:,1] # ∂x∂ξ
  ∂x∂ξ[:,2] = master.∇ϕ[:,:,2]' * nodes[:,1] # ∂x∂η
  ∂x∂ξ[:,3] = master.∇ϕ[:,:,1]' * nodes[:,2] # ∂y∂ξ
  ∂x∂ξ[:,4] = master.∇ϕ[:,:,2]' * nodes[:,2] # ∂y∂η

  jac = ∂x∂ξ[:,1].*∂x∂ξ[:,4] - ∂x∂ξ[:,2].*∂x∂ξ[:,3]

  for ii in 1:length(jac)
    jcw[ii,ii] = master.gwts[ii] * jac[ii]
  end

  ∂ξ∂x[:,1] =  1./jac .* ∂x∂ξ[:,4]
  ∂ξ∂x[:,2] = -1./jac .* ∂x∂ξ[:,2]
  ∂ξ∂x[:,3] = -1./jac .* ∂x∂ξ[:,3]
  ∂ξ∂x[:,4] =  1./jac .* ∂x∂ξ[:,1]

end

"""
    getderbfel( master::Master2D, ∂ξ∂x::Matrix{Float64} )

Returns ∇ϕ in the global coordinate system for 2D meshes.
"""

function getderbfel( master::Master2D, ∂ξ∂x::Matrix{Float64} )

  ∇ϕc = fill(0.0, size(master.∇ϕ))

  ∇ϕc[:,:,1] = master.∇ϕ[:,:,1] * diagm( ∂ξ∂x[:,1] ) + master.∇ϕ[:,:,2] * diagm( ∂ξ∂x[:,3] )
  ∇ϕc[:,:,2] = master.∇ϕ[:,:,1] * diagm( ∂ξ∂x[:,2] ) + master.∇ϕ[:,:,2] * diagm( ∂ξ∂x[:,4] )

  return ∇ϕc

end

"""
    compJacobFace( mesh::Mesh2D, master::Master2D, el::Int64, face::Int64 )

Returns Jacobian on on the `face` in element `el` for a 2D mesh.
"""
function compJacobFace( mesh::Mesh2D, master::Master2D, el::Int64, face::Int64 )

  indnod = 1
  rotdir = false

  if mesh.t2f[el,face] < 0
    # face DOES NOT follow counter-clockwise rotation
    indnod = 2
    rotdir = true
  end

  nod    = master.perm[:,face,indnod]

  ∂x₁∂ξ₁ = master.∇ϕ1d' * mesh.nodes[nod,1,el]
  ∂x₂∂ξ₁ = master.∇ϕ1d' * mesh.nodes[nod,2,el]

  p1d  = master.ϕ1d'  * mesh.nodes[nod,:,el]

  jac  = sqrt.( ∂x₁∂ξ₁.^2 + ∂x₂∂ξ₁.^2 )
  jcw  = master.gwts1d .* jac

  normal = [ ∂x₂∂ξ₁ -∂x₁∂ξ₁ ] ./ [jac jac]

  if rotdir
    normal *= -1.
  end

  tangent = [ normal[:,2] -normal[:,1] ]

  return (master.ϕ1d, p1d, nod, normal, jcw)

end

"""
    compJacob( master::Master2D, ∇ϕ::Array{Float64,3}, gwts::Vector{Float64},
               nodes::Matrix{Float64} )

Returns Jacobian of a 2D mesh.
"""
function compJacob( master::Master2D, nodes::Matrix{Float64} )

  ∂x∂ξ = master.∇ϕ[:,:,1]' * nodes[:,1]
  ∂x∂η = master.∇ϕ[:,:,2]' * nodes[:,1]
  ∂y∂ξ = master.∇ϕ[:,:,1]' * nodes[:,2]
  ∂y∂η = master.∇ϕ[:,:,2]' * nodes[:,2]

  jac = ∂x∂ξ.*∂y∂η - ∂x∂η.*∂y∂ξ
  jcw = diagm( master.gwts ) * jac

  ∂ξ∂x =  1./jac .* ∂y∂η
  ∂η∂x = -1./jac .* ∂y∂ξ
  ∂ξ∂y = -1./jac .* ∂x∂η
  ∂η∂y =  1./jac .* ∂x∂ξ

  ∂ξ∂x_vec = fill(0.0, size(∂x∂ξ,1), 4)
  ∂ξ∂x_vec[:,1] = ∂ξ∂x
  ∂ξ∂x_vec[:,2] = ∂ξ∂y
  ∂ξ∂x_vec[:,3] = ∂η∂x
  ∂ξ∂x_vec[:,4] = ∂η∂y

  return (jcw, ∂ξ∂x_vec)

end

"""
    compJacob!( mesh::Mesh3D, master::Master3D )

Computes Jacobians for 3D meshes and saves them in `mesh`.
"""
function compJacob!( mesh::Mesh3D, master::Master3D )

  # http://www.csun.edu/~lcaretto/me692/Coordinate%20transformations.pdf

  ∂x∂ξ = fill( 0.0, size(master.∇ϕ,2), size(mesh.nodes,3), 3, 3 )
  ∂ξ∂x2 = fill( 0.0, size(master.∇ϕ,2), size(mesh.nodes,3), 3, 3 )

  ∂x∂ξ[:,:,1,1] = master.∇ϕ[:,:,1]' * mesh.nodes[:,1,:]
  ∂x∂ξ[:,:,1,2] = master.∇ϕ[:,:,2]' * mesh.nodes[:,1,:]
  ∂x∂ξ[:,:,1,3] = master.∇ϕ[:,:,3]' * mesh.nodes[:,1,:]
  ∂x∂ξ[:,:,2,1] = master.∇ϕ[:,:,1]' * mesh.nodes[:,2,:]
  ∂x∂ξ[:,:,2,2] = master.∇ϕ[:,:,2]' * mesh.nodes[:,2,:]
  ∂x∂ξ[:,:,2,3] = master.∇ϕ[:,:,3]' * mesh.nodes[:,2,:]
  ∂x∂ξ[:,:,3,1] = master.∇ϕ[:,:,1]' * mesh.nodes[:,3,:]
  ∂x∂ξ[:,:,3,2] = master.∇ϕ[:,:,2]' * mesh.nodes[:,3,:]
  ∂x∂ξ[:,:,3,3] = master.∇ϕ[:,:,3]' * mesh.nodes[:,3,:]

  jac = ∂x∂ξ[:,:,1,1].*∂x∂ξ[:,:,2,2].*∂x∂ξ[:,:,3,3] + ∂x∂ξ[:,:,2,1].*∂x∂ξ[:,:,3,2].*∂x∂ξ[:,:,1,3] + ∂x∂ξ[:,:,3,1].*∂x∂ξ[:,:,1,2].*∂x∂ξ[:,:,2,3] -
        ∂x∂ξ[:,:,3,1].*∂x∂ξ[:,:,2,2].*∂x∂ξ[:,:,1,3] - ∂x∂ξ[:,:,2,1].*∂x∂ξ[:,:,1,2].*∂x∂ξ[:,:,3,3] - ∂x∂ξ[:,:,1,1].*∂x∂ξ[:,:,3,2].*∂x∂ξ[:,:,2,3]
  jcw = diagm( master.gwts ) * jac

  for ii in 1:3, jj in 1:3
    iip1 = ii + 1
    jjp1 = jj + 1
    iip2 = ii + 2
    jjp2 = jj + 2
    if jjp1 != 3
      jjp1 = rem(jjp1,3)
    end
    if jjp2 != 3
      jjp2 = rem(jjp2,3)
    end
    if iip1 != 3
      iip1 = rem(iip1,3)
    end
    if iip2 != 3
      iip2 = rem(iip2,3)
    end
    ∂ξ∂x2[:,:,ii,jj] = 1./jac .* ( ∂x∂ξ[:,:,jjp1,iip1].*∂x∂ξ[:,:,jjp2,iip2] - ∂x∂ξ[:,:,jjp1,iip2].*∂x∂ξ[:,:,jjp2,iip1] )
  end

  ∂ξ∂x = fill(0.0, size(master.∇ϕ,2), size(mesh.nodes,3), 9)
  ∂ξ∂x[:,:,1] = ∂ξ∂x2[:,:,1,1]#∂ξ₁∂x₁
  ∂ξ∂x[:,:,2] = ∂ξ∂x2[:,:,1,2]#∂ξ₁∂x₂
  ∂ξ∂x[:,:,3] = ∂ξ∂x2[:,:,1,3]#∂ξ₁∂x₃

  ∂ξ∂x[:,:,4] = ∂ξ∂x2[:,:,2,1]#∂ξ₂∂x₁
  ∂ξ∂x[:,:,5] = ∂ξ∂x2[:,:,2,2]#∂ξ₂∂x₂
  ∂ξ∂x[:,:,6] = ∂ξ∂x2[:,:,2,3]#∂ξ₂∂x₃

  ∂ξ∂x[:,:,7] = ∂ξ∂x2[:,:,3,1]#∂ξ₃∂x₁
  ∂ξ∂x[:,:,8] = ∂ξ∂x2[:,:,3,2]#∂ξ₃∂x₂
  ∂ξ∂x[:,:,9] = ∂ξ∂x2[:,:,3,3]#∂ξ₃∂x₃

  mesh.jcw  = jcw
  mesh.∂ξ∂x = ∂ξ∂x

end

"""
    getderbfel( master::Master3D, ∂ξ∂x::Matrix{Float64} )

Returns ∇ϕ in the global coordinate system for 3D meshes.
"""

function getderbfel( master::Master3D, ∂ξ∂x::Matrix{Float64} )

  ∇ϕc = fill(0.0, size(master.∇ϕ))

  ∇ϕc[:,:,1] = master.∇ϕ[:,:,1] * diagm( ∂ξ∂x[:,1] ) + master.∇ϕ[:,:,2] * diagm( ∂ξ∂x[:,4] ) + master.∇ϕ[:,:,3] * diagm( ∂ξ∂x[:,7] )
  ∇ϕc[:,:,2] = master.∇ϕ[:,:,1] * diagm( ∂ξ∂x[:,2] ) + master.∇ϕ[:,:,2] * diagm( ∂ξ∂x[:,5] ) + master.∇ϕ[:,:,3] * diagm( ∂ξ∂x[:,8] )
  ∇ϕc[:,:,3] = master.∇ϕ[:,:,1] * diagm( ∂ξ∂x[:,3] ) + master.∇ϕ[:,:,2] * diagm( ∂ξ∂x[:,6] ) + master.∇ϕ[:,:,3] * diagm( ∂ξ∂x[:,9] )

  return ∇ϕc

end

"""
    compJacobFace( mesh::Mesh3D, master::Master3D, el::Int64, face::Int64 )

Returns Jacobian on the `face` in element `el` for a 3D mesh.
"""
function compJacobFace( mesh::Mesh3D, master::Master3D, el::Int64, face::Int64 )

  nod  = master.perm[ :, face, abs.(mesh.t2f[el,face+4]) ]
  #nod  = master.perm[ :, 1, 1 ]

  p2d  = master.ϕ2D'  * mesh.nodes[nod,:,el]

  ∂x₁∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,1,el]
  ∂x₁∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,1,el]
  ∂x₂∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,2,el]
  ∂x₂∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,2,el]
  ∂x₃∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,3,el]
  ∂x₃∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,3,el]

  # cross product to find normal vector, normal = ∂x∂ξ₁ × ∂x∂ξ₂
  normal = hcat( ( ∂x₂∂ξ₁ .* ∂x₃∂ξ₂ - ∂x₂∂ξ₂ .* ∂x₃∂ξ₁ ),
                -( ∂x₁∂ξ₁ .* ∂x₃∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₃∂ξ₁ ),
                 ( ∂x₁∂ξ₁ .* ∂x₂∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₂∂ξ₁ ) )
  # normalize the normal vector
  jac      = sqrt.( normal[:,1].^2 + normal[:,2].^2 + normal[:,3].^2 )
  normal ./= jac * [1,1,1]'

  if mesh.t2f[el,face+4] < 0
    normal *= -1.0
  end

  jcw = master.gwts2D .* jac

  return (master.ϕ2D, p2d, nod, normal, jcw)

end

"""
    compJacob( master::Master3D, nodes::Matrix{Float64},
               ∂ξ∂x::Matrix{Float64}, jcw::Vector{Float64},
               ∂x∂ξ::Matrix{Float64})

Returns Jacobian of a 3D cell.
"""
function compJacob!( master::Master3D, nodes::Matrix{Float64},
                     ∂ξ∂x::Matrix{Float64}, jcw::Matrix{Float64},
                     ∂x∂ξ::Matrix{Float64} )

  # http://www.csun.edu/~lcaretto/me692/Coordinate%20transformations.pdf

  ∂x∂ξ[:,1] = master.∇ϕ[:,:,1]' * nodes[:,1] # 1 1 - 1
  ∂x∂ξ[:,2] = master.∇ϕ[:,:,2]' * nodes[:,1] # 1 2 - 2
  ∂x∂ξ[:,3] = master.∇ϕ[:,:,3]' * nodes[:,1] # 1 3 - 3
  ∂x∂ξ[:,4] = master.∇ϕ[:,:,1]' * nodes[:,2] # 2 1 - 4
  ∂x∂ξ[:,5] = master.∇ϕ[:,:,2]' * nodes[:,2] # 2 2 - 5
  ∂x∂ξ[:,6] = master.∇ϕ[:,:,3]' * nodes[:,2] # 2 3 - 6
  ∂x∂ξ[:,7] = master.∇ϕ[:,:,1]' * nodes[:,3] # 3 1 - 7
  ∂x∂ξ[:,8] = master.∇ϕ[:,:,2]' * nodes[:,3] # 3 2 - 8
  ∂x∂ξ[:,9] = master.∇ϕ[:,:,3]' * nodes[:,3] # 3 3 - 9

  jac = ∂x∂ξ[:,1].*∂x∂ξ[:,5].*∂x∂ξ[:,9] + ∂x∂ξ[:,4].*∂x∂ξ[:,8].*∂x∂ξ[:,3] + ∂x∂ξ[:,7].*∂x∂ξ[:,2].*∂x∂ξ[:,6] -
        ∂x∂ξ[:,7].*∂x∂ξ[:,5].*∂x∂ξ[:,3] - ∂x∂ξ[:,4].*∂x∂ξ[:,2].*∂x∂ξ[:,9] - ∂x∂ξ[:,1].*∂x∂ξ[:,8].*∂x∂ξ[:,6]

  for ii in 1:length(jac)
    jcw[ii,ii] = master.gwts[ii] * jac[ii]
  end

  for ii in 1:3, jj in 1:3
    kk = 3 * (ii - 1) + jj

    iip1 = ii + 1
    jjp1 = jj + 1

    iip2 = ii + 2
    jjp2 = jj + 2

    if jjp1 != 3
      jjp1 = rem(jjp1,3)
    end
    if jjp2 != 3
      jjp2 = rem(jjp2,3)
    end
    if iip1 != 3
      iip1 = rem(iip1,3)
    end
    if iip2 != 3
      iip2 = rem(iip2,3)
    end

    kk11 = 3 * (jjp1 - 1) + iip1
    kk22 = 3 * (jjp2 - 1) + iip2
    kk12 = 3 * (jjp1 - 1) + iip2
    kk21 = 3 * (jjp2 - 1) + iip1

    ∂ξ∂x[:,kk] = 1./jac .* ( ∂x∂ξ[:,kk11].*∂x∂ξ[:,kk22] - ∂x∂ξ[:,kk12].*∂x∂ξ[:,kk21] )
  end

end
