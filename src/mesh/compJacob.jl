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
    compJacob!( mesh::Mesh2D, master::Master2D )

Computes Jacobians for 2D meshes and saves them in `mesh`.
"""
function compJacob!( mesh::Mesh2D, master::Master2D )

  ∂x∂ξ = master.∇ϕ[:,:,1]' * mesh.nodes[:,1,:]
  ∂x∂η = master.∇ϕ[:,:,2]' * mesh.nodes[:,1,:]
  ∂y∂ξ = master.∇ϕ[:,:,1]' * mesh.nodes[:,2,:]
  ∂y∂η = master.∇ϕ[:,:,2]' * mesh.nodes[:,2,:]

  jac = ∂x∂ξ.*∂y∂η - ∂x∂η.*∂y∂ξ
  jcw = diagm( master.gwts ) * jac

  # [∂ξ∂x ∂η∂x; ∂ξ∂y ∂η∂y] = [∂x∂ξ ∂y∂ξ; ∂x∂η ∂y∂η]^-1

  ∂ξ∂x =  1./jac .* ∂y∂η
  ∂η∂x = -1./jac .* ∂y∂ξ
  ∂ξ∂y = -1./jac .* ∂x∂η
  ∂η∂y =  1./jac .* ∂x∂ξ

  ∂ξ∂x_vec = fill(0.0, size(∂x∂ξ,1), size(∂x∂ξ,2), 4)
  ∂ξ∂x_vec[:,:,1] = ∂ξ∂x
  ∂ξ∂x_vec[:,:,2] = ∂ξ∂y
  ∂ξ∂x_vec[:,:,3] = ∂η∂x
  ∂ξ∂x_vec[:,:,4] = ∂η∂y

  mesh.jcw  = jcw
  mesh.∂ξ∂x = ∂ξ∂x_vec

end

"""
    compJacobFace( mesh::Mesh3D, master::Master3D, el::Int64, face::Int64 )

Returns Jacobian on on the `face` in element `el` for a 3D mesh.
"""
function compJacobFace( mesh::Mesh3D, master::Master3D, el::Int64, face::Int64 )

  nod  = master.perm[ :, face, abs.(mesh.t2f[el,face+4]) ]

  p2d  = master.ϕ2D'  * mesh.nodes[nod,:,el]

  ∂x₁∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,1,el]
  ∂x₁∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,1,el]
  ∂x₂∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,2,el]
  ∂x₂∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,2,el]
  ∂x₃∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,3,el]
  ∂x₃∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,3,el]

  # cross product to find normal vector, normal = ∂x∂ξ₁ × ∂x∂ξ₂
  normal = hcat( ( ∂x₂∂ξ₁ .* ∂x₃∂ξ₂ - ∂x₂∂ξ₂ .* ∂x₃∂ξ₁ ),
              -(   ∂x₁∂ξ₁ .* ∂x₃∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₃∂ξ₁ ),
               (   ∂x₁∂ξ₁ .* ∂x₂∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₂∂ξ₁ ) )
  # normalize the normal vector
  jac    = sqrt.( normal[:,1].^2 + normal[:,2].^2 + normal[:,3].^2 )
  normal = normal ./ ( jac * [1,1,1]' )

  if mesh.t2f[el,face+4] < 0
    normal = -normal
  end

  jcw = master.gwts2D .* jac

  return (master.ϕ2D, p2d, nod, normal, jcw)

end

"""
    compJacobFace( mesh::Mesh2D, master::Master2D, el::Int64, face::Int64 )

Returns Jacobian on on the `face` in element `el` for a 2D mesh.
"""
function compJacobFace( mesh::Mesh2D, master::Master2D, el::Int64, face::Int64 )

  if mesh.t2f[el,face] < 0
    # face DOES NOT follow counter-clockwise rotation
    nod    = master.perm[:,face,2]
    rotdir = true
  else
    # face follows counter-clockwise rotation
    nod    = master.perm[:,face,1]
    rotdir = false
  end

  ∂x₁∂ξ₁ = master.∇ϕ1d' * mesh.nodes[nod,1,el]
  ∂x₂∂ξ₁ = master.∇ϕ1d' * mesh.nodes[nod,2,el]

  p1d  = master.ϕ1d'  * mesh.nodes[nod,:,el]

  jac  = sqrt( ∂x₁∂ξ₁.^2 + ∂x₂∂ξ₁.^2 )
  jcw  = master.gwts1d .* jac

  if rotdir
    normal = -[ ∂x₂∂ξ₁ -∂x₁∂ξ₁ ] ./ [jac jac]
  else
    normal =  [ ∂x₂∂ξ₁ -∂x₁∂ξ₁ ] ./ [jac jac]
  end

  tangent = [ normal[:,2] -normal[:,1] ]

  return (master.ϕ1d, p1d, nod, normal, jcw)

end

"""
    compJacob( ::Type{Val{2}}, ∇ϕ::Array{Float64,3}, gwts::Vector{Float64},
               nodes::Matrix{Float64} )

Returns Jacobian of a 2D mesh.
"""
function compJacob( ::Type{Val{2}}, ∇ϕ::Array{Float64,3}, gwts::Vector{Float64}, nodes::Matrix{Float64} )

  ∂x∂ξ = ∇ϕ[:,:,1]' * nodes[:,1]
  ∂x∂η = ∇ϕ[:,:,2]' * nodes[:,1]
  ∂y∂ξ = ∇ϕ[:,:,1]' * nodes[:,2]
  ∂y∂η = ∇ϕ[:,:,2]' * nodes[:,2]

  jac = ∂x∂ξ.*∂y∂η - ∂x∂η.*∂y∂ξ
  jcw = diagm( gwts ) * jac

  # [∂ξ∂x ∂η∂x; ∂ξ∂y ∂η∂y] = [∂x∂ξ ∂y∂ξ; ∂x∂η ∂y∂η]^-1

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

  ∂x₁∂ξ₁ = master.∇ϕ[:,:,1]' * mesh.nodes[:,1,:]
  ∂x₁∂ξ₂ = master.∇ϕ[:,:,2]' * mesh.nodes[:,1,:]
  ∂x₁∂ξ₃ = master.∇ϕ[:,:,3]' * mesh.nodes[:,1,:]
  ∂x₂∂ξ₁ = master.∇ϕ[:,:,1]' * mesh.nodes[:,2,:]
  ∂x₂∂ξ₂ = master.∇ϕ[:,:,2]' * mesh.nodes[:,2,:]
  ∂x₂∂ξ₃ = master.∇ϕ[:,:,3]' * mesh.nodes[:,2,:]
  ∂x₃∂ξ₁ = master.∇ϕ[:,:,1]' * mesh.nodes[:,3,:]
  ∂x₃∂ξ₂ = master.∇ϕ[:,:,2]' * mesh.nodes[:,3,:]
  ∂x₃∂ξ₃ = master.∇ϕ[:,:,3]' * mesh.nodes[:,3,:]

  jac = ∂x₁∂ξ₁.*∂x₂∂ξ₂.*∂x₃∂ξ₃ + ∂x₂∂ξ₁.*∂x₃∂ξ₂.*∂x₁∂ξ₃ + ∂x₃∂ξ₁.*∂x₁∂ξ₂.*∂x₂∂ξ₃ -
        ∂x₃∂ξ₁.*∂x₂∂ξ₂.*∂x₁∂ξ₃ - ∂x₂∂ξ₁.*∂x₁∂ξ₂.*∂x₃∂ξ₃ - ∂x₁∂ξ₁.*∂x₃∂ξ₂.*∂x₂∂ξ₃
  jcw = diagm( master.gwts ) * jac

  ∂ξ₁∂x₁ =  1./jac .* ( ∂x₂∂ξ₂ .* ∂x₃∂ξ₃ - ∂x₃∂ξ₂ .* ∂x₂∂ξ₃ )
  ∂ξ₁∂x₂ =  1./jac .* ( ∂x₁∂ξ₃ .* ∂x₃∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₃∂ξ₃ )
  ∂ξ₁∂x₃ =  1./jac .* ( ∂x₁∂ξ₂ .* ∂x₂∂ξ₃ - ∂x₁∂ξ₃ .* ∂x₂∂ξ₂ )

  ∂ξ₂∂x₁ =  1./jac .* ( ∂x₂∂ξ₃ .* ∂x₃∂ξ₁ - ∂x₂∂ξ₁ .* ∂x₃∂ξ₃ )
  ∂ξ₂∂x₂ =  1./jac .* ( ∂x₁∂ξ₁ .* ∂x₃∂ξ₃ - ∂x₁∂ξ₃ .* ∂x₃∂ξ₁ )
  ∂ξ₂∂x₃ =  1./jac .* ( ∂x₁∂ξ₁ .* ∂x₂∂ξ₃ - ∂x₁∂ξ₃ .* ∂x₂∂ξ₁ )

  ∂ξ₃∂x₁ =  1./jac .* ( ∂x₂∂ξ₁ .* ∂x₃∂ξ₂ - ∂x₂∂ξ₂ .* ∂x₃∂ξ₁ )
  ∂ξ₃∂x₂ =  1./jac .* ( ∂x₁∂ξ₂ .* ∂x₃∂ξ₁ - ∂x₁∂ξ₁ .* ∂x₃∂ξ₂ )
  ∂ξ₃∂x₃ =  1./jac .* ( ∂x₁∂ξ₁ .* ∂x₂∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₂∂ξ₁ )

  ∂ξ∂x = fill(0.0, size(∂x₁∂ξ₁,1), size(∂x₁∂ξ₁,2), 9)
  ∂ξ∂x[:,:,1] = ∂ξ₁∂x₁
  ∂ξ∂x[:,:,2] = ∂ξ₁∂x₂
  ∂ξ∂x[:,:,3] = ∂ξ₁∂x₃

  ∂ξ∂x[:,:,4] = ∂ξ₂∂x₁
  ∂ξ∂x[:,:,5] = ∂ξ₂∂x₂
  ∂ξ∂x[:,:,6] = ∂ξ₂∂x₃

  ∂ξ∂x[:,:,7] = ∂ξ₃∂x₁
  ∂ξ∂x[:,:,8] = ∂ξ₃∂x₂
  ∂ξ∂x[:,:,9] = ∂ξ₃∂x₃

  mesh.jcw  = jcw
  mesh.∂ξ∂x = ∂ξ∂x

end

"""
    compJacob( ::Type{Val{3}}, ∇ϕ::Array{Float64,3}, gwts::Vector{Float64},
               nodes::Matrix{Float64} )

Returns Jacobian of a 3D mesh.
"""
function compJacob( ::Type{Val{3}}, ∇ϕ::Array{Float64,3}, gwts::Vector{Float64}, nodes::Matrix{Float64} )

  # http://www.csun.edu/~lcaretto/me692/Coordinate%20transformations.pdf

  ∂x₁∂ξ₁ = ∇ϕ[:,:,1]' * nodes[:,1]
  ∂x₁∂ξ₂ = ∇ϕ[:,:,2]' * nodes[:,1]
  ∂x₁∂ξ₃ = ∇ϕ[:,:,3]' * nodes[:,1]
  ∂x₂∂ξ₁ = ∇ϕ[:,:,1]' * nodes[:,2]
  ∂x₂∂ξ₂ = ∇ϕ[:,:,2]' * nodes[:,2]
  ∂x₂∂ξ₃ = ∇ϕ[:,:,3]' * nodes[:,2]
  ∂x₃∂ξ₁ = ∇ϕ[:,:,1]' * nodes[:,3]
  ∂x₃∂ξ₂ = ∇ϕ[:,:,2]' * nodes[:,3]
  ∂x₃∂ξ₃ = ∇ϕ[:,:,3]' * nodes[:,3]

  jac = ∂x₁∂ξ₁.*∂x₂∂ξ₂.*∂x₃∂ξ₃ + ∂x₂∂ξ₁.*∂x₃∂ξ₂.*∂x₁∂ξ₃ + ∂x₃∂ξ₁.*∂x₁∂ξ₂.*∂x₂∂ξ₃ -
        ∂x₃∂ξ₁.*∂x₂∂ξ₂.*∂x₁∂ξ₃ - ∂x₂∂ξ₁.*∂x₁∂ξ₂.*∂x₃∂ξ₃ - ∂x₁∂ξ₁.*∂x₃∂ξ₂.*∂x₂∂ξ₃
  jcw = diagm( gwts ) * jac

  ∂ξ₁∂x₁ =  1./jac .* ( ∂x₂∂ξ₂ .* ∂x₃∂ξ₃ - ∂x₃∂ξ₂ .* ∂x₂∂ξ₃ )
  ∂ξ₁∂x₂ =  1./jac .* ( ∂x₁∂ξ₃ .* ∂x₃∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₃∂ξ₃ )
  ∂ξ₁∂x₃ =  1./jac .* ( ∂x₁∂ξ₂ .* ∂x₂∂ξ₃ - ∂x₁∂ξ₃ .* ∂x₂∂ξ₂ )

  ∂ξ₂∂x₁ =  1./jac .* ( ∂x₂∂ξ₃ .* ∂x₃∂ξ₁ - ∂x₂∂ξ₁ .* ∂x₃∂ξ₃ )
  ∂ξ₂∂x₂ =  1./jac .* ( ∂x₁∂ξ₁ .* ∂x₃∂ξ₃ - ∂x₁∂ξ₃ .* ∂x₃∂ξ₁ )
  ∂ξ₂∂x₃ =  1./jac .* ( ∂x₁∂ξ₁ .* ∂x₂∂ξ₃ - ∂x₁∂ξ₃ .* ∂x₂∂ξ₁ )

  ∂ξ₃∂x₁ =  1./jac .* ( ∂x₂∂ξ₁ .* ∂x₃∂ξ₂ - ∂x₂∂ξ₂ .* ∂x₃∂ξ₁ )
  ∂ξ₃∂x₂ =  1./jac .* ( ∂x₁∂ξ₂ .* ∂x₃∂ξ₁ - ∂x₁∂ξ₁ .* ∂x₃∂ξ₂ )
  ∂ξ₃∂x₃ =  1./jac .* ( ∂x₁∂ξ₁ .* ∂x₂∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₂∂ξ₁ )

  ∂ξ∂x = fill(0.0, size(∂x₁∂ξ₁,1), 9)
  ∂ξ∂x[:,:,1] = ∂ξ₁∂x₁
  ∂ξ∂x[:,:,2] = ∂ξ₁∂x₂
  ∂ξ∂x[:,:,3] = ∂ξ₁∂x₃

  ∂ξ∂x[:,:,4] = ∂ξ₂∂x₁
  ∂ξ∂x[:,:,5] = ∂ξ₂∂x₂
  ∂ξ∂x[:,:,6] = ∂ξ₂∂x₃

  ∂ξ∂x[:,:,7] = ∂ξ₃∂x₁
  ∂ξ∂x[:,:,8] = ∂ξ₃∂x₂
  ∂ξ∂x[:,:,9] = ∂ξ₃∂x₃

  return (jcw, ∂ξ∂x)

end
