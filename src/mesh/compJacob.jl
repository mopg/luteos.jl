function compJacob!( mesh::Mesh2D, master::Master2D )

  ∂x∂ξ = master.∇ϕ[:,:,1]' * mesh.nodes[:,1,:]
  ∂x∂η = master.∇ϕ[:,:,2]' * mesh.nodes[:,1,:]
  ∂y∂ξ = master.∇ϕ[:,:,1]' * mesh.nodes[:,2,:]
  ∂y∂η = master.∇ϕ[:,:,2]' * mesh.nodes[:,2,:]

  jac = ∂x∂ξ.*∂y∂η - ∂x∂η.*∂y∂ξ
  jcw = diagm( master.gwts ) * jac

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
