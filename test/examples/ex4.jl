# Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.33)

P = 3 # Polynomial order of solution

mesh   = Mesh3D( "cube", P, N = 5)
master = Master3D( P )

compJacob!( mesh, master )

# source function
function funcS( p::Array{Float64} )
  return -3*pi^2 * sin.(pi * (p[:,1]-0.5)) .* sin.(pi * (p[:,2]-0.5)) .* sin.(pi * (p[:,3]-0.5))
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 1)
end

function funcB2( p::Array{Float64} )
  return cos.( π * (p[:,1] - 0.5) ).^2 .* cos.( π * (p[:,3] - 0.5) ).^2
end

function funcU( p::Array{Float64} )
  return sin.(pi * (p[:,1]-0.5)) .* sin.(pi * (p[:,2]-0.5)) .* sin.(pi * (p[:,3]-0.5))
end

prob = Problem( "Example 3 - Poission", funcS, [1,1,1,1,1,1], 1, [funcB, funcB, funcB, funcB, funcB, funcB2] )

(uhath, uh, qh, uhathTri, A, B, N, D, H, M, K, L, C, E, R, G, F ) = hdgSolveCD( master, mesh, mat, prob )

# write solution
writeTecplotCD( "blaCD3D.dat", prob, mesh, uh, qh )


# NOTE 1: derivatives still look weird

# Derivative other method
∂x∂ξ  = fill( 0.0, size(master.∇ϕnod,2), size(mesh.nodes,3), 3, 3 )
∂ξ∂x2 = fill( 0.0, size(master.∇ϕnod,2), size(mesh.nodes,3), 3, 3 )

∂x∂ξ[:,:,1,1] = master.∇ϕnod[:,:,1]' * mesh.nodes[:,1,:]
∂x∂ξ[:,:,1,2] = master.∇ϕnod[:,:,2]' * mesh.nodes[:,1,:]
∂x∂ξ[:,:,1,3] = master.∇ϕnod[:,:,3]' * mesh.nodes[:,1,:]
∂x∂ξ[:,:,2,1] = master.∇ϕnod[:,:,1]' * mesh.nodes[:,2,:]
∂x∂ξ[:,:,2,2] = master.∇ϕnod[:,:,2]' * mesh.nodes[:,2,:]
∂x∂ξ[:,:,2,3] = master.∇ϕnod[:,:,3]' * mesh.nodes[:,2,:]
∂x∂ξ[:,:,3,1] = master.∇ϕnod[:,:,1]' * mesh.nodes[:,3,:]
∂x∂ξ[:,:,3,2] = master.∇ϕnod[:,:,2]' * mesh.nodes[:,3,:]
∂x∂ξ[:,:,3,3] = master.∇ϕnod[:,:,3]' * mesh.nodes[:,3,:]

jac = ∂x∂ξ[:,:,1,1].*∂x∂ξ[:,:,2,2].*∂x∂ξ[:,:,3,3] + ∂x∂ξ[:,:,2,1].*∂x∂ξ[:,:,3,2].*∂x∂ξ[:,:,1,3] + ∂x∂ξ[:,:,3,1].*∂x∂ξ[:,:,1,2].*∂x∂ξ[:,:,2,3] -
      ∂x∂ξ[:,:,3,1].*∂x∂ξ[:,:,2,2].*∂x∂ξ[:,:,1,3] - ∂x∂ξ[:,:,2,1].*∂x∂ξ[:,:,1,2].*∂x∂ξ[:,:,3,3] - ∂x∂ξ[:,:,1,1].*∂x∂ξ[:,:,3,2].*∂x∂ξ[:,:,2,3]

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

∂ξ∂x = fill(0.0, size(master.∇ϕnod,2), size(mesh.nodes,3), 9)
∂ξ∂x[:,:,1] = ∂ξ∂x2[:,:,1,1]#∂ξ₁∂x₁
∂ξ∂x[:,:,2] = ∂ξ∂x2[:,:,1,2]#∂ξ₁∂x₂
∂ξ∂x[:,:,3] = ∂ξ∂x2[:,:,1,3]#∂ξ₁∂x₃

∂ξ∂x[:,:,4] = ∂ξ∂x2[:,:,2,1]#∂ξ₂∂x₁
∂ξ∂x[:,:,5] = ∂ξ∂x2[:,:,2,2]#∂ξ₂∂x₂
∂ξ∂x[:,:,6] = ∂ξ∂x2[:,:,2,3]#∂ξ₂∂x₃

∂ξ∂x[:,:,7] = ∂ξ∂x2[:,:,3,1]#∂ξ₃∂x₁
∂ξ∂x[:,:,8] = ∂ξ∂x2[:,:,3,2]#∂ξ₃∂x₂
∂ξ∂x[:,:,9] = ∂ξ∂x2[:,:,3,3]#∂ξ₃∂x₃

qh2 = similar( qh )

for pp in 1:size(mesh.nodes,3) # Loop over all elements

    ∇ϕc = fill(0.0, size(master.∇ϕnod))

    ∇ϕc[:,:,1] = master.∇ϕnod[:,:,1] * diagm( ∂ξ∂x[:,pp,1] ) + master.∇ϕnod[:,:,2] * diagm( ∂ξ∂x[:,pp,4] ) + master.∇ϕnod[:,:,3] * diagm( ∂ξ∂x[:,pp,7] )
    ∇ϕc[:,:,2] = master.∇ϕnod[:,:,1] * diagm( ∂ξ∂x[:,pp,2] ) + master.∇ϕnod[:,:,2] * diagm( ∂ξ∂x[:,pp,5] ) + master.∇ϕnod[:,:,3] * diagm( ∂ξ∂x[:,pp,8] )
    ∇ϕc[:,:,3] = master.∇ϕnod[:,:,1] * diagm( ∂ξ∂x[:,pp,3] ) + master.∇ϕnod[:,:,2] * diagm( ∂ξ∂x[:,pp,6] ) + master.∇ϕnod[:,:,3] * diagm( ∂ξ∂x[:,pp,9] )

    qh2[:,1,pp] = ∇ϕc[:,:,1]' * uh[:,1,pp]
    qh2[:,2,pp] = ∇ϕc[:,:,2]' * uh[:,1,pp]
    qh2[:,3,pp] = ∇ϕc[:,:,3]' * uh[:,1,pp]

end

writeTecplotCD( "blaCD3D_v2.dat", prob, mesh, uh, qh2 )
