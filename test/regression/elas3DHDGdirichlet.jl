using luteos
using SymPy

# let # limit scope

Ps = 1:3         # Range of polynomial order
Ns = [7, 13]#[9, 17, 33] # Range of grid size

dim = 3

# Material
mat = Material( ν=0.33, E=1.0 )

## Set up problem
#   Get functions for exact solution
function ExactSol3D( Cstiff )

  dim = 3

  @syms x1 x2 x3

  u1 = 10 * (x2 - x2.^2) .* (x3 - x3.^2) .* sin.(pi*x1) .* (1 - x1) .* (1 - x2/2) .* (1 - x3/4)
  u2 =  2 * (x1 - x1.^2) .* (x3 - x3.^2) .* sin.(pi*x2) .* (1 - x2) .* (1 - x1/2) .* (1 - x3/4)
  u3 =  5 * (x1 - x1.^2) .* (x2 - x2.^2) .* sin.(pi*x3) .* (1 - x3) .* (1 - x1/2) .* (1 - x2/4)

  ∂u1∂x1 = diff( u1, x1 )
  ∂u1∂x2 = diff( u1, x2 )
  ∂u1∂x3 = diff( u1, x3 )
  ∂u2∂x1 = diff( u2, x1 )
  ∂u2∂x2 = diff( u2, x2 )
  ∂u2∂x3 = diff( u2, x3 )
  ∂u3∂x1 = diff( u3, x1 )
  ∂u3∂x2 = diff( u3, x2 )
  ∂u3∂x3 = diff( u3, x3 )

  ∇u = [ ∂u1∂x1 ∂u2∂x1 ∂u3∂x1; ∂u1∂x2 ∂u2∂x2 ∂u3∂x2; ∂u1∂x3 ∂u2∂x3 ∂u3∂x3 ]

  ϵ = 1/2 * ( ∇u + ∇u.' )

  σ = fill( 0*x1, dim, dim )

  for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
    σ[ii,jj] += (Cstiff[ii,jj,kk,ll] * ϵ[kk,ll])[1,1,1,1]
  end

  F  = fill( 0*x1, dim )

  ∇σ = fill( 0*x1, dim, dim, dim )

  for ii in 1:dim, jj in 1:dim
    ∇σ[ii,jj,1] = diff( σ[ii,jj], x1 )
    ∇σ[ii,jj,2] = diff( σ[ii,jj], x2 )
    ∇σ[ii,jj,3] = diff( σ[ii,jj], x3 )
  end

  for ii in 1:dim, jj in 1:dim
    F[ii] -= ∇σ[ii,jj,jj]
  end

  u1func_org = lambdify( u1, [x1, x2, x3] )
  u2func_org = lambdify( u2, [x1, x2, x3] )
  u3func_org = lambdify( u3, [x1, x2, x3] )

  F1func_org  = lambdify( F[1], [x1, x2, x3] )
  F2func_org  = lambdify( F[2], [x1, x2, x3] )
  F3func_org  = lambdify( F[3], [x1, x2, x3] )

  σ11func_org = lambdify( σ[1,1], [x1, x2, x3] )
  σ12func_org = lambdify( σ[1,2], [x1, x2, x3] )
  σ13func_org = lambdify( σ[1,3], [x1, x2, x3] )
  σ22func_org = lambdify( σ[2,2], [x1, x2, x3] )
  σ23func_org = lambdify( σ[2,3], [x1, x2, x3] )
  σ33func_org = lambdify( σ[3,3], [x1, x2, x3] )

  ϵ11func_org = lambdify( ϵ[1,1], [x1, x2, x3] )
  ϵ12func_org = lambdify( ϵ[1,2], [x1, x2, x3] )
  ϵ13func_org = lambdify( ϵ[1,3], [x1, x2, x3] )
  ϵ22func_org = lambdify( ϵ[2,2], [x1, x2, x3] )
  ϵ23func_org = lambdify( ϵ[2,3], [x1, x2, x3] )
  ϵ33func_org = lambdify( ϵ[3,3], [x1, x2, x3] )

  # make sure you can call it with matrices
  u1func = x -> u1func_org.(x[:,1], x[:,2], x[:,3])
  u2func = x -> u2func_org.(x[:,1], x[:,2], x[:,3])
  u3func = x -> u3func_org.(x[:,1], x[:,2], x[:,3])

  F1func = x -> F1func_org.(x[:,1], x[:,2], x[:,3])
  F2func = x -> F2func_org.(x[:,1], x[:,2], x[:,3])
  F3func = x -> F3func_org.(x[:,1], x[:,2], x[:,3])

  σ11func = x -> σ11func_org.(x[:,1], x[:,2], x[:,3])
  σ12func = x -> σ12func_org.(x[:,1], x[:,2], x[:,3])
  σ13func = x -> σ13func_org.(x[:,1], x[:,2], x[:,3])
  σ22func = x -> σ22func_org.(x[:,1], x[:,2], x[:,3])
  σ23func = x -> σ23func_org.(x[:,1], x[:,2], x[:,3])
  σ33func = x -> σ33func_org.(x[:,1], x[:,2], x[:,3])

  ϵ11func = x -> ϵ11func_org.(x[:,1], x[:,2], x[:,3])
  ϵ12func = x -> ϵ12func_org.(x[:,1], x[:,2], x[:,3])
  ϵ13func = x -> ϵ13func_org.(x[:,1], x[:,2], x[:,3])
  ϵ22func = x -> ϵ22func_org.(x[:,1], x[:,2], x[:,3])
  ϵ23func = x -> ϵ23func_org.(x[:,1], x[:,2], x[:,3])
  ϵ33func = x -> ϵ33func_org.(x[:,1], x[:,2], x[:,3])

  return( u1func, u2func, u3func,
          F1func, F2func, F3func,
          σ11func, σ12func, σ13func, σ22func, σ23func, σ33func,
          ϵ11func, ϵ12func, ϵ13func, ϵ22func, ϵ23func, ϵ33func )

end

( u1func, u2func, u3func,
  F1func, F2func, F3func,
  σ11func, σ12func, σ13func, σ22func, σ23func, σ33func,
  ϵ11func, ϵ12func, ϵ13func, ϵ22func, ϵ23func, ϵ33func ) = ExactSol3D( mat.Cstiff[dim] )

#   Setup boundary conditions
function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 3)
end
bctype = [1,1,1,1,1,1] # All Dirichlet
#   Setup source function
source = (p) -> [F1func(p) F2func(p) F3func(p)]

## Compute error for different polynomial orders and grid sizes
#   Initialize arrays
Err_uh1 = fill( 0.0, length(Ps), length(Ns) )
Err_uh2 = fill( 0.0, length(Ps), length(Ns) )
Err_uh3 = fill( 0.0, length(Ps), length(Ns) )
Err_σh1 = fill( 0.0, length(Ps), length(Ns) )
Err_σh2 = fill( 0.0, length(Ps), length(Ns) )
Err_σh3 = fill( 0.0, length(Ps), length(Ns) )
Err_σh5 = fill( 0.0, length(Ps), length(Ns) )
Err_σh6 = fill( 0.0, length(Ps), length(Ns) )
Err_σh9 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh1 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh2 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh3 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh5 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh6 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh9 = fill( 0.0, length(Ps), length(Ns) )

#   Loop over polynomial order and grid size
for ii in 1:length(Ps), jj in 1:length(Ns)

  P = Ps[ii]; N = Ns[jj]

  @printf( " %6i %6i\n", P, N )

  mesh   = Mesh3D( "cube", P, N = N)
  master = Master3D( P )

  prob = Problem( @sprintf("Reg %i %i", P, N), source, bctype, 0, [funcB, funcB, funcB, funcB, funcB, funcB] )

  (uhathTri, uh, σh, ϵh) = hdgSolveElas( master, mesh, mat, prob )

  #   Initialize arrays
  err_uh = fill( 0.0, dim )
  err_σh = fill( 0.0, dim^2 )
  err_ϵh = fill( 0.0, dim^2 )

  # preallocate
  jcw  = fill( 0.0, size(master.∇ϕ,2) )
  ∂ξ∂x = fill( 0.0, size(master.∇ϕ,2), dim^2 )
  ∂x∂ξ = fill( 0.0, size(master.∇ϕ,2), dim^2 )

  for kk in 1:size(mesh.t,1)

    # Compute Jacobians
    compJacob!( master, mesh.nodes[:,:,kk], ∂ξ∂x, jcw, ∂x∂ξ )

    jcwd = diagm( jcw )

    # u
    Δuh1       = master.ϕ' * ( uh[:,1,kk] - u1func( mesh.nodes[:,:,kk] ) )
    Δuh2       = master.ϕ' * ( uh[:,2,kk] - u2func( mesh.nodes[:,:,kk] ) )
    Δuh3       = master.ϕ' * ( uh[:,3,kk] - u3func( mesh.nodes[:,:,kk] ) )
    err_uh[1] += Δuh1' * jcwd * Δuh1
    err_uh[2] += Δuh2' * jcwd * Δuh2
    err_uh[3] += Δuh3' * jcwd * Δuh3

    # σ
    Δσh1       = master.ϕ' * ( σh[:,1,kk] - σ11func( mesh.nodes[:,:,kk] ) )
    Δσh2       = master.ϕ' * ( σh[:,2,kk] - σ12func( mesh.nodes[:,:,kk] ) )
    Δσh3       = master.ϕ' * ( σh[:,3,kk] - σ13func( mesh.nodes[:,:,kk] ) )
    Δσh5       = master.ϕ' * ( σh[:,5,kk] - σ22func( mesh.nodes[:,:,kk] ) )
    Δσh6       = master.ϕ' * ( σh[:,6,kk] - σ23func( mesh.nodes[:,:,kk] ) )
    Δσh9       = master.ϕ' * ( σh[:,9,kk] - σ33func( mesh.nodes[:,:,kk] ) )
    err_σh[1] += Δσh1' * jcwd * Δσh1
    err_σh[2] += Δσh2' * jcwd * Δσh2
    err_σh[3] += Δσh3' * jcwd * Δσh3
    err_σh[5] += Δσh5' * jcwd * Δσh5
    err_σh[6] += Δσh6' * jcwd * Δσh6
    err_σh[9] += Δσh9' * jcwd * Δσh9

    # ϵ
    Δϵh1       = master.ϕ' * ( ϵh[:,1,kk] - ϵ11func( mesh.nodes[:,:,kk] ) )
    Δϵh2       = master.ϕ' * ( ϵh[:,2,kk] - ϵ12func( mesh.nodes[:,:,kk] ) )
    Δϵh3       = master.ϕ' * ( ϵh[:,3,kk] - ϵ13func( mesh.nodes[:,:,kk] ) )
    Δϵh5       = master.ϕ' * ( ϵh[:,5,kk] - ϵ22func( mesh.nodes[:,:,kk] ) )
    Δϵh6       = master.ϕ' * ( ϵh[:,6,kk] - ϵ23func( mesh.nodes[:,:,kk] ) )
    Δϵh9       = master.ϕ' * ( ϵh[:,9,kk] - ϵ33func( mesh.nodes[:,:,kk] ) )
    err_ϵh[1] += Δϵh1' * jcwd * Δϵh1
    err_ϵh[2] += Δϵh2' * jcwd * Δϵh2
    err_ϵh[3] += Δϵh3' * jcwd * Δϵh3
    err_ϵh[5] += Δϵh5' * jcwd * Δϵh5
    err_ϵh[6] += Δϵh6' * jcwd * Δϵh6
    err_ϵh[9] += Δϵh9' * jcwd * Δϵh9

  end

  Err_uh1[ii,jj] = sqrt(err_uh[1])
  Err_uh2[ii,jj] = sqrt(err_uh[2])
  Err_uh3[ii,jj] = sqrt(err_uh[3])

  Err_σh1[ii,jj] = sqrt(err_σh[1])
  Err_σh2[ii,jj] = sqrt(err_σh[2])
  Err_σh3[ii,jj] = sqrt(err_σh[3])
  Err_σh5[ii,jj] = sqrt(err_σh[5])
  Err_σh6[ii,jj] = sqrt(err_σh[6])
  Err_σh9[ii,jj] = sqrt(err_σh[9])

  Err_ϵh1[ii,jj] = sqrt(err_ϵh[1])
  Err_ϵh2[ii,jj] = sqrt(err_ϵh[2])
  Err_ϵh3[ii,jj] = sqrt(err_ϵh[3])
  Err_ϵh5[ii,jj] = sqrt(err_ϵh[5])
  Err_ϵh6[ii,jj] = sqrt(err_ϵh[6])
  Err_ϵh9[ii,jj] = sqrt(err_ϵh[9])

end

# Compute convergence rates
h = 1 ./ ( Ns - 1 )
conv_uh1 = (log.( Err_uh1[:,end-1]) - log.( Err_uh1[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_uh2 = (log.( Err_uh2[:,end-1]) - log.( Err_uh2[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_uh3 = (log.( Err_uh3[:,end-1]) - log.( Err_uh3[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));

conv_σh1 = (log.( Err_σh1[:,end-1]) - log.( Err_σh1[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_σh2 = (log.( Err_σh2[:,end-1]) - log.( Err_σh2[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_σh3 = (log.( Err_σh3[:,end-1]) - log.( Err_σh3[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_σh5 = (log.( Err_σh5[:,end-1]) - log.( Err_σh5[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_σh6 = (log.( Err_σh6[:,end-1]) - log.( Err_σh6[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_σh9 = (log.( Err_σh9[:,end-1]) - log.( Err_σh9[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));

conv_ϵh1 = (log.( Err_ϵh1[:,end-1]) - log.( Err_ϵh1[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_ϵh2 = (log.( Err_ϵh2[:,end-1]) - log.( Err_ϵh2[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_ϵh3 = (log.( Err_ϵh3[:,end-1]) - log.( Err_ϵh3[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_ϵh5 = (log.( Err_ϵh5[:,end-1]) - log.( Err_ϵh5[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_ϵh6 = (log.( Err_ϵh6[:,end-1]) - log.( Err_ϵh6[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_ϵh9 = (log.( Err_ϵh9[:,end-1]) - log.( Err_ϵh9[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));

# Output to terminal
@printf("\n")
@printf("   Convergence rates for 3D Dirichlet problem\n\n")
@printf("   ------------------------------------------\n\n")
@printf( "P    ")
for jj in 1:size(Ps,1)
  @printf( " %6i", Ps[jj] )
end
@printf( "\n" )

#   u
@printf( "u₁   ")
for jj in 1:size(Ps,1)
  @printf( " %6.4f", conv_uh1[jj] )
end
@printf( "\n" )
@printf( "u₂   ")
for jj in 1:size(Ps,1)
  @printf( " %6.4f", conv_uh2[jj] )
end
@printf( "\n" )
@printf( "u₃   ")
for jj in 1:size(Ps,1)
  @printf( " %6.4f", conv_uh3[jj] )
end
@printf( "\n" )

#   σ
@printf( "σ₁₁  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh1[jj] )
end
@printf( "\n" )
@printf( "σ₁₂  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh2[jj] )
end
@printf( "\n" )
@printf( "σ₁₃  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh3[jj] )
end
@printf( "\n" )
@printf( "σ₂₂  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh5[jj] )
end
@printf( "\n" )
@printf( "σ₂₃  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh6[jj] )
end
@printf( "\n" )
@printf( "σ₃₃  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh9[jj] )
end
@printf( "\n" )

#   ϵ
@printf( "ϵ₁₁  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh1[jj] )
end
@printf( "\n" )
@printf( "ϵ₁₂  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh2[jj] )
end
@printf( "\n" )
@printf( "ϵ₁₃  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh3[jj] )
end
@printf( "\n" )
@printf( "ϵ₂₂  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh5[jj] )
end
@printf( "\n" )
@printf( "ϵ₂₃  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh6[jj] )
end
@printf( "\n" )
@printf( "ϵ₃₃  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh9[jj] )
end
@printf( "\n" )

open("errors_Elas_Dirichlet3D.dat", "w") do f
  @printf(f, "P \t N \t E_uh1 \t E_uh2 \t E_uh3 \t E_σh1 \t E_σh2 \t E_σh3 \t E_σh5 \t E_σh6 \t E_σh9 \t E_ϵh1 \t E_ϵh2 \t E_ϵh3 \t E_ϵh5 \t E_ϵh6 \t E_ϵh9 \n")#\t E_J\n")
  for ii in 1:length(Ps), jj in 1:length(Ns)
    @printf(f, "%i \t %i \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e\n",
      Ps[ii], Ns[jj], Err_uh1[ii,jj], Err_uh2[ii,jj], Err_uh3[ii,jj],
      Err_σh1[ii,jj], Err_σh2[ii,jj], Err_σh3[ii,jj], Err_σh5[ii,jj], Err_σh6[ii,jj], Err_σh9[ii,jj],
      Err_ϵh1[ii,jj], Err_ϵh2[ii,jj], Err_ϵh3[ii,jj], Err_ϵh5[ii,jj], Err_ϵh6[ii,jj], Err_ϵh9[ii,jj] )
  end
end

# end # limit scope
