using luteos
using SymPy

# let # limit scope

Ps = [P1(), P2(), P3()] # Range of polynomial order
Ns = [9, 17, 33]        # Range of grid size

dim = 2

# Material
mat = Material( ν=0.33, E=1.0 )

## Set up problem
#   Get functions for exact solution
function ExactSol( Cstiff )

  dim = 2

  @syms x1 x2

  u1 = 10 * (x2 - x2.^2) .* sin.(pi*x1) .* (1 - x1) .* (1 - x2/2)
  u2 =  2 * (x1 - x1.^2) .* sin.(pi*x2) .* (1 - x2) .* (1 - x1/2)

  ∂u1∂x1 = diff( u1, x1 )
  ∂u1∂x2 = diff( u1, x2 )
  ∂u2∂x1 = diff( u2, x1 )
  ∂u2∂x2 = diff( u2, x2 )

  ∇u = [ ∂u1∂x1 ∂u2∂x1; ∂u1∂x2 ∂u2∂x2 ]

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
  end

  for ii in 1:dim, jj in 1:dim
    F[ii] -= ∇σ[ii,jj,jj]
  end

  u1func_org = lambdify( u1 )
  u2func_org = lambdify( u2 )

  F1func_org = lambdify( F[1] )
  F2func_org = lambdify( F[2] )
  σ1func_org = lambdify( σ[1] )
  σ2func_org = lambdify( σ[2] )
  σ4func_org = lambdify( σ[4] )
  ϵ1func_org = lambdify( ϵ[1] )
  ϵ2func_org = lambdify( ϵ[2] )
  ϵ4func_org = lambdify( ϵ[4] )

  # make sure you can call it with matrices
  u1func = x -> u1func_org.(x[:,1], x[:,2])
  u2func = x -> u2func_org.(x[:,1], x[:,2])

  F1func = x -> F1func_org.(x[:,1], x[:,2])
  F2func = x -> F2func_org.(x[:,1], x[:,2])
  σ1func = x -> σ1func_org.(x[:,1], x[:,2])
  σ2func = x -> σ2func_org.(x[:,1], x[:,2])
  σ4func = x -> σ4func_org.(x[:,1], x[:,2])
  ϵ1func = x -> ϵ1func_org.(x[:,1], x[:,2])
  ϵ2func = x -> ϵ2func_org.(x[:,1], x[:,2])
  ϵ4func = x -> ϵ4func_org.(x[:,1], x[:,2])

  return( u1func, u2func, F1func, F2func, σ1func, σ2func, σ4func, ϵ1func, ϵ2func, ϵ4func )

end

( u1func, u2func, F1func, F2func, σ1func, σ2func, σ4func, ϵ1func, ϵ2func, ϵ4func ) = ExactSol( mat.Cstiff[2] )

#   Setup boundary conditions
function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 2)
end
bctype = [1,1,1,1] # All Dirichlet
#   Setup source function
source = (p) -> [F1func(p) F2func(p)]

## Compute error for different polynomial orders and grid sizes
#   Initialize arrays
Err_uh1 = fill( 0.0, length(Ps), length(Ns) )
Err_uh2 = fill( 0.0, length(Ps), length(Ns) )
Err_σh1 = fill( 0.0, length(Ps), length(Ns) )
Err_σh2 = fill( 0.0, length(Ps), length(Ns) )
Err_σh4 = fill( 0.0, length(Ps), length(Ns) )
# Err_ϵh1 = fill( 0.0, length(Ps), length(Ns) )
# Err_ϵh2 = fill( 0.0, length(Ps), length(Ns) )
# Err_ϵh4 = fill( 0.0, length(Ps), length(Ns) )

#   Loop over polynomial order and grid size
for ii in 1:length(Ps), jj in 1:length(Ns)

  P = Ps[ii]; N = Ns[jj]

  mesh   = Mesh2D( "square", P, N = N)
  master = Master2D( P )

  prob = Problem( @sprintf("Reg %i %i", P.p, N), source, bctype, 0, [funcB, funcB, funcB, funcB] )

  (uhathTri, uh, σh) = hdgSolveElas( master, mesh, mat, prob )

  #   Initialize arrays
  err_uh = fill( 0.0, dim )
  err_σh = fill( 0.0, dim^2 )
  # err_ϵh = fill( 0.0, dim^2 )

  # preallocate
  jcwd = fill( 0.0, size(master.∇ϕ,2), size(master.∇ϕ,2) )
  ∂ξ∂x = fill( 0.0, size(master.∇ϕ,2), dim^2 )
  ∂x∂ξ = fill( 0.0, size(master.∇ϕ,2), dim^2 )

  for kk in 1:size(mesh.t,1)

    # Compute Jacobians
    luteos.compJacob!( master, mesh.nodes[:,:,kk], ∂ξ∂x, jcwd, ∂x∂ξ )

    # u
    Δuh1       = master.ϕ' * ( uh[:,1,kk] - u1func( mesh.nodes[:,:,kk] ) )
    Δuh2       = master.ϕ' * ( uh[:,2,kk] - u2func( mesh.nodes[:,:,kk] ) )
    err_uh[1] += Δuh1' * jcwd * Δuh1
    err_uh[2] += Δuh2' * jcwd * Δuh2

    # σ
    Δσh1       = master.ϕ' * ( σh[:,1,kk] - σ1func( mesh.nodes[:,:,kk] ) )
    Δσh2       = master.ϕ' * ( σh[:,2,kk] - σ2func( mesh.nodes[:,:,kk] ) )
    Δσh4       = master.ϕ' * ( σh[:,4,kk] - σ4func( mesh.nodes[:,:,kk] ) )
    err_σh[1] += Δσh1' * jcwd * Δσh1
    err_σh[2] += Δσh2' * jcwd * Δσh2
    err_σh[4] += Δσh4' * jcwd * Δσh4

    # # ϵ
    # Δϵh1       = master.ϕ' * ( ϵh[:,1,kk] - ϵ1func( mesh.nodes[:,:,kk] ) )
    # Δϵh2       = master.ϕ' * ( ϵh[:,2,kk] - ϵ2func( mesh.nodes[:,:,kk] ) )
    # Δϵh4       = master.ϕ' * ( ϵh[:,4,kk] - ϵ4func( mesh.nodes[:,:,kk] ) )
    # err_ϵh[1] += Δϵh1' * jcwd * Δϵh1
    # err_ϵh[2] += Δϵh2' * jcwd * Δϵh2
    # err_ϵh[4] += Δϵh4' * jcwd * Δϵh4

  end

  Err_uh1[ii,jj] = sqrt(err_uh[1])
  Err_uh2[ii,jj] = sqrt(err_uh[2])
  Err_σh1[ii,jj] = sqrt(err_σh[1])
  Err_σh2[ii,jj] = sqrt(err_σh[2])
  Err_σh4[ii,jj] = sqrt(err_σh[4])
  # Err_ϵh1[ii,jj] = sqrt(err_ϵh[1])
  # Err_ϵh2[ii,jj] = sqrt(err_ϵh[2])
  # Err_ϵh4[ii,jj] = sqrt(err_ϵh[4])

end

# Compute convergence rates
h = 1 ./ ( Ns - 1 )
conv_uh1 = (log.( Err_uh1[:,end-2]) - log.( Err_uh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_uh2 = (log.( Err_uh2[:,end-2]) - log.( Err_uh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

conv_σh1 = (log.( Err_σh1[:,end-2]) - log.( Err_σh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_σh2 = (log.( Err_σh2[:,end-2]) - log.( Err_σh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_σh4 = (log.( Err_σh4[:,end-2]) - log.( Err_σh4[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

# conv_ϵh1 = (log.( Err_ϵh1[:,end-2]) - log.( Err_ϵh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
# conv_ϵh2 = (log.( Err_ϵh2[:,end-2]) - log.( Err_ϵh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
# conv_ϵh4 = (log.( Err_ϵh4[:,end-2]) - log.( Err_ϵh4[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

# Output to terminal
@printf("\n")
@printf("   Convergence rates for 2D Dirichlet problem\n\n")
@printf("   ------------------------------------------\n\n")
@printf( "P   ")
for jj in 1:size(Ps,1)
  @printf( " %6i", Ps[jj].p )
end
@printf( "\n" )

#   u
@printf( "u₁  ")
for jj in 1:size(Ps,1)
  @printf( " %6.4f", conv_uh1[jj] )
end
@printf( "\n" )
@printf( "u₂  ")
for jj in 1:size(Ps,1)
  @printf( " %6.4f", conv_uh2[jj] )
end
@printf( "\n" )

#   σ
@printf( "σ₁  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh1[jj] )
end
@printf( "\n" )
@printf( "σ₂  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh2[jj] )
end
@printf( "\n" )
@printf( "σ₄  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_σh4[jj] )
end
@printf( "\n" )

# #   ϵ
# @printf( "ϵ₁  ")
# for jj in 1:length(Ps)
#   @printf( " %6.4f", conv_ϵh1[jj] )
# end
# @printf( "\n" )
# @printf( "ϵ₂  ")
# for jj in 1:length(Ps)
#   @printf( " %6.4f", conv_ϵh2[jj] )
# end
# @printf( "\n" )
# @printf( "ϵ₄  ")
# for jj in 1:length(Ps)
#   @printf( " %6.4f", conv_ϵh4[jj] )
# end
# @printf( "\n" )

open("errors_Elas_Dirichlet2D.dat", "w") do f
  @printf(f, "P \t N \t E_uh1 \t E_uh2 \t E_σh1 \t E_σh2 \t E_σh4 \n")#\t E_ϵh1 \t E_ϵh2 \t E_ϵh4 \n")#\t E_J\n")
  for ii in 1:length(Ps), jj in 1:length(Ns)
    @printf(f, "%i \t %i \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e \n",#\t %16.15e \t %16.15e \t %16.15e\n",
      Ps[ii].p, Ns[jj], Err_uh1[ii,jj], Err_uh2[ii,jj],
      Err_σh1[ii,jj], Err_σh2[ii,jj], Err_σh4[ii,jj] )#,
      #Err_ϵh1[ii,jj], Err_ϵh2[ii,jj], Err_ϵh4[ii,jj] )
  end
end

# end # limit scope
