using luteos
using SymPy

# let # limit scope

Ps = 1:3         # Range of polynomial order
Ns = [9, 17, 33] # Range of grid size

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
  u1func = x -> u1func_org(x[:,1], x[:,2])
  u2func = x -> u2func_org(x[:,1], x[:,2])

  F1func = x -> F1func_org(x[:,1], x[:,2])
  F2func = x -> F2func_org(x[:,1], x[:,2])
  σ1func = x -> σ1func_org(x[:,1], x[:,2])
  σ2func = x -> σ2func_org(x[:,1], x[:,2])
  σ4func = x -> σ4func_org(x[:,1], x[:,2])
  ϵ1func = x -> ϵ1func_org(x[:,1], x[:,2])
  ϵ2func = x -> ϵ2func_org(x[:,1], x[:,2])
  ϵ4func = x -> ϵ4func_org(x[:,1], x[:,2])

  return( u1func, u2func, F1func, F2func, σ1func, σ2func, σ4func, ϵ1func, ϵ2func, ϵ4func )

end

( u1func, u2func, F1func, F2func, σ1func, σ2func, σ4func, ϵ1func, ϵ2func, ϵ4func ) = ExactSol( mat.Cstiff[2] )

#   Setup boundary conditions
function funcB( p::Array{Float64} )
  return fill(1.0, size(p,1), 2)
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
Err_ϵh1 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh2 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh4 = fill( 0.0, length(Ps), length(Ns) )

#   Loop over polynomial order and grid size
for jj in 1:length(Ns), ii in 1:length(Ps)

  P = Ps[ii]; N = Ns[jj]

  mesh   = Mesh2D( "square", P, N = N)
  master = Master2D( P )

  compJacob!( mesh, master )

  prob = Problem( @sprintf("Reg %i %i", P, N), source, bctype, 0, [funcB, funcB, funcB, funcB] )

  (uhathTri, uh, ϵh, σh) = hdgSolveElas( master, mesh, mat, prob )

  @printf( "%i %i %i\n", P, N, size( uh, 1 )*size( uh, 3 ) )

  #   Initialize arrays
  err_uh = fill( 0.0, dim )
  err_σh = fill( 0.0, dim^2 )
  err_ϵh = fill( 0.0, dim^2 )

  for kk in 1:size(mesh.t,1)

    jcwd = diagm(mesh.jcw[:,kk])

    # u
    diff_uh1   = master.ϕ' * ( uh[:,1,kk] - u1func( mesh.nodes[:,:,kk] ) )
    diff_uh2   = master.ϕ' * ( uh[:,2,kk] - u2func( mesh.nodes[:,:,kk] ) )
    err_uh[1] += diff_uh1' * jcwd * diff_uh1
    err_uh[1] += diff_uh2' * jcwd * diff_uh2

    # σ
    diff_σh1   = master.ϕ' * ( σh[:,1,kk] - σ1func( mesh.nodes[:,:,kk] ) )
    diff_σh2   = master.ϕ' * ( σh[:,2,kk] - σ2func( mesh.nodes[:,:,kk] ) )
    diff_σh4   = master.ϕ' * ( σh[:,4,kk] - σ4func( mesh.nodes[:,:,kk] ) )
    err_σh[1] += diff_σh1' * jcwd * diff_σh1
    err_σh[2] += diff_σh2' * jcwd * diff_σh2
    err_σh[4] += diff_σh4' * jcwd * diff_σh4

    # ϵ
    diff_ϵh1   = master.ϕ' * ( ϵh[:,1,kk] - ϵ1func( mesh.nodes[:,:,kk] ) )
    diff_ϵh2   = master.ϕ' * ( ϵh[:,2,kk] - ϵ2func( mesh.nodes[:,:,kk] ) )
    diff_ϵh4   = master.ϕ' * ( ϵh[:,4,kk] - ϵ4func( mesh.nodes[:,:,kk] ) )
    err_ϵh[1] += diff_ϵh1' * jcwd * diff_ϵh1
    err_ϵh[2] += diff_ϵh2' * jcwd * diff_ϵh2
    err_ϵh[4] += diff_ϵh4' * jcwd * diff_ϵh4

  end

  Err_uh1[ii,jj] += err_uh[1]
  Err_uh2[ii,jj] += err_uh[2]
  Err_σh1[ii,jj] += err_σh[1]
  Err_σh2[ii,jj] += err_σh[2]
  Err_σh4[ii,jj] += err_σh[4]
  Err_ϵh1[ii,jj] += err_ϵh[1]
  Err_ϵh2[ii,jj] += err_ϵh[2]
  Err_ϵh4[ii,jj] += err_ϵh[4]

end

# Compute convergence rates
h = 1 ./ ( Ns - 1 )
conv_uh1 = (log.( Err_uh1[:,end-2]) - log.( Err_uh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_uh2 = (log.( Err_uh2[:,end-2]) - log.( Err_uh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

conv_σh1 = (log.( Err_σh1[:,end-2]) - log.( Err_σh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_σh2 = (log.( Err_σh2[:,end-2]) - log.( Err_σh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_σh4 = (log.( Err_σh4[:,end-2]) - log.( Err_σh4[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

conv_ϵh1 = (log.( Err_ϵh1[:,end-2]) - log.( Err_ϵh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_ϵh2 = (log.( Err_ϵh2[:,end-2]) - log.( Err_ϵh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_ϵh4 = (log.( Err_ϵh4[:,end-2]) - log.( Err_ϵh4[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

# Output to terminal
@printf("   Convergence rates for 2D Dirichlet problem\n\n")
@printf( "P   ")
for jj in 1:size(Ps,1)
  @printf( " %6i", Ps[jj] )
end
@printf( "\n" )

#   u
@printf( "u₁  ")
for jj in 1:size(Ps,1)
  @printf( " %6.4f", conv_uh1[jj] )
end
@printf( "\n" )
@printf( "u₁  ")
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

#   ϵ
@printf( "ϵ₁  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh1[jj] )
end
@printf( "\n" )
@printf( "ϵ₂  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh2[jj] )
end
@printf( "\n" )
@printf( "ϵ₄  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_ϵh4[jj] )
end
@printf( "\n" )

# end # limit scope
