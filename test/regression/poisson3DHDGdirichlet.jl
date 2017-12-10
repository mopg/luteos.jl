using luteos
using SymPy

# let # limit scope

Ps = [P1(), P2()]#, P3()]  # Range of polynomial order
Ns = [9, 17]             # Range of grid size

dim = 3

# Material
mat = Material( ν=0.33, E=1.0 )

## Set up problem
#   Get functions for exact solution
function ExactSol3D(  )

  dim = 3

  @syms x1 x2 x3

  u = sin.(pi*x1) .* sin.(pi*x2) .* sin.(pi*x3)

  g = (x1^2 - x1) .* (x2^2 - x2) .* (x3^2 - x3)

  J = integrate( integrate( integrate( u .* g , x1, 0, 1 ), x2, 0, 1 ), x3, 0, 1 )

  ∂u∂x1 = diff( u, x1 )
  ∂u∂x2 = diff( u, x2 )
  ∂u∂x3 = diff( u, x3 )

  F  = diff( -∂u∂x1, x1 ) + diff( -∂u∂x2, x2 ) + diff( -∂u∂x3, x3 )

  ufunc_org = lambdify( u, [x1, x2, x3] )

  Ffunc_org  = lambdify( F, [x1, x2, x3] )

  q1func_org = lambdify( ∂u∂x1, [x1, x2, x3] )
  q2func_org = lambdify( ∂u∂x2, [x1, x2, x3] )
  q3func_org = lambdify( ∂u∂x3, [x1, x2, x3] )

  gfunc_org  = lambdify( g, [x1, x2, x3] )

  # make sure you can call it with matrices
  ufunc = x -> ufunc_org.(x[:,1], x[:,2], x[:,3])

  Ffunc = x -> Ffunc_org.(x[:,1], x[:,2], x[:,3])

  q1func = x -> q1func_org.(x[:,1], x[:,2], x[:,3])
  q2func = x -> q2func_org.(x[:,1], x[:,2], x[:,3])
  q3func = x -> q3func_org.(x[:,1], x[:,2], x[:,3])

  gfunc  = x -> gfunc_org.(x[:,1], x[:,2], x[:,3])

  return( ufunc,
          q1func,q2func,q3func,
          Ffunc, gfunc, J )

end

( ufunc, q1func, q2func, q3func, Ffunc, gfunc, J ) = ExactSol3D(  )

#   Setup boundary conditions
function funcB( p::Array{Float64} )
  return fill( 0.0, size(p,1) )
end
bctype = [1,1,1,1,1,1] # All Dirichlet
#   Setup source function
source = (p) -> Ffunc(p)

## Compute error for different polynomial orders and grid sizes
#   Initialize arrays
Err_uh  = fill( 0.0, length(Ps), length(Ns) )
Err_qh1 = fill( 0.0, length(Ps), length(Ns) )
Err_qh2 = fill( 0.0, length(Ps), length(Ns) )
Err_qh3 = fill( 0.0, length(Ps), length(Ns) )
Err_J   = fill( 0.0, length(Ps), length(Ns) )

#   Loop over polynomial order and grid size
for ii in 1:length(Ps), jj in 1:length(Ns)

  P = Ps[ii]; N = Ns[jj]

  @printf( " %6i %6i\n", P.p, N )

  mesh   = Mesh3D( "cube", P, N = N)
  master = Master3D( P )

  prob = Problem( @sprintf("Poisson - Reg %i %i", P.p, N), source, bctype, 0, [funcB, funcB, funcB, funcB, funcB, funcB] )

  (uhath, uh, qh, uhathTri) = hdgSolveCD( master, mesh, mat, prob )

  #   Initialize arrays
  err_uh = 0.0
  err_qh = fill( 0.0, dim )
  Jcomp  = 0.0

  # preallocate
  jcwd = fill( 0.0, size(master.∇ϕ,2), size(master.∇ϕ,2) )
  ∂ξ∂x = fill( 0.0, size(master.∇ϕ,2), dim^2 )
  ∂x∂ξ = fill( 0.0, size(master.∇ϕ,2), dim^2 )

  for kk in 1:size(mesh.t,1)

    # Compute Jacobians
    luteos.compJacob!( master, mesh.nodes[:,:,kk], ∂ξ∂x, jcwd, ∂x∂ξ )

    # u
    Δuh     = master.ϕ' * ( uh[:,1,kk] - ufunc( mesh.nodes[:,:,kk] ) )
    err_uh += Δuh' * jcwd * Δuh

    # σ
    Δqh1       = master.ϕ' * ( qh[:,1,kk] - q1func( mesh.nodes[:,:,kk] ) )
    Δqh2       = master.ϕ' * ( qh[:,2,kk] - q2func( mesh.nodes[:,:,kk] ) )
    Δqh3       = master.ϕ' * ( qh[:,3,kk] - q3func( mesh.nodes[:,:,kk] ) )

    # output functional
    gquad  = master.ϕ' * gfunc( mesh.nodes[:,:,kk] )
    uquad  = master.ϕ' * uh[:,1,kk]
    Jcomp += gquad' * jcwd * uquad

    err_qh[1] += Δqh1' * jcwd * Δqh1
    err_qh[2] += Δqh2' * jcwd * Δqh2
    err_qh[3] += Δqh3' * jcwd * Δqh3

  end

  Err_uh[ii,jj]  = sqrt(err_uh)

  Err_qh1[ii,jj] = sqrt(err_qh[1])
  Err_qh2[ii,jj] = sqrt(err_qh[2])
  Err_qh3[ii,jj] = sqrt(err_qh[3])

  Err_J[ii,jj]   = abs(Jcomp - J)

end

# Compute convergence rates
h = 1 ./ ( Ns - 1 )
conv_uh  = (log.( Err_uh[:,end-1]) - log.( Err_uh[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));

conv_qh1 = (log.( Err_qh1[:,end-1]) - log.( Err_qh1[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_qh2 = (log.( Err_qh2[:,end-1]) - log.( Err_qh2[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));
conv_qh3 = (log.( Err_qh3[:,end-1]) - log.( Err_qh3[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));

conv_J   = (log.( Err_J[:,end-1])   - log.( Err_J[:,end] ) ) / (log.( h[end-1]) - log.( h[end] ));

# Output to terminal
@printf("\n")
@printf("   Convergence rates for Poisson 3D Dirichlet problem\n\n")
@printf("   ---------------------------------------------------\n\n")
@printf( "P   ")
for jj in 1:size(Ps,1)
  @printf( " %6i", Ps[jj].p )
end
@printf( "\n" )

#   u
@printf( "u   ")
for jj in 1:size(Ps,1)
  @printf( " %6.4f", conv_uh[jj] )
end
@printf( "\n" )

#   q
@printf( "q₁  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_qh1[jj] )
end
@printf( "\n" )
@printf( "q₂  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_qh2[jj] )
end
@printf( "\n" )
@printf( "q₃  ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_qh3[jj] )
end
@printf( "\n" )

# J
@printf( "J   ")
for jj in 1:length(Ps)
  @printf( " %6.4f", conv_J[jj] )
end
@printf( "\n" )

open("errors_Poisson_Dirichlet3D.dat", "w") do f
  @printf(f, "P \t N \t E_uh \t E_qh1 \t E_qh2 \t E_qh3 \t E_J\n")
  for ii in 1:length(Ps), jj in 1:length(Ns)
    @printf(f, "%i \t %i \t %16.15e \t %16.15e \t %16.15e \t %16.15e \t %16.15e\n",
    Ps[ii].p,Ns[jj],Err_uh[ii,jj],Err_qh1[ii,jj],Err_qh2[ii,jj],Err_qh3[ii,jj],Err_J[ii,jj])
  end
end

# end # limit scope
