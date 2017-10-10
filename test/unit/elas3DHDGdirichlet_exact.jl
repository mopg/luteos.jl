using luteos
using SymPy

# let # limit scope

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

  u1func_org = lambdify( u1 )
  u2func_org = lambdify( u2 )
  u3func_org = lambdify( u3 )

  F1func_org  = lambdify( F[1] )
  F2func_org  = lambdify( F[2] )
  F3func_org  = lambdify( F[3] )

  σ11func_org = lambdify( σ[1,1] )
  σ12func_org = lambdify( σ[1,2] )
  σ13func_org = lambdify( σ[1,3] )
  σ22func_org = lambdify( σ[2,2] )
  σ23func_org = lambdify( σ[2,3] )
  σ33func_org = lambdify( σ[3,3] )

  ϵ11func_org = lambdify( ϵ[1,1] )
  ϵ12func_org = lambdify( ϵ[1,2] )
  ϵ13func_org = lambdify( ϵ[1,3] )
  ϵ22func_org = lambdify( ϵ[2,2] )
  ϵ23func_org = lambdify( ϵ[2,3] )
  ϵ33func_org = lambdify( ϵ[3,3] )

  # make sure you can call it with matrices
  u1func = x -> u1func_org(x[:,1], x[:,2], x[:,3])
  u2func = x -> u2func_org(x[:,1], x[:,2], x[:,3])
  u3func = x -> u3func_org(x[:,1], x[:,2], x[:,3])

  F1func = x -> F1func_org(x[:,1], x[:,2], x[:,3])
  F2func = x -> F2func_org(x[:,1], x[:,2], x[:,3])
  F3func = x -> F3func_org(x[:,1], x[:,2], x[:,3])

  σ11func = x -> σ11func_org(x[:,1], x[:,2], x[:,3])
  σ12func = x -> σ12func_org(x[:,1], x[:,2], x[:,3])
  σ13func = x -> σ13func_org(x[:,1], x[:,2], x[:,3])
  σ22func = x -> σ22func_org(x[:,1], x[:,2], x[:,3])
  σ23func = x -> σ23func_org(x[:,1], x[:,2], x[:,3])
  σ33func = x -> σ33func_org(x[:,1], x[:,2], x[:,3])

  ϵ11func = x -> ϵ11func_org(x[:,1], x[:,2], x[:,3])
  ϵ12func = x -> ϵ12func_org(x[:,1], x[:,2], x[:,3])
  ϵ13func = x -> ϵ13func_org(x[:,1], x[:,2], x[:,3])
  ϵ22func = x -> ϵ22func_org(x[:,1], x[:,2], x[:,3])
  ϵ23func = x -> ϵ23func_org(x[:,1], x[:,2], x[:,3])
  ϵ33func = x -> ϵ33func_org(x[:,1], x[:,2], x[:,3])

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

P = 1 # Polynomial order
N = 3 # Grid size

@time mesh   = Mesh3D( "cube", P, N = N)
@time master = Master3D( P )

@time compJacob!( mesh, master )

prob = Problem( @sprintf("Reg %i %i", P, N), source, bctype, 0,
        [funcB, funcB, funcB, funcB, funcB, funcB] )

(uhathTri, uh, σh, ϵh) = hdgSolveElas( master, mesh, mat, prob )

# end # limit scope
