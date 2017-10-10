# Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.33)

P = 3 # Polynomial order of solution

mesh   = Mesh3D( "cube", P, N = 3)
master = Master3D( P )

compJacob!( mesh, master )

# source function
function funcS( p::Array{Float64} )
  return -3*pi^2 * sin.(pi * p[:,1]) .* sin.(pi * p[:,2]) .* sin.(pi * p[:,3])
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 1)
end

function funcU( p::Array{Float64} )
  return sin.(pi * p[:,1]) .* sin.(pi * p[:,2]) .* sin.(pi * p[:,3])
end

prob = Problem( "Example 3 - Poission", funcS, [1,1,1,1,1,1], 1, [funcB, funcB, funcB, funcB, funcB, funcB] )

(uhath, uh, qh, uhathTri, A, B, N, D, H, M, K, L, C, E, R, G, F ) = hdgSolveCD( master, mesh, mat, prob )

# write solution
writeTecplotCD( "blaCD3D.dat", prob, mesh, uh, qh )


# NOTE 1: derivatives still look weird
