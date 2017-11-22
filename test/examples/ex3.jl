# Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.33)

P = P3() # Polynomial order of solution

mesh   = Mesh2D( "square", P, N = 5, M = 3)
master = Master2D( P )

# source function
function funcS( p::Array{Float64} )
  return fill(1.0, size(p,1), 1)
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 1)
end

prob = Problem( "Example 3 - Poission", funcS, [1,1,1,1], 1, [funcB, funcB, funcB, funcB] )

(uhath, uh, qh, uhathTri ) = hdgSolveCD( master, mesh, mat, prob )

# write solution
writeTecplotCD( "blaCD.dat", prob, mesh, uh, qh )
