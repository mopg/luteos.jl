# Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.33)

P = P3() # Polynomial order of solution

mesh   = Mesh2D( "square", P, N = 11)
master = Master2D( P )

# source function
function funcS( p::Array{Float64} )
  return fill(1.0, size(p,1), 2)
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 2)
end

prob = Problem( "Example 1", funcS, [1,1,1,1], 1, [funcB, funcB, funcB, funcB] )

(uhath, uh, σh, ϵh, uhathTri ) = hdgSolveElas( master, mesh, mat, prob )

# write solution
writeTecplot( "bla.dat", prob, mesh, uh, σh, ϵh )
