# Toy problem to test the solver

using luteos

P = P3() # Polynomial order of solution

mesh   = Mesh2D( "square", P, N = 30, M = 30)
master = Master2D( P )

# source function
function funcS( p::Array{Float64} )
  return fill(1.0, size(p,1), 1)
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 1)
end

prob = CDR( "Example 3 - Poission", 1.0, [0.0,0.0], funcS, [1,1,1,1], true, [funcB, funcB, funcB, funcB] )

@time (uhath, uh, qh, uhathTri ) = hdgSolve( master, mesh, prob )

# write solution
writeTecplot( "blaCD.dat", prob, mesh, uh, qh )
