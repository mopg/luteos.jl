# Toy problem to test the solver

using luteos

P = P1() # Polynomial order of solution

mesh   = Mesh2D( "square", P, N = 11)
master = Master2D( P )

mat = Material( E = 1, ν = 0.0, dim = mesh.dim )

# source function
function funcS( p::Array{Float64} )
  return fill(1.0, size(p,1), 2)
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 2)
end

prob = Elas( "Example 1", mat, funcS, [1,1,1,1], true, [funcB, funcB, funcB, funcB] )

(uhath, uh, σh ) = hdgSolve( master, mesh, prob )

# write solution
# writeTecplot( "bla.dat", prob, mesh, uh, σh )
