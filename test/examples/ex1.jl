# Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.33)

P = 1 # Polynomial order of solution

mesh   = Mesh2D( "square", P, N = 5)
master = Master2D( P, pgauss=5 )

compJacob!( mesh, master )

# source function
function func( p::Array{Float64} )
  return fill(1.0, size(p,1), 2)
end

function funcB( p::Array{Float64} )
  return fill(1.0, size(p,1), 2)
end

prob = Problem( "Example 1", func, [1,1,1,1], 1, [funcB, funcB, funcB, funcB] )

hdgSolveElas( master, mesh, mat, prob )
