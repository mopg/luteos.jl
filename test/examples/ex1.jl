# 3D Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.33)

P = 1 # Polynomial order of solution

mesh   = Mesh3D( "cube", P, N = 3)
master = Master3D( P )

compJacob!( mesh, master )

# source function
function funcS( p::Array{Float64} )
  return fill(1.0, size(p,1), 3)
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 3)
end

prob = Problem( "Example 1", funcS, [1,1,1,1,1,1], 0, [funcB, funcB, funcB, funcB, funcB, funcB] )

(uhath, uh, σh, ϵh ) = hdgSolveElas( master, mesh, mat, prob )

# write solution
# writeTecplot( "bla.dat", prob, mesh, uh, σh, ϵh )
