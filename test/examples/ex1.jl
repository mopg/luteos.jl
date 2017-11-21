# 3D Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.0)

P = 3 # Polynomial order of solution

mesh   = Mesh3D( "cube", P, N = 2)
master = Master3D( P )

compJacob!( mesh, master )

# source function
function funcS( p::Array{Float64} )
  return fill(1.0, size(p,1), 3)
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 3)
end
function funcBNH( p::Array{Float64} )
  return fill(1.0, size(p,1), 3)
end

prob = Problem( funcS, [funcB, funcB, funcB, funcB, funcB, funcB], [1,1,1,1,1,1],
                name="Example 1", bcnorm=false )

(uhath, uh, σh, ϵh, uhathTri ) = hdgSolveElas( master, mesh, mat, prob )

# write solution
writeTecplot( "elas3D.dat", prob, mesh, uh, σh, ϵh )
