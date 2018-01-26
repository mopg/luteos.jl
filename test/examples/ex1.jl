# 3D Toy problem to test the solver

using luteos

P = P1() # Polynomial order of solution

mesh   = Mesh3D( "cube", P, N = 7 )#17)
master = Master3D( P )

mat = Material( E = 1, ν = 0.0, dim = mesh.dim )

# source function
function funcS( p::Array{Float64,2} )
  return fill(1.0, size(p,1), 3)
end

function funcB( p::Array{Float64,2} )
  return fill(0.0, size(p,1), 3)
end
function funcBNH( p::Array{Float64,2} )
  return fill(1.0, size(p,1), 3)
end

# prob = Problem( funcS, [funcB, funcB, funcB, funcB, funcB, funcB], [1,1,1,1,1,1],
#                 name="Example 1", bcnorm=false )
prob = Elas( "Example 1", mat, funcS,  [1,1,1,1,1,1], false,[funcB, funcB, funcB, funcB, funcB, funcBNH] )

@time (uhath, uh, σh ) = hdgSolve( master, mesh, prob )

# write solution
writeTecplot( "elas3D.dat", prob, mesh, uh, σh )

@printf("No. elements: %i", size(σh,3))
