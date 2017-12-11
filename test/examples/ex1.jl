# 3D Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.0)

P = P1() # Polynomial order of solution

mesh   = Mesh3D( "cube", P, N = 3 )#17)
master = Master3D( P )

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
prob = Problem( "Example 1", funcS,  [1,1,1,1,1,1], false,[funcB, funcB, funcB, funcB, funcB, funcBNH] )

@time (uhath, uh, σh ) = hdgSolveElas( master, mesh, mat, prob )

# write solution
writeTecplot( "elas3D.dat", prob, mesh, uh, σh )

@printf("No. elements: %i", size(σh,3))
