# Toy problem to test the solver

using luteos

P = P1() # Polynomial order of solution

println("Generate mesh")
@time mesh   = Mesh3D( "cube", P, N = 13)#9)
println("Generate master")
@time master = Master3D( P )

# source function
function funcS( p::Array{Float64} )
  return -3*pi^2 * sin.(pi * (p[:,1]-0.5)) .* sin.(pi * (p[:,2]-0.5)) .* sin.(pi * (p[:,3]-0.5))
end

function funcB( p::Array{Float64} )
  return fill(0.0, size(p,1), 1)
end

function funcB2( p::Array{Float64} )
  return cos.( π * (p[:,1] - 0.5) ).^2 .* cos.( π * (p[:,3] - 0.5) ).^2
end

function funcU( p::Array{Float64} )
  return sin.(pi * (p[:,1]-0.5)) .* sin.(pi * (p[:,2]-0.5)) .* sin.(pi * (p[:,3]-0.5))
end

println("Setup problem")
@time prob = Problem( "Example 4 - Poisson", 1.0, [0.0,0.0,0.0], funcS, [1,1,1,2,1,1], true, [funcB, funcB, funcB, funcB, funcB, funcB] )

println("Solve problem")
@time (uhath, uh, qh, uhathTri ) = hdgSolve( master, mesh, prob )

# write solution
println("Write solution")
@time writeTecplot( "blaCD3D.dat", prob, mesh, uh, qh )
