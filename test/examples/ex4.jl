# Toy problem to test the solver

using luteos

mat = Material(E = 1, ν = 0.33)

P = P3() # Polynomial order of solution

println("Generate mesh")
@time mesh   = Mesh3D( "cube", P, N = 5)
println("Generate master")
@time master = Master3D( P )

# println("Compute Jacobians")
# @time compJacob!( mesh, master )

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
@time prob = Problem( "Example 3 - Poisson", funcS, [1,1,1,2,1,1], 1, [funcB, funcB, funcB, funcB, funcB, funcB] )

println("Solve problem")
@time (uhath, uh, qh, uhathTri ) = hdgSolveCD( master, mesh, mat, prob )

# write solution
println("Write solution")
@time writeTecplotCD( "blaCD3D.dat", prob, mesh, uh, qh )
