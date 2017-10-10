tol = 1e-13

pmax = 4

### Lagrangian basis functions
println( "Lagrangian basis functions")
for ii = 1:pmax

  ploc,   = luteos.genlocal(ii)
  ploc1D = fill( 0.0, ii+1, 1 )
  for jj = 1:ii+1
    ploc1D[jj] = (jj-1)/ii
  end
  temp = ploc1D[2:ii]
  ploc1D[1:2] = ploc1D[ [1,ii+1] ]
  ploc1D[3:ii+1] = temp

  master = Master2D( ii; typeb = "lag" )

  phi1D, dphi1D = luteos.basisFuncLineLag( Val{ii}, ploc1D )
  sz    = size( phi1D )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(phi1D[jj,kk] - 1.0) < tol
    else
      @test abs(phi1D[jj,kk]) < tol
    end
  end

  phi, dphi = luteos.basisFuncTriangleLag( Val{ii}, ploc[:,2], ploc[:,3] )

  sz  = size( phi )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(phi[jj,kk] - 1.0) < tol
    else
      @test abs(phi[jj,kk]) < tol
    end
  end

end

# test interpolation

for ii = 1:pmax

  (ploc,) = luteos.genlocal( ii )

  (ϕ,∇ϕ) = luteos.basisFuncTriangleLag( Val{ii}, ploc[:,2], ploc[:,3] )

  sz  = size( ϕ )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(ϕ[jj,kk] - 1.0) < tol
    else
      @test abs(ϕ[jj,kk]) < tol
    end
  end

  # check evaluation for random point
  xloc   = [ 0.11698  0.2678;
             0.06905  0.1376 ]
  (ϕ,∇ϕ) = luteos.basisFuncTriangleLag( Val{ii}, xloc[:,1], xloc[:,2] )
  for jj in 1:2, dd in 1:2
    xloctest = ϕ[:,jj]' * ploc[:,dd+1]
    @test abs( xloctest - xloc[jj,dd] ) < tol

  end

end

# ### Legendre basis functions
# println( "Legendre basis functions")
# for ii = 1:pmax
#
#   @printf( "  Test 2D master element at %i\n", ii )
#
#   master = Master2D( ii; typeb = "leg" )
#
#   println( "      1D")
#
#   res1D = master.ϕ1d' * diagm(master.gwts1d) * master.ϕ1d
#   sz    = size( res1D )
#   for jj = 1:sz[1], kk = 1:sz[2]
#     if jj == kk
#       @test abs(res1D[jj,kk] - 1.0) < tol
#     else
#       @test abs(res1D[jj,kk]) < tol
#     end
#   end
#
#   println( "      2D")
#
#   res = master.ϕ' * diagm(master.gwts) * master.ϕ
#   sz  = size( res )
#   for jj = 1:sz[1], kk = 1:sz[2]
#     if jj == kk
#       @test abs(res[jj,kk] - 1.0) < tol
#     else
#       @test abs(res[jj,kk]) < tol
#     end
#   end
#
# end
