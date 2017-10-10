tol = 1e-13

pmax = 3

for ii = 1:pmax

  (ploc,) = luteos.genlocal3D( ii )

  (ϕ,∇ϕ) = luteos.basisFuncTetLag( Val{ii}, ploc[:,2], ploc[:,3], ploc[:,4] )
  println( "      3D")

  sz  = size( ϕ )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(ϕ[jj,kk] - 1.0) < tol
    else
      @test abs(ϕ[jj,kk]) < tol
    end
  end

  # check evaluation for random point
  xloc   = [0.11698  0.2678  0.05263;
            0.06905  0.1376  0.08579]
  (ϕ,∇ϕ) = luteos.basisFuncTetLag( Val{ii}, xloc[:,1], xloc[:,2], xloc[:,3] )

  for jj in 1:2, dd in 1:3
    xloctest = ϕ[:,jj]' * ploc[:,dd+1]
    @test abs( xloctest - xloc[jj,dd] ) < tol

  end

end
