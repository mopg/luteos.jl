# ---------------------------------------------------------------------------- #
#
#   master3D.jl
#
#   Type for 3D master element
#   Inherits from the abstract "master" type
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

type Master3D <: Master

  porder::Int64         # Polynomial order of mesh
  pgauss::Int64         # Polynomial order to be integrated exactly

  gpts::Array{Float64}    # Gauss points  3D
  gwts::Array{Float64}    # Gauss weights 3D

  gpts2D::Array{Float64}  # Gauss points  2D
  gwts2D::Array{Float64}  # Gauss weights 2D

  gpts1D::Array{Float64}  # Gauss points  1D
  gwts1D::Array{Float64}  # Gauss weights 1D

  ϕ::Array{Float64}     # Shape functions in 3D
  ∇ϕ::Array{Float64}    # Derivatives of shape functions in 3D

  ϕ2D::Array{Float64}   # Shape functions in 2D
  ∇ϕ2D::Array{Float64}  # Derivatives of shape functions in 2D

  ϕ1D::Array{Float64}   # Shape functions in 1D
  ∇ϕ1D::Array{Float64}  # Derivatives of shape functions in 1D

  perm::Array{Int64}    # Node numbers on faces

end

function Master3D( porder::Int64; pgauss::Int64 = 3*porder, typeb = "lag" )

  if porder > 3
    error(" P > 3 not supported for 3D ")
  end

  (go1D, go2D, go3D) = compOrder3D( pgauss )

  (gpts_,   gwts_)   = quadratureTet( Val{go3D} )
  (gpts2D_, gwts2D_) = quadratureTriangle( Val{go2D} )
  (gpts1D_, gwts1D_) = quadratureLine( Val{go1D} )

  if typeb == "lag"
    (ϕ_,   ∇ϕ_)   = basisFuncTetLag( Val{porder}, gpts_[:,1], gpts_[:,2],  gpts_[:,3] )
    (ϕ2D_, ∇ϕ2D_) = basisFuncTriangleLag( Val{porder}, gpts2D_[:,1], gpts2D_[:,2] )
    (ϕ1D_, ∇ϕ1D_) = basisFuncLineLag( Val{porder}, gpts1D_ )
  elseif typeb == "leg"
    (ϕ_,   ∇ϕ_)   = basisFuncTetLeg( Val{porder}, gpts_[:,1], gpts_[:,2],  gpts_[:,3] )
    (ϕ2D_, ∇ϕ2D_) = basisFuncTriangleLeg( Val{porder}, gpts2D_[:,1], gpts2D_[:,2] )
    (ϕ1D_, ∇ϕ1D_) = basisFuncLineLeg( Val{porder}, gpts1D_ )
  end

  perm_ = findPerm3D( porder )

  Master3D( porder, pgauss, gpts_, gwts_, gpts2D_, gwts2D_, gpts1D_, gwts1D_,
    ϕ_, ∇ϕ_, ϕ2D_, ∇ϕ2D_, ϕ1D_, ∇ϕ1D_, perm_ )

end

function compOrder3D( orderRq )

  # 1D

  if orderRq <= 1
    order1D = 1
  elseif orderRq <= 3
    order1D = 3
  elseif orderRq <= 5
    order1D = 5
  elseif orderRq <= 7
    order1D = 7
  elseif orderRq <= 9
    order1D = 9
  elseif orderRq <= 11
    order1D = 11
  elseif orderRq <= 13
    order1D = 13
  elseif orderRq <= 15
    order1D = 15
  elseif orderRq <= 17
    order1D = 17
  elseif orderRq <= 19
    order1D = 19
  elseif orderRq <= 21
    order1D = 21
  elseif orderRq <= 23
    order1D = 23
  end

  # 2D

  if orderRq <= 2
    order2D = 1
  elseif orderRq <= 3
    order2D = 2
  elseif orderRq <= 4
    order2D = 3
  elseif orderRq <= 5
    order2D = 4
  elseif orderRq <= 8
    order2D = 5
  # elseif orderRq <= 10
  #   order2D = 6
  # elseif orderRq <= 11
  #   order2D = 7
  elseif orderRq <= 12
    order2D = 8
  elseif orderRq <= 13
    order2D = 9
  elseif orderRq <= 14
    order2D = 10
  end

  # 3D

  if orderRq <= 2
    order3D = 1
  elseif orderRq <= 3
    order3D = 2
  elseif orderRq <= 5
    order3D = 3
  elseif orderRq <= 6
    order3D = 4
  elseif orderRq <= 7
    order3D = 5
  elseif orderRq <= 8
    order3D = 6
  elseif orderRq <= 9
    order3D = 7
  elseif orderRq <= 10
    order3D = 8
  elseif orderRq <= 11
    order3D = 9
  elseif orderRq <= 12
    order3D = 10
  elseif orderRq <= 13
    order3D = 11
  elseif orderRq <= 14
    order3D = 12
  end

  return order1D, order2D, order3D

end

function findPerm3D( p )

  # 1:  1 2 3
  # 2:  1 3 2
  # 3:  2 1 3
  # 4:  2 3 1
  # 5:  3 1 2
  # 6:  3 2 1

  orderind = [1 2 3;
              1 3 2;
              2 1 3;
              2 3 1;
              3 1 2;
              3 2 1]
  orderflip = [1, -1, -1, 1, 1, -1]

  n2d = Int( round( 0.5*(p+1)*(p+2) ) )

  indfaces = fill( 0::Int64, 4, n2d )

  perm = fill( 0::Int64, n2d, 4, 6 )

  if p == 1

    # 2 3 4 - 1 4 3  - 1 2 4 - 1 3 2 : First ind is the vertices of face which does not included current vertices

    indfaces = [2 3 4;
                1 4 3;
                1 2 4;
                1 3 2]

    for jj in 1:4, ii in 1:6
      perm[:,jj,ii] = indfaces[ jj, orderind[ii,:] ]
    end

  elseif p == 2

    indfaces = [2  3  4  5  6  7
                3  1  4  9  5  8
                1  2  4  6  9  10
                2  1  3  8  7  10]

    for jj in 1:4
      ifac = indfaces[jj,:]
      for ii in 1:6
        # corners
        ifac[1:3] = indfaces[ jj,     orderind[ii,:] ]
        # edges
        ifac[4:6] = indfaces[ jj, 3 + orderind[ii,:] ]

        perm[:,jj,ii] = ifac
      end
    end

  elseif p == 3

    indfaces = [2  3  4   5  6  7  8  9 10 17
                3  1  4  13 14  6  5 11 12 18
                1  2  4   8  7 14 13 15 16 19
                2  1  3  12 11 10  9 16 15 20]
    for jj in 1:4
      ifac = indfaces[jj,:]
      perm[:,jj,1] = ifac
      for ii in [4,5]
        ## not flipped direction
        # corners
        ifac[1:3] = indfaces[ jj, orderind[ii,:] ]
        # edges
        oind = [ orderind[ii,:]'; orderind[ii,:]' ]
        oind = 2*oind
        oind[1,:] -= 1
        oind = oind[:]
        ifac[4:9] = indfaces[ jj, 3 + oind ]
        # middle point unchanged

        perm[:,jj,ii] = ifac
      end
      for ii in [2,3,6]
        # flipped direction
        ifac[1:3] = indfaces[ jj, orderind[ii,:] ]
        # edges
        oind = [ orderind[ii,:]'; orderind[ii,:]' ]
        oind = 2*oind
        oind[2,:] -= 1
        oind = oind[:]
        ifac[4:9] = indfaces[ jj, 3 + oind ]
        # middle point unchanged

        perm[:,jj,ii] = ifac
      end

    end

  end

  return perm

end
