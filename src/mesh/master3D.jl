include("../integration/quadratureLine.jl")
include("../integration/quadratureTriangle.jl")
include("../integration/quadratureTet.jl")
include("../integration/basisFuncLine.jl")
include("../integration/basisFuncTriangle.jl")
include("../integration/basisFuncTet.jl")

type Master3D

  porder::Int64         # Polynomial order of mesh
  pgauss::Int64         # Polynomial order to be integrated exactly

  gpts::Array{Float64}    # Gauss points  3D
  gwts::Array{Float64}    # Gauss weights 3D

  gpts2D::Array{Float64}  # Gauss points  2D
  gwts2D::Array{Float64}  # Gauss weights 2D

  gpts1D::Array{Float64}  # Gauss points  1D
  gwts1D::Array{Float64}  # Gauss weights 1D

  phi::Array{Float64}     # Shape functions in 3D
  dphi::Array{Float64}    # Derivatives of shape functions in 3D

  phi2D::Array{Float64}   # Shape functions in 2D
  dphi2D::Array{Float64}  # Derivatives of shape functions in 2D

  phi1D::Array{Float64}   # Shape functions in 1D
  dphi1D::Array{Float64}  # Derivatives of shape functions in 1D

end

function Master3D( porder::Int64, pgauss::Int64 = 3*porder )

  (go1D, go2D, go3D) = compOrder( pgauss )

  (gpts_,   gwts_)   = quadratureTet( Val{go3D} )
  (gpts2D_, gwts2D_) = quadratureTriangle( Val{go2D} )
  (gpts1D_, gwts1D_) = quadratureLine( Val{go1D} )

  (phi_,   dphi_)   = basisFuncTet( Val{porder}, gpts_[:,1], gpts_[:,2],  gpts_[:,3] )
  (phi2D_, dphi2D_) = basisFuncTriangle( Val{porder}, gpts2D_[:,1], gpts2D_[:,2] )
  (phi1D_, dphi1D_) = basisFuncLine( Val{porder}, gpts1D_ )

  Master3D( porder, pgauss, gpts_, gwts_, gpts2D_, gwts2D_, gpts1D_, gwts1D_,
    phi_, dphi_, phi2D_, dphi2D_, phi1D_, dphi1D_ )

end

function compOrder( orderRq )

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
