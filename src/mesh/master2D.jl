include("../integration/quadratureTriangle.jl")
include("../integration/quadratureLine.jl")
include("../integration/basisFuncLineLeg.jl")
include("../integration/basisFuncTriangleLeg.jl")
include("../integration/basisFuncLineLag.jl")
include("../integration/basisFuncTriangleLag.jl")

type Master2D

  porder::Int64         # Polynomial order of mesh
  pgauss::Int64         # Polynomial order to be integrated exactly

  gpts::Array{Float64}    # Gauss points  2D
  gwts::Array{Float64}    # Gauss weights 2D

  gpts1D::Array{Float64}  # Gauss points  1D
  gwts1D::Array{Float64}  # Gauss weights 1D

  phi::Array{Float64}     # Shape functions in 2D
  dphi::Array{Float64}    # Derivatives of shape functions in 2D

  phi1D::Array{Float64}   # Shape functions in 1D
  dphi1D::Array{Float64}  # Derivatives of shape functions in 1D

end

function Master2D( porder::Int64; pgauss::Int64 = 3*porder, typeb = "lag" )

  (go1D, go2D) = compOrder( pgauss )

  (gpts_,   gwts_)   = quadratureTriangle( Val{go2D} )
  (gpts1D_, gwts1D_) = quadratureLine( Val{go1D} )

  if typeb == "lag"
    (phi_,   dphi_)   = basisFuncTriangleLag( Val{porder}, gpts_[:,1], gpts_[:,2] )
    (phi1D_, dphi1D_) = basisFuncLineLag( Val{porder}, gpts1D_ )
  elseif typeb == "leg"
    (phi_,   dphi_)   = basisFuncTriangleLeg( Val{porder}, gpts_[:,1], gpts_[:,2] )
    (phi1D_, dphi1D_) = basisFuncLineLeg( Val{porder}, gpts1D_ )
  else
    error("Unknown type of basis function")
  end

  Master2D( porder, pgauss, gpts_, gwts_, gpts1D_, gwts1D_,
    phi_, dphi_, phi1D_, dphi1D_ )

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

  return order1D, order2D

end
