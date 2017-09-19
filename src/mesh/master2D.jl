# include("../integration/quadratureTriangle.jl")
# include("../integration/quadratureLine.jl")
# include("../integration/basisFuncLineLeg.jl")
# include("../integration/basisFuncTriangleLeg.jl")
# include("../integration/basisFuncLineLag.jl")
# include("../integration/basisFuncTriangleLag.jl")

type Master2D

  dim::Int64              # Dimension of the problem

  porder::Int64           # Polynomial order of mesh
  pgauss::Int64           # Polynomial order to be integrated exactly

  gpts::Array{Float64}    # Gauss points  2D
  gwts::Array{Float64}    # Gauss weights 2D

  gpts1d::Array{Float64}  # Gauss points  1D
  gwts1d::Array{Float64}  # Gauss weights 1D

  ϕ::Array{Float64}       # Shape functions in 2D
  ∇ϕ::Array{Float64}      # Derivatives of shape functions in 2D

  ϕ1d::Array{Float64}     # Shape functions in 1D
  ∇ϕ1d::Array{Float64}    # Derivatives of shape functions in 1D

  perm::Array{Int64}      # Node numbers on faces

end

function Master2D( porder::Int64; pgauss::Int64 = 3*porder, typeb = "lag" )

  (go1D, go2D) = compOrder( pgauss )

  (gpts_,   gwts_)   = quadratureTriangle( Val{go2D} )
  (gpts1d_, gwts1d_) = quadratureLine( Val{go1D} )

  if typeb == "lag"
    (ϕ_,   ∇ϕ_)   = basisFuncTriangleLag( Val{porder}, gpts_[:,1], gpts_[:,2] )
    (ϕ1d_, ∇ϕ1d_) = basisFuncLineLag( Val{porder}, gpts1d_ )
  elseif typeb == "leg"
    (ϕ_,   ∇ϕ_)   = basisFuncTriangleLeg( Val{porder}, gpts_[:,1], gpts_[:,2] )
    (ϕ1d_, ∇ϕ1d_) = basisFuncLineLeg( Val{porder}, gpts1d_ )
  else
    error("Unknown type of basis function")
  end

  perm_ = findPerm( porder )

  Master2D( 2, porder, pgauss, gpts_, gwts_, gpts1d_, gwts1d_,
    ϕ_, ∇ϕ_, ϕ1d_, ∇ϕ1d_, perm_ )

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

function findPerm( p )

  permt = fill( 0::Int64, p+1, 3 )

  szF = p + 1

  for qq in 1:3

    # corners
    c1 = qq + 1
    c2 = qq + 2
    c1 = c1 - (c1 > 3)*3
    c2 = c2 - (c2 > 3)*3

    permt[1,qq] = c1
    permt[2,qq] = c2

    # faces
    permt[3:p+1,qq] = 3 + ( (1+(qq-1)*(szF-2)):(qq*(szF-2)) )

  end

  perm = fill( 0::Int64, p+1, 3, 2 )
  perm[:,:,1]   = permt
  perm[1:2,:,2] = permt[[2,1],:]
  perm[3:(p+1),:,2] = permt[(p+1):-1:3,:]

  return perm

end
