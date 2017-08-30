include("../integration/quadratureTriangle.jl")
include("../integration/quadratureLine.jl")
include("../integration/basisFuncLine.jl")
include("../integration/basisFuncTriangle.jl")


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

function Master2D( porder::Int64, pgauss::Int64 = 3*porder )

  (gpts_,   gwts_)   = quadratureTriangle( pgauss )
  (gpts1D_, gwts1D_) = quadratureLine( pgauss )

  (phi_,   dphi_)   = basisFuncTriangle( porder, gpts_[:,1], gpts_[:,2] )
  (phi1D_, dphi1D_) = basisFuncLine( porder, gpts1D_ )

  Master2D( porder, pgauss, gpts_, gwts_, gpts1D_, gwts1D_,
    phi_, dphi_, phi1D_, dphi1D_ )

end
