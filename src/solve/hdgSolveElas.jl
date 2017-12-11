# ---------------------------------------------------------------------------- #
#
#   hdgSolveElas.jl
#
#   Solve linear elasticity equations (n-dimensional)
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    hdgSolveElas( master::Master, mesh::Mesh, material::Material,
                  problem::Problem)

Solves the linear elasticity equations for n-dimensional problems.
"""
function hdgSolveElas( master::Master, mesh::Mesh, material::Material, problem::Problem)

dim     = mesh.dim
nelem   = size( mesh.nodes, 3 ) # Mumber of elements in mesh
nnodes  = size( mesh.nodes, 1 ) # Number of nodes in one element
nodfac  = master.nodfac # Number of nodes on a face

unkUhat = nodfac * dim          # Total of uhat unknowns per face
nfaces  = dim + 1               # Number of faces per element

# Compute the stiffness tensor
Cstiff = material.Cstiff

# Stability parameter
τ = Cstiff

### Initialize quantities

# Equations of motion (Newton's second law)
A = fill(0.0, nnodes*dim    , nnodes*dim^2      )    # (v_{i,j}, σ^h_{ij})_{T^h} - <v_i, σ^h_{ij} n_j>_{∂T^h}
B = fill(0.0, nnodes*dim    , nnodes*dim        )    # <v_i, func(u^h_k) >_{∂T^h}
C = fill(0.0, nnodes*dim    , nodfac*nfaces*dim )    # <v_i, func(\hat{u}^h_k) >_{∂T^h}

F = fill(0.0, nnodes*dim    , 1                 )    # (v_{i}, b_{i})_{T^h}

# strain-displacement
D = fill(0.0, nnodes*dim^2  , nnodes*dim^2      )    # (w_{ij}  , ϵ^h_{ij})_{T^h}
H = fill(0.0, nnodes*dim^2  , nnodes*dim        )    # <w_{ij}  , (\hat{u}^h_i*n_j+ \hat{u}^h_j*n_i)>_{∂T^h}
J = fill(0.0, nnodes*dim^2  , nodfac*nfaces*dim )    # (w_{ij,j}, u^h_i)_{T^h} + (w_{ij,i}, u^h_j*n_i)_{T^h}

# Hooke's law
K = fill(0.0, nnodes*dim^2  , nnodes*dim^2      )    # (z_{ij}  , σ^h_{ij})_{T^h}
M = fill(0.0, nnodes*dim^2  , nnodes*dim^2      )    # (z_{ij}  , C_{ijkl}ϵ^h_{kl})_{T^h}

# Flux jumps
N = fill(0.0, nfaces*nodfac*dim, nnodes*dim^2      ) # <μ_{i} , σ^h_{ij}*n_j>_{∂T^h\∂Ω_D}
P = fill(0.0, nfaces*nodfac*dim, nnodes*dim        ) # <μ_{i} , func(u^h_k)>_{∂T^h\∂Ω_D}
Q = fill(0.0, nfaces*nodfac*dim, nodfac*nfaces*dim ) # <μ_{i} , func(\hat{u}^h_k)>_{∂T^h\∂Ω_D} or <μ_{i} , \hat{u}^h_{i})_{∂Ω_D}

# Dirichlet BC
#   BC RHS
G = fill(0.0, nfaces*nodfac*dim, 1                ) # (μ_{i} , t_i)_{∂T^h\∂Ω_D}) or (μ_{i} , \bar{u}_i)_{∂Ω_D}

# Matrix with unknowns for uhath
A_UHATH = fill(0.0, nfaces*nodfac*dim, nodfac*nfaces*dim)
B_UHATH = fill(0.0, nfaces*nodfac*dim, 1                )

# Matrices that need to be saved in order to recover u, σ and ϵ at the end
sysUUH_11_12 = fill(0.0, nnodes*dim    , nodfac*nfaces*dim, nelem )
sysUUH_11_F  = fill(0.0, nnodes*dim    , 1,                 nelem )
KMDH         = fill(0.0, nnodes*dim^2  , nnodes*dim,        nelem )
KMDJ         = fill(0.0, nnodes*dim^2  , nodfac*nfaces*dim, nelem )

# Set up column and row indices for sparse matrix
indRow = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )
indCol = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )
indUnk = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )

Rfull = fill( 0.0, size(mesh.f,1)*unkUhat ) # RHS vector

rr = 1 # iterator for unknowns in sparsity matrix

# preallocate
jcwd = fill( 0.0, size(master.∇ϕ,2), size(master.∇ϕ,2) )
∂ξ∂x = fill( 0.0, size(master.∇ϕ,2), dim^2 )
∂x∂ξ = fill( 0.0, size(master.∇ϕ,2), dim^2 )

@time for pp in 1:nelem # Loop over all elements

  # zero out all matrices
  A .*= 0.0
  B .*= 0.0
  C .*= 0.0
  F .*= 0.0
  D .*= 0.0
  H .*= 0.0
  J .*= 0.0
  K .*= 0.0
  M .*= 0.0
  N .*= 0.0
  P .*= 0.0
  Q .*= 0.0
  G .*= 0.0
  A_UHATH .*= 0.0
  B_UHATH .*= 0.0

  # Compute Jacobians
  compJacob!( master, mesh.nodes[:,:,pp], ∂ξ∂x, jcwd, ∂x∂ξ )

  pLoc = master.ϕ' * mesh.nodes[:,:,pp]

  ∇ϕc = getderbfel( master, ∂ξ∂x )

  # ------------------------- Volume integrals ------------------------------- #
  ## A (first part)
  # (v_{i,j}, σ^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    A[1+(ii-1)*nnodes:ii*nnodes,1+(ij-1)*nnodes:ij*nnodes] = ∇ϕc[:,:,jj] * jcwd * master.ϕ'
  end

  ## D
  # (w_{ij}  , ϵ^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    D[1+(ij-1)*nnodes:ij*nnodes,1+(ij-1)*nnodes:ij*nnodes] = master.ϕ * jcwd * master.ϕ'
  end

  ## H
  # (w_{ij,j}, u^h_i)_{T^h} + (w_{ij,i}, u^h_j)_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    H[1+(ij-1)*nnodes:ij*nnodes,1+(ii-1)*nnodes:ii*nnodes] += 0.5 * ∇ϕc[:,:,jj] * jcwd * master.ϕ'
    H[1+(ij-1)*nnodes:ij*nnodes,1+(jj-1)*nnodes:jj*nnodes] += 0.5 * ∇ϕc[:,:,ii] * jcwd * master.ϕ'
  end

  ## K
  # (z_{ij}  , σ^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    K[1+(ij-1)*nnodes:ij*nnodes,1+(ij-1)*nnodes:ij*nnodes] = master.ϕ * jcwd * master.ϕ'
  end

  ## M
  # (z_{ij}  , C_{ijkl}ϵ^h_{kl})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    for kk in 1:dim, ll in 1:dim
      kl = (kk-1)*dim + ll
      M[1+(ij-1)*nnodes:ij*nnodes,1+(kl-1)*nnodes:kl*nnodes] -= Cstiff[ii,jj,kk,ll] * master.ϕ * jcwd * master.ϕ'
    end
  end

  ## F
  # (v_{i}, b_{i})_{T^h}
  src = problem.source(pLoc)
  for ii in 1:dim
      F[1+(ii-1)*nnodes:ii*nnodes,1] = master.ϕ * jcwd * src[:,ii]
  end

  # -------------------------------------------------------------------------- #

  # ------------------------- Boundary integrals ----------------------------- #
  for qq in 1:nfaces # loop over number of faces

    (ϕdm, pdm, nod, normal, jcwdm) = compJacobFace( mesh, master, pp, qq )

    jcwddm = diagm( jcwdm )

    faceInd = dim * (qq-1) * nodfac + ( 1:nodfac )
    # We need to premultiply with "dim" because those are the number of unknowns per face

    ## A (second part)
    # -<v_i, σ^h_{ij} n_j>_{∂T^h}
    for ii in 1:dim, jj in 1:dim
      ij = (ii-1)*dim + jj
      A[(ii-1)*nnodes + nod,(ij-1)*nnodes + nod] -= ϕdm * jcwddm * ( ϕdm * diagm(normal[:,jj]) )'
    end

    ## B
    # <v_i, τ_{ijkl} u^h_k * n_l * n_j >_{∂T^h}
    for ll in 1:dim, kk in 1:dim, jj in 1:dim, ii in 1:dim # column-major ordering
      B[(ii-1)*nnodes + nod,(kk-1)*nnodes + nod] +=
        τ[ii,jj,kk,ll] * ϕdm * jcwddm * ( ϕdm * diagm(normal[:,ll].*normal[:,jj]) )'
    end

    ## C
    # <v_i, τ_{ijkl} \hat{u}^h_k * n_l * n_j >_{∂T^h}
    for ll in 1:dim, kk in 1:dim, jj in 1:dim, ii in 1:dim # column-major ordering
      C[(ii-1)*nnodes + nod, (kk-1)*nodfac + faceInd] -=
          τ[ii,jj,kk,ll] * ϕdm * jcwddm * ( ϕdm * diagm(normal[:,ll] .* normal[:,jj]) )'
    end

    ## J
    # <w_{ij}  , (\hat{u}^h_i*n_j + \hat{u}^h_j*n_i)>_{∂T^h}
    for ii in 1:dim, jj in 1:dim
      ij = (ii-1)*dim + jj
      J[(ij-1)*nnodes + nod, (ii-1)*nodfac + faceInd] -= 0.5 * ϕdm * jcwddm * ( ϕdm * diagm(normal[:,jj]) )'
      J[(ij-1)*nnodes + nod, (jj-1)*nodfac + faceInd] -= 0.5 * ϕdm * jcwddm * ( ϕdm * diagm(normal[:,ii]) )'
    end

    indF = abs( mesh.t2f[pp,qq] )

    if ( mesh.f[ indF, end ] < 0 )
      # At boundary

      bNo = - mesh.f[ indF, end ]

      if problem.bctype[bNo] == 1
        # Dirichlet boundary condition

        ## Q
        # (μ_{i} , \hat{u}^h_{i})_{∂Ω_D}
        for ii in 1:dim     #i
          Q[(ii-1)*nodfac + faceInd,(ii-1)*nodfac + faceInd] = ϕdm * jcwddm * ϕdm'
        end

        # BC RHS
        ## G
        # (μ_{i} , \bar{u}_i)_{∂Ω_D}
        bcout = problem.bcfunc[bNo](pdm)
        for ii in 1:dim     #i
          G[(ii-1)*nodfac + faceInd, 1] = ϕdm * jcwddm * bcout[:,ii]
        end

      elseif problem.bctype[bNo] == 2
        # Neumann boundary condition

        ## N
        # <μ_{i} , σ^h_{ij}*n_j>_{∂T^h\∂Ω_D}
        for ii in 1:dim, jj in 1:dim   #j
          ij = (ii-1)*dim + jj
          N[(ii-1)*nodfac + faceInd, (ij-1)*nnodes + nod] +=
              ϕdm * jcwddm * ( ϕdm * diagm(normal[:,jj]) )'
        end

        ## P
        # <μ_{i} , -τ_{ijkl} u^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
        for ll in 1:dim, kk in 1:dim, jj in 1:dim, ii in 1:dim # column-major ordering
          P[(ii-1)*nodfac + faceInd,(kk-1)*nnodes + nod] -=
            τ[ii,jj,kk,ll] * ϕdm * jcwddm * ( ϕdm * diagm(normal[:,ll] .* normal[:,jj]) )'
        end

        ## Q
        # <μ_{i} , τ_{ijkl} \hat{u}^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
        for ll in 1:dim, kk in 1:dim, jj in 1:dim, ii in 1:dim # column-major ordering
          Q[(ii-1)*nodfac + faceInd,(kk-1)*nodfac + faceInd] +=
            τ[ii,jj,kk,ll] * ϕdm * jcwddm * ( ϕdm * diagm(normal[:,ll] .* normal[:,jj]) )'
        end

        # BC RHS
        ## G
        # (μ_{i} , t_i)_{∂T^h\∂Ω_D})
        if problem.bcnorm
          # Boundary conditions are defined normal to the surface
          bcout = problem.bcfunc[bNo](pdm)
          for ii in 1:dim     #i
            temp = normal[:,ii] .* bcout[:,1]# + tL[:,ii] .* bcout[:,2] # NOTE: No tangential component, because hard to define in 3D
            G[(ii-1)*nodfac + faceInd, 1] = ϕdm * jcwddm * temp
          end
        else
          # Boundary conditions are defined in the general coordinate system
          bcout = problem.bcfunc[bNo](pdm)
          for ii in 1:dim     #i
            G[(ii-1)*nodfac + faceInd, 1] = ϕdm * jcwddm * bcout[:,ii]
          end
        end
      else
          error("hdgSolveElas:: BC type not recognized. 1 = Dirichlet, 2 = Neumann")
      end

    else
      # In interior

      ## N
      # <μ_{i} , σ^h_{ij}*n_j>_{∂T^h\∂Ω_D}
      for ii in 1:dim, jj in 1:dim
        ij = (ii-1)*dim + jj
        N[(ii-1)*nodfac + faceInd, (ij-1)*nnodes + nod] +=
          ϕdm * jcwddm * ( ϕdm * diagm(normal[:,jj]) )'
      end

      ## P
      # <μ_{i} , -τ_{ijkl} u^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
      for ll in 1:dim, kk in 1:dim, jj in 1:dim, ii in 1:dim # column-major ordering
        P[(ii-1)*nodfac + faceInd,(kk-1)*nnodes + nod] -=
          τ[ii,jj,kk,ll] * ϕdm * jcwddm * ( ϕdm * diagm(normal[:,ll] .* normal[:,jj]) )'
      end

      ## Q
      # <μ_{i} , τ_{ijkl} \hat{u}^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
      for ll in 1:dim, kk in 1:dim, jj in 1:dim, ii in 1:dim # column-major ordering
        Q[(ii-1)*nodfac + faceInd,(kk-1)*nodfac + faceInd] +=
          τ[ii,jj,kk,ll] * ϕdm * jcwddm * ( ϕdm * diagm(normal[:,ll] .* normal[:,jj]) )'
      end

    end # boundary if-statement

  end # end loop over element faces
  # -------------------------------------------------------------------------- #

  # ------------------------- Form elemental quantities ---------------------- #
  invKM = K \ M
  invDH = D \ H
  invDJ = D \ J

  SYS_U_UHATH_11 = B + A * invKM * invDH
  SYS_U_UHATH_12 = C + A * invKM * invDJ
  SYS_U_UHATH_21 = P + N * invKM * invDH
  SYS_U_UHATH_22 = Q + N * invKM * invDJ

  tempMat =  SYS_U_UHATH_11 \ ( [ -SYS_U_UHATH_12 F ] )
  A_UHATH =  SYS_U_UHATH_21 * tempMat[:,1:end-1] + SYS_U_UHATH_22
  B_UHATH = -SYS_U_UHATH_21 * tempMat[:,end]     + G

  sysUUH_11_12[:,:,pp] = tempMat[:,1:end-1]
  sysUUH_11_F[:,:,pp]  = tempMat[:,end]
  KMDH[:,:,pp]         = invKM * invDH
  KMDJ[:,:,pp]         = invKM * invDJ
  # -------------------------------------------------------------------------- #

  # ------------------------- Fill up complete H and R matrices -------------- #
  for qq = 1:nfaces
    indJ = [i for i=1:nfaces]
    deleteat!(indJ,qq)

    indFall = abs.(mesh.t2f[pp,1:nfaces])
    deleteat!(indFall,qq)

    indF = abs(mesh.t2f[pp,qq])
    Rfull[1+(indF-1)*unkUhat:indF*unkUhat] += B_UHATH[1+(qq-1)*unkUhat:qq*unkUhat,1]

    indRow[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, unkUhat )
    indCol[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, 1, unkUhat )'[:]
    indUnk[rr:rr-1+unkUhat^2] = A_UHATH[ 1+(qq-1)*unkUhat:qq*unkUhat, 1+(qq-1)*unkUhat:qq*unkUhat ][:]

    rr = rr + unkUhat^2

    for tt = 1:(nfaces-1)
      indRow[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, unkUhat )
      indCol[rr:rr-1+unkUhat^2] = repmat(  1+(indFall[tt]-1)*unkUhat:indFall[tt]*unkUhat, 1, unkUhat )'[:]
      indUnk[rr:rr-1+unkUhat^2] = A_UHATH[ 1+(qq-1)*unkUhat:qq*unkUhat, 1+(indJ[tt]-1)*unkUhat:indJ[tt]*unkUhat ][:]

      rr = rr + unkUhat^2
    end
  end
  # -------------------------------------------------------------------------- #
end # end element loop

### Compute approximate trace
Hfull = sparse( indRow, indCol, indUnk, size(mesh.f,1) * unkUhat, size(mesh.f,1) * unkUhat )

# @time fact = ILU.crout_ilu( Hfull, τ = 0.1 )
@time uhath = Hfull \ Rfull #IterativeSolvers.gmres( Hfull, Rfull, Pl=fact )
# ---------------------------------------------------------------------------- #

## Compute approximate scalar value and flux
uhathTri = fill( 0.0, unkUhat*nfaces, 1 )
uh       = fill( 0.0, nnodes,         dim,   nelem )
σh       = fill( 0.0, nnodes,         dim^2, nelem )

for pp in 1:nelem
    # ----------- Find uhath corresponding to this element ------------------- #
    indF = abs.( mesh.t2f[pp,1:nfaces] )
    for qq in 1:nfaces
      uhathTri[ 1+(qq-1)*unkUhat:qq*unkUhat, 1 ] = uhath[ 1+(indF[qq]-1)*unkUhat:indF[qq]*unkUhat ]
    end
    # ------------------------------------------------------------------------ #
    # ----------- Compute approximate displacement value --------------------- #
    uTemp = sysUUH_11_F[:,:,pp] + sysUUH_11_12[:,:,pp] * uhathTri
    for ii in 1:dim
        uh[:,ii,pp] = uTemp[ 1 + (ii-1)*nnodes:nnodes*ii ]
    end
    # ------------------------------------------------------------------------ #
    # ----------- Compute approximate stress field --------------------------- #
    sTemp = KMDH[:,:,pp] * uTemp + KMDJ[:,:,pp] * uhathTri
    for ii in 1:dim^2
        σh[ :,ii,pp ] = sTemp[ 1 + (ii-1)*nnodes:nnodes*ii ]
    end
    # ------------------------------------------------------------------------ #
end

return (uhath, uh, σh)

end # end hdgSolveElas
