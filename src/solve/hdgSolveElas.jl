function hdgSolveElas( master::Master2D, mesh::Mesh2D, material::Material, problem::Problem)

dim     = 2 # Should grab this from Master
nelem   = size( mesh.nodes, 3 ) # Mumber of elements in mesh
nnodes  = size( mesh.nodes, 1 ) # Number of nodes in one element
szF     = (mesh.porder+1)       # Number of nodes per face
unkUhat = szF * dim             # Total of uhat unknowns per face

# Compute the stiffness tensor
Cstiff = material.Cstiff[2] #compCstiff( material, dim )

# Stability parameter
τ = Cstiff

### Initialize quantities

# Equations of motion (Newton's second law)
A = fill(0.0, nnodes*dim    , nnodes*dim^2  , nelem) # (v_{i,j}, σ^h_{ij})_{T^h} - <v_i, σ^h_{ij} n_j>_{∂T^h}
B = fill(0.0, nnodes*dim    , nnodes*dim    , nelem) # <v_i, func(u^h_k) >_{∂T^h}
C = fill(0.0, nnodes*dim    , szF*3*dim     , nelem) # <v_i, func(\hat{u}^h_k) >_{∂T^h}

F = fill(0.0, nnodes*dim    , 1            , nelem) # (v_{i}, b_{i})_{T^h}

# strain-displacement
D = fill(0.0, nnodes*dim^2  , nnodes*dim^2  , nelem) # (w_{ij}  , ϵ^h_{ij})_{T^h}
H = fill(0.0, nnodes*dim^2  , nnodes*dim    , nelem) # <w_{ij}  , (\hat{u}^h_i*n_j+ \hat{u}^h_j*n_i)>_{∂T^h}
J = fill(0.0, nnodes*dim^2  , szF*3*dim     , nelem) # (w_{ij,j}, u^h_i)_{T^h} + (w_{ij,i}, u^h_j*n_i)_{T^h}

# Hooke's law
K = fill(0.0, nnodes*dim^2  , nnodes*dim^2  , nelem) # (z_{ij}  , σ^h_{ij})_{T^h}
M = fill(0.0, nnodes*dim^2  , nnodes*dim^2  , nelem) # (z_{ij}  , C_{ijkl}ϵ^h_{kl})_{T^h}

# Flux jumps
N = fill(0.0, 3*szF*dim     , nnodes*dim^2  , nelem) # <μ_{i} , σ^h_{ij}*n_j>_{∂T^h\∂Ω_D}
P = fill(0.0, 3*szF*dim     , nnodes*dim    , nelem) # <μ_{i} , func(u^h_k)>_{∂T^h\∂Ω_D}
Q = fill(0.0, 3*szF*dim     , szF*3*dim     , nelem) # <μ_{i} , func(\hat{u}^h_k)>_{∂T^h\∂Ω_D} or <μ_{i} , \hat{u}^h_{i})_{∂Ω_D}

# Dirichlet BC
#   BC RHS
G = fill(0.0, 3*szF*dim     , 1            , nelem) # (μ_{i} , t_i)_{∂T^h\∂Ω_D}) or (μ_{i} , \bar{u}_i)_{∂Ω_D}

# Matrix with unknowns for uhath
A_UHATH = fill(0.0, 3*szF*dim     , szF*3*dim     , nelem)
B_UHATH = fill(0.0, 3*szF*dim     , 1            , nelem)

# Set up column and row indices for sparse matrix
indRow = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )
indCol = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )
indUnk = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )

Rfull = fill( 0.0, size(mesh.f,1)*unkUhat ) # RHS vector

rr = 1 # iterator for unknowns in sparsity matrix

∇ϕc = fill(0.0, size(master.∇ϕ))

for pp in 1:nelem # Loop over all elements

  pLoc = master.ϕ' * mesh.nodes[:,:,pp]

  ∇ϕc[:,:,1] = master.∇ϕ[:,:,1] * diagm( mesh.∂ξ∂x[:,pp,1] ) + master.∇ϕ[:,:,2] * diagm( mesh.∂ξ∂x[:,pp,3] )
  ∇ϕc[:,:,2] = master.∇ϕ[:,:,1] * diagm( mesh.∂ξ∂x[:,pp,2] ) + master.∇ϕ[:,:,2] * diagm( mesh.∂ξ∂x[:,pp,4] )

  jcwd = diagm(mesh.jcw[:,pp])

  # ------------------------- Volume integrals ------------------------------- #
  ## A (first part)
  # (v_{i,j}, σ^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    A[1+(ii-1)*nnodes:ii*nnodes,1+(ij-1)*nnodes:ij*nnodes,pp] = ∇ϕc[:,:,jj] * jcwd * master.ϕ'
  end

  ## D
  # (w_{ij}  , ϵ^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    D[1+(ij-1)*nnodes:ij*nnodes,1+(ij-1)*nnodes:ij*nnodes,pp] = master.ϕ * jcwd * master.ϕ'
  end

  ## H
  # (w_{ij,j}, u^h_i)_{T^h} + (w_{ij,i}, u^h_j)_{T^h}
  for ii = 1:dim, jj = 1:dim
    ij = (ii-1)*dim + jj
    H[1+(ij-1)*nnodes:ij*nnodes,1+(ii-1)*nnodes:ii*nnodes,pp] += 0.5*∇ϕc[:,:,jj] * jcwd * master.ϕ'
    H[1+(ij-1)*nnodes:ij*nnodes,1+(jj-1)*nnodes:jj*nnodes,pp] += 0.5*∇ϕc[:,:,ii] * jcwd * master.ϕ'
  end

  ## K
  # (z_{ij}  , σ^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    K[1+(ij-1)*nnodes:ij*nnodes,1+(ij-1)*nnodes:ij*nnodes,pp] = master.ϕ * jcwd * master.ϕ'
  end

  ## M
  # (z_{ij}  , C_{ijkl}ϵ^h_{kl})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    for kk in 1:dim, ll in 1:dim
      kl = (kk-1)*dim + ll
      M[1+(ij-1)*nnodes:ij*nnodes,1+(kl-1)*nnodes:kl*nnodes,pp] -= Cstiff[ii,jj,kk,ll] * master.ϕ * jcwd * master.ϕ'
    end
  end

  ## F
  # (v_{i}, b_{i})_{T^h}
  src = problem.source(pLoc)
  for ii in 1:dim
      F[1+(ii-1)*nnodes:ii*nnodes,1,pp] = master.ϕ * jcwd * src[:,ii]
  end

  # -------------------------------------------------------------------------- #

  # ------------------------- Boundary integrals ----------------------------- #
  for qq in 1:(dim+1) # loop over number of faces

    if mesh.t2f[pp,qq] < 0
      # face DOES NOT follow counter-clockwise rotation
      nod    = master.perm[:,qq,2]
      rotdir = true
    else
      # face follows counter-clockwise rotation
      nod    = master.perm[:,qq,1]
      rotdir = false
    end

    p1d  = master.ϕ1d'  * mesh.nodes[nod,:,pp]
    dp1d = master.∇ϕ1d' * mesh.nodes[nod,:,pp]

    jac1d  = sqrt( dp1d[:,1].^2 + dp1d[:,2].^2 )
    jcw1d  = master.gwts1d .* jac1d
    jcw1dd = diagm( jcw1d )

    if rotdir
      nL = -[ dp1d[:,2] -dp1d[:,1] ] ./ [jac1d jac1d]
    else
      nL =  [ dp1d[:,2] -dp1d[:,1] ] ./ [jac1d jac1d]
    end

    tL = [ nL[:,2] -nL[:,1] ]

    faceInd = dim * (qq-1) * (mesh.porder + 1) + ( 1:(mesh.porder+1) )
    # We need to premultiply with "dim" because those are the number of unknowns per face

    ## A (second part)
    # -<v_i, σ^h_{ij} n_j>_{∂T^h}
    for ii in 1:dim, jj in 1:dim
      ij = (ii-1)*dim + jj
      A[(ii-1)*nnodes + nod,(ij-1)*nnodes + nod,pp] -= master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,jj]) )'
    end

    ## B
    # -<v_i, τ_{ijkl} u^h_k * n_l * n_j >_{∂T^h}
    for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim #l
      B[(ii-1)*nnodes + nod,(kk-1)*nnodes + nod,pp] +=
        τ[ii,jj,kk,ll] * master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,ll].*nL[:,jj]) )'
    end

    ## C
    # <v_i, τ_{ijkl} \hat{u}^h_k * n_l * n_j >_{∂T^h}
    for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
      C[(ii-1)*nnodes + nod, (kk-1)*szF + faceInd,pp] -=
          τ[ii,jj,kk,ll] * master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,ll] .* nL[:,jj]) )'
    end

    ## J
    # <w_{ij}  , (\hat{u}^h_i*n_j + \hat{u}^h_j*n_i)>_{∂T^h}
    for ii in 1:dim, jj in 1:dim
      ij = (ii-1)*dim + jj
      J[(ij-1)*nnodes + nod, (ii-1)*szF + faceInd,pp] -= 0.5*master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,jj]) )'
      J[(ij-1)*nnodes + nod, (jj-1)*szF + faceInd,pp] -= 0.5*master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,ii]) )'
    end

    indF = abs( mesh.t2f[pp,qq] )

    if ( mesh.f[ indF, 4 ] < 0 )
      # At boundary

      bNo = - mesh.f[ indF, 4 ]

      if problem.bctype[bNo] == 1
        # Dirichlet boundary condition

        ## Q
        # (μ_{i} , \hat{u}^h_{i})_{∂Ω_D}
        for ii in 1:dim     #i
          Q[(ii-1)*szF + faceInd,(ii-1)*szF + faceInd,pp] = master.ϕ1d * jcw1dd * master.ϕ1d'
        end

        # BC RHS
        ## G
        # (μ_{i} , \bar{u}_i)_{∂Ω_D}
        bcout = problem.bcfunc[bNo](p1d)
        for ii in 1:dim     #i
          G[(ii-1)*szF + faceInd, 1,pp] = master.ϕ1d * jcw1dd * bcout[:,ii]
        end

      elseif problem.bctype[bNo] == 2
        # Neumann boundary condition

        ## N
        # <μ_{i} , σ^h_{ij}*n_j>_{∂T^h\∂Ω_D}
        for ii in 1:dim, jj in 1:dim   #j
          ij = (ii-1)*dim + jj
          N[(ii-1)*szF + faceInd, (ij-1)*nnodes + nod,pp] +=
              master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,jj]) )'
        end

        ## P
        # <μ_{i} , -τ_{ijkl} u^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
        for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
          P[(ii-1)*szF + faceInd,(kk-1)*nnodes + nod,pp] -=
            τ[ii,jj,kk,ll] * master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,ll] .* nL[:,jj]) )'
        end

        ## Q
        # <μ_{i} , τ_{ijkl} \hat{u}^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
        for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
          Q[(ii-1)*szF + faceInd,(kk-1)*szF + faceInd,pp] +=
            τ[ii,jj,kk,ll] * master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,ll] .* nL[:,jj]) )'
        end

        # BC RHS
        ## G
        # (μ_{i} , t_i)_{∂T^h\∂Ω_D})
        if problem.bcnorm
          # Boundary conditions are defined normal to the surface
          bcout = problem.bcfunc[bNo](p1d)
          for ii in 1:dim     #i
            temp = nL[:,ii] .* bcout[:,1] + tL[:,ii] .* bcout[:,2]
            G[(ii-1)*szF + faceInd, 1,pp] = master.ϕ1d * jcw1dd * temp
          end
        else
          # Boundary conditions are defined in the general coordinate system
          bcout = problem.bcfunc[bNo](p1d)
          for ii in 1:dim     #i
            G[(ii-1)*szF + faceInd, 1,pp] = master.ϕ1d * jcw1dd * bcout[:,ii]
          end
        end
      else
          error("BC type not recognized. 1 = Dirichlet, 2 = Neumann")
      end

    else
      # In interior

      ## N
      # <μ_{i} , σ^h_{ij}*n_j>_{∂T^h\∂Ω_D}
      for ii in 1:dim, jj in 1:dim
        ij = (ii-1)*dim + jj
        N[(ii-1)*szF + faceInd, (ij-1)*nnodes + nod,pp] +=
          master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,jj]) )'
      end

      ## P
      # <μ_{i} , -τ_{ijkl} u^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
      for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
        P[(ii-1)*szF + faceInd,(kk-1)*nnodes + nod,pp] -=
          τ[ii,jj,kk,ll] * master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,ll] .* nL[:,jj]) )'
      end

      ## Q
      # <μ_{i} , τ_{ijkl} \hat{u}^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
      for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
        Q[(ii-1)*szF + faceInd,(kk-1)*szF + faceInd,pp] +=
          τ[ii,jj,kk,ll] * master.ϕ1d * jcw1dd * ( master.ϕ1d * diagm(nL[:,ll] .* nL[:,jj]) )'
      end

    end # boundary if-statement

  end # end loop over element faces
  # -------------------------------------------------------------------------- #

  # ------------------------- Form elemental quantities ---------------------- #
  invKM = K[:,:,pp] \ M[:,:,pp]
  invDH = D[:,:,pp] \ H[:,:,pp]
  invDJ = D[:,:,pp] \ J[:,:,pp]

  SYS_U_UHATH_11 = B[:,:,pp] + A[:,:,pp] * invKM * invDH
  SYS_U_UHATH_12 = C[:,:,pp] + A[:,:,pp] * invKM * invDJ
  SYS_U_UHATH_21 = P[:,:,pp] + N[:,:,pp] * invKM * invDH
  SYS_U_UHATH_22 = Q[:,:,pp] + N[:,:,pp] * invKM * invDJ

  tempMat = SYS_U_UHATH_11 \ ( [ -SYS_U_UHATH_12 F[:,:,pp] ] )
  A_UHATH[:,:,pp] =  SYS_U_UHATH_21 * tempMat[:,1:end-1] + SYS_U_UHATH_22
  B_UHATH[:,:,pp] = -SYS_U_UHATH_21 * tempMat[:,end]     + G[:,:,pp]
  # -------------------------------------------------------------------------- #

  # ------------------------- Fill up complete H and R matrices -------------- #
  for qq = 1:(dim+1)
    indJ = [i for i=1:(dim+1)]
    deleteat!(indJ,qq)

    indFall = abs.(mesh.t2f[pp,:])
    deleteat!(indFall,qq)

    indF = abs(mesh.t2f[pp,qq])
    Rfull[1+(indF-1)*unkUhat:indF*unkUhat] += B_UHATH[1+(qq-1)*unkUhat:qq*unkUhat,1,pp]

    indRow[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, unkUhat )
    indCol[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, 1, unkUhat )'[:]
    indUnk[rr:rr-1+unkUhat^2] = A_UHATH[ 1+(qq-1)*unkUhat:qq*unkUhat, 1+(qq-1)*unkUhat:qq*unkUhat, pp ][:]

    rr = rr + unkUhat^2

    for tt = 1:2
      indRow[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, unkUhat )
      indCol[rr:rr-1+unkUhat^2] = repmat(  1+(indFall[tt]-1)*unkUhat:indFall[tt]*unkUhat, 1, unkUhat )'[:]
      indUnk[rr:rr-1+unkUhat^2] = A_UHATH[ 1+(qq-1)*unkUhat:qq*unkUhat, 1+(indJ[tt]-1)*unkUhat:indJ[tt]*unkUhat, pp ][:]

      rr = rr + unkUhat^2
    end
  end
  # -------------------------------------------------------------------------- #
end # end element loop

### Compute approximate trace
Hfull = sparse( indRow, indCol, indUnk, size(mesh.f,1) * unkUhat, size(mesh.f,1) * unkUhat )

uhath = Hfull \ Rfull
# ---------------------------------------------------------------------------- #

## Compute approximate scalar value and flux
uhathTri = fill( 0.0, unkUhat*(dim+1), 1,    nelem )
uh       = fill( 0.0, nnodes,         dim,   nelem )
σh       = fill( 0.0, nnodes,         dim^2, nelem )
ϵh       = fill( 0.0, nnodes,         dim^2, nelem )

for pp in 1:nelem
    # ----------- Find uhath corresponding to this element ------------------- #
    indF = abs.( mesh.t2f[pp,:] )
    uhathTri[ 1          :  unkUhat, 1, pp ] = uhath[ 1+(indF[1]-1)*unkUhat:indF[1]*unkUhat, 1]
    uhathTri[ 1+  unkUhat:2*unkUhat, 1, pp ] = uhath[ 1+(indF[2]-1)*unkUhat:indF[2]*unkUhat, 1]
    uhathTri[ 1+2*unkUhat:3*unkUhat, 1, pp ] = uhath[ 1+(indF[3]-1)*unkUhat:indF[3]*unkUhat, 1]
    # ------------------------------------------------------------------------ #
    # ----------- Compute approximate displacement value --------------------- #
    invKM = K[:,:,pp] \ M[:,:,pp]
    invDH = D[:,:,pp] \ H[:,:,pp]
    invDJ = D[:,:,pp] \ J[:,:,pp]

    SYS_U_UHATH_11 = B[:,:,pp] + A[:,:,pp] * invKM * invDH
    SYS_U_UHATH_12 = C[:,:,pp] + A[:,:,pp] * invKM * invDJ

    uTemp = SYS_U_UHATH_11 \ ( F[:,:,pp] - SYS_U_UHATH_12 * uhathTri[:,:,pp] )

    for ii in 1:dim
        uh[:,ii,pp] = uTemp[ 1 + (ii-1)*nnodes:nnodes*ii ]
    end
    # ------------------------------------------------------------------------ #
    # ----------- Compute approximate strain field --------------------------- #
    eTemp = - invDH*uTemp - invDJ * uhathTri[:,:,pp]
    for ii in 1:dim^2
        ϵh[:,ii,pp] = eTemp[ 1 + (ii-1)*nnodes:nnodes*ii ]
    end
    # ------------------------------------------------------------------------ #
    # ----------- Compute approximate stress field --------------------------- #
    sTemp = - invKM * eTemp
    for ii in 1:dim^2
        σh[ :,ii,pp ] = sTemp[ 1 + (ii-1)*nnodes:nnodes*ii ]
    end
    # ------------------------------------------------------------------------ #
end

return (uhath, uh, σh, ϵh )

end # end function
