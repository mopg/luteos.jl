function hdgSolveElas( master::Master2D, mesh::Mesh2D, material::Material, problem::Problem)

dim     = 2 # Should grab this from Master
nelem   = size( mesh.nodes, 3 ) # Mumber of elements in mesh
nnodes  = size( mesh.nodes, 1 ) # Number of nodes in one element
szF     = (mesh.porder+1)       # Number of nodes per face
unkUhat = szF * dim             # Total of uhat unknowns per face

# Compute the stiffness tensor
Cstiff = material.Cstiff[2] #compCstiff( material, dim )

# Stability parameter
tau = Cstiff

### Initialize quantities

# Equations of motion (Newton's second law)
A = fill(0.0, nnodes*dim    , nnodes*dim^2  , nelem) # (v_{i,j}, \sigma^h_{ij})_{T^h} - <v_i, \sigma^h_{ij} n_j>_{dT^h}
B = fill(0.0, nnodes*dim    , nnodes*dim    , nelem) # <v_i, func(u^h_k) >_{dT^h}
C = fill(0.0, nnodes*dim    , szF*3*dim     , nelem) # <v_i, func(\hat{u}^h_k) >_{dT^h}

F = fill(0.0, nnodes*dim    , 1            , nelem) # (v_{i}, b_{i})_{T^h}

# strain-displacement
D = fill(0.0, nnodes*dim^2  , nnodes*dim^2  , nelem) # (w_{ij}  , \epsilon^h_{ij})_{T^h}
H = fill(0.0, nnodes*dim^2  , nnodes*dim    , nelem) # <w_{ij}  , (\hat{u}^h_i*n_j+ \hat{u}^h_j*n_i)>_{dT^h}
J = fill(0.0, nnodes*dim^2  , szF*3*dim     , nelem) # (w_{ij,j}, u^h_i)_{T^h} + (w_{ij,i}, u^h_j*n_i)_{T^h}

# Hooke's law
K = fill(0.0, nnodes*dim^2  , nnodes*dim^2  , nelem) # (z_{ij}  , \sigma^h_{ij})_{T^h}
M = fill(0.0, nnodes*dim^2  , nnodes*dim^2  , nelem) # (z_{ij}  , C_{ijkl}\epsilon^h_{kl})_{T^h}

# Flux jumps
N = fill(0.0, 3*szF*dim     , nnodes*dim^2  , nelem) # <\mu_{i} , \sigma^h_{ij}*n_j>_{dT^h\dOm_D}
P = fill(0.0, 3*szF*dim     , nnodes*dim    , nelem) # <\mu_{i} , func(u^h_k)>_{dT^h\dOm_D}
Q = fill(0.0, 3*szF*dim     , szF*3*dim     , nelem) # <\mu_{i} , func(\hat{u}^h_k)>_{dT^h\dOm_D} or <\mu_{i} , \hat{u}^h_{i})_{dOm_D}

# Dirichlet BC
#   BC RHS
G = fill(0.0, 3*szF*dim     , 1            , nelem) # (\mu_{i} , t_i)_{dT^h\dOm_D}) or (\mu_{i} , \bar{u}_i)_{dOm_D}

# Matrix with unknowns for uhath
A_UHATH = fill(0.0, 3*szF*dim     , szF*3*dim     , nelem)
B_UHATH = fill(0.0, 3*szF*dim     , 1            , nelem)

Rfull = fill(0.0, size(mesh.f,1)*unkUhat,1);

rr = 1 # iterator for unknowns in sparsity matrix

dphic = fill(0.0, size(master.dphi))

for pp in 1:nelem # Loop over all elements

  pLoc = master.phi' * mesh.nodes[:,:,pp]

  dphic[:,:,1] = master.dphi[:,:,1] * diagm( mesh.xix[:,pp,1] ) + master.dphi[:,:,2] * diagm( mesh.xix[:,pp,3] )
  dphic[:,:,2] = master.dphi[:,:,1] * diagm( mesh.xix[:,pp,2] ) + master.dphi[:,:,2] * diagm( mesh.xix[:,pp,4] )

  jcwd = diagm(mesh.jcw[:,pp])

  # ------------------------- Volume integrals ------------------------------- #
  ## A (first part)
  # (v_{i,j}, \sigma^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    A[1+(ii-1)*nnodes:ii*nnodes,1+(ij-1)*nnodes:ij*nnodes,pp] = dphic[:,:,jj] * jcwd * master.phi'
  end

  ## D
  # (w_{ij}  , \epsilon^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    D[1+(ij-1)*nnodes:ij*nnodes,1+(ij-1)*nnodes:ij*nnodes,pp] = master.phi * jcwd * master.phi'
  end

  ## H
  # (w_{ij,j}, u^h_i)_{T^h} + (w_{ij,i}, u^h_j)_{T^h}
  for ii = 1:dim, jj = 1:dim
    ij = (ii-1)*dim + jj
    H[1+(ij-1)*nnodes:ij*nnodes,1+(ii-1)*nnodes:ii*nnodes,pp] += 0.5*dphic[:,:,jj] * jcwd * master.phi'
    H[1+(ij-1)*nnodes:ij*nnodes,1+(jj-1)*nnodes:jj*nnodes,pp] += 0.5*dphic[:,:,ii] * jcwd * master.phi'
  end

  ## K
  # (z_{ij}  , \sigma^h_{ij})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    K[1+(ij-1)*nnodes:ij*nnodes,1+(ij-1)*nnodes:ij*nnodes,pp] = master.phi * jcwd * master.phi'
  end

  ## M
  # (z_{ij}  , C_{ijkl}\epsilon^h_{kl})_{T^h}
  for ii in 1:dim, jj in 1:dim
    ij = (ii-1)*dim + jj
    for kk in 1:dim, ll in 1:dim
      kl = (kk-1)*dim + ll
      M[1+(ij-1)*nnodes:ij*nnodes,1+(kl-1)*nnodes:kl*nnodes,pp] -= Cstiff[ii,jj,kk,ll] * master.phi * jcwd * master.phi'
    end
  end

  ## F
  # (v_{i}, b_{i})_{T^h}
  src = problem.source(pLoc)
  for ii in 1:dim
      F[1+(ii-1)*nnodes:ii*nnodes,1,pp] = master.phi * jcwd * src[:,ii]
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

    println( size(mesh.nodes) )
    println(nod)
    println(pp)
    p1d  = master.phi1d'  * mesh.nodes[nod,:,pp]
    dp1d = master.dphi1d' * mesh.nodes[nod,:,pp]

    jac1d  = sqrt( dp1d[:,1].^2 + dp1d[:,2].^2 )
    jcw1d  = master.gwts1d .* jac1d
    jcw1dd = diagm( jcw1d )

    if rotdir
      nL =  [ dp1d[:,2] -dp1d[:,1] ] ./ [jac1d jac1d]
    else
      nL = -[ dp1d[:,2] -dp1d[:,1] ] ./ [jac1d jac1d]
    end

    tL = [ nL[:,2] -nL[:,1] ]

    faceInd = (1 + (qq-1) * (mesh.porder + 1) ):qq*(mesh.porder+1)

    ## A (second part)
    # -<v_i, \sigma^h_{ij} n_j>_{dT^h}
    for ii in 1:dim, jj in 1:dim
      ij = (ii-1)*dim + jj
      A[(ii-1)*nnodes + nod,(ij-1)*nnodes + nod,pp] -= master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,jj]) )' # This last part was flipped in MATLAB
    end

    ## B
    # -<v_i, tau_{ijkl} u^h_k * n_l * n_j >_{dT^h}
    for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim #l
      B[(ii-1)*nnodes + nod,(kk-1)*nnodes + nod,pp] +=
        tau[ii,jj,kk,ll] * master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,ll].*nL[:,jj]) )' # This last part was flipped in MATLAB
    end

    ## C
    # <v_i, tau_{ijkl} \hat{u}^h_k * n_l * n_j >_{dT^h}
    for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
      C[(ii-1)*nnodes + nod,(qq-1)*szF + (kk-1)*szF + faceInd,pp] -=
          tau[ii,jj,kk,ll] * master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,ll] .* nL[:,jj]) )' # This last part was flipped in MATLAB
    end

    ## J
    # <w_{ij}  , (\hat{u}^h_i*n_j+ \hat{u}^h_j*n_i)>_{dT^h}
    for ii in 1:dim, jj in 1:dim
      ij = (ii-1)*dim + jj
      J[(ij-1)*nnodes + nod,(qq-1)*szF + (ii-1)*szF + faceInd,pp] -= 0.5*master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,jj]) )'
      J[(ij-1)*nnodes + nod,(qq-1)*szF + (jj-1)*szF + faceInd,pp] -= 0.5*master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,ii]) )'
    end

    indF = abs( mesh.t2f[pp,qq] )

    if ( mesh.f[ indF, 4 ] < 0 )
      # At boundary

      bNo = - mesh.f[ indF, 4 ]

      if problem.bctype[bNo] == 1
        # Dirichlet boundary condition

        ## Q
        # (\mu_{i} , \hat{u}^h_{i})_{dOm_D}
        for ii in 1:dim     #i
          Q[(ii-1)*szF + (qq-1)*szF + faceInd,(ii-1)*szF + (qq-1)*szF + faceInd,pp] = master.phi1d * jcw1dd * master.phi1d'
        end

        # BC RHS
        ## G
        # (\mu_{i} , \bar{u}_i)_{dOm_D}
        bcout = problem.bcfunc[bNo](p1d)
        for ii in 1:dim     #i
          G[(ii-1)*szF + (qq-1)*szF + faceInd, 1,pp] = master.phi1d * jcw1dd * bcout[:,ii]
        end

      elseif problem.bctype[bNo] == 2
        # Neumann boundary condition

        ## N
        # <\mu_{i} , \sigma^h_{ij}*n_j>_{dT^h\dOm_D}
        for ii in 1:dim, jj in 1:dim   #j
          ij = (ii-1)*dim + jj
          N[(ii-1)*szF + (qq-1)*szF + faceInd, (ij-1)*nnodes + nod,pp] +=
              master.phi1d * jcw1dd * (diag(nL[:,jj]) * master.phi1d)'
        end

        ## P
        # <\mu_{i} , -tau_{ijkl} u^h_k n_l n_j ) >_{dT^h\dOm_D}
        for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
          P[(ii-1)*szF + (qq-1)*szF + faceInd,(kk-1)*nnodes + nod,pp] -=
            tau[ii,jj,kk,ll] * master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,ll] .* nL[:,jj]) )'
        end

        ## Q
        # <\mu_{i} , tau_{ijkl} \hat{u}^h_k n_l n_j ) >_{dT^h\dOm_D}
        for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
          Q[(ii-1)*szF + (qq-1)*szF + faceInd,(kk-1)*szF + (qq-1)*szF + faceInd,pp] +=
            tau[ii,jj,kk,ll] * master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,ll] .* nL[:,jj]) )'
        end

        # BC RHS
        ## G
        # (\mu_{i} , t_i)_{dT^h\dOm_D})
        if problem.bcnorm
          # Boundary conditions are defined normal to the surface
          bcout = problem.bcfunc[bNo](p1d)
          for ii in 1:dim     #i
              temp = nL[:,ii] .* bcout[:,1] + tL[:,ii] .* bcout[:,2]
              G[(ii-1)*szF + (qq-1)*szF + faceInd, 1,pp] = master.phi1d * jcw1dd * temp
          end
        else
          # Boundary conditions are defined in the general coordinate system
          for ii in 1:dim     #i
              G[(ii-1)*szF + (qq-1)*szF + faceInd, 1,pp] = master.phi1d * jcw1dd * problem.bcfunc[bNo,ii](p1d)
          end
        end

      else
          error("BC type not recognized. 1 = Dirichlet, 2 = Neumann")

      end

    else
      # In interior

      ## N
      # <\mu_{i} , \sigma^h_{ij}*n_j>_{dT^h\dOm_D}
      for ii in 1:dim, jj in 1:dim
        ij = (ii-1)*dim + jj
        N[(ii-1)*szF + (qq-1)*szF + faceInd, (ij-1)*nnodes + nod,pp] += master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,jj]) )'
      end

      ## P
      # <\mu_{i} , -tau_{ijkl} u^h_k n_l n_j ) >_{dT^h\dOm_D}
      for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
        P[(ii-1)*szF + (qq-1)*szF + faceInd,(kk-1)*nnodes + nod,pp] -=
          tau[ii,jj,kk,ll] * master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,ll] .* nL[:,jj]) )'
      end

      ## Q
      # <\mu_{i} , tau_{ijkl} \hat{u}^h_k n_l n_j ) >_{dT^h\dOm_D}
      for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
        Q[(ii-1)*szF + (qq-1)*szF + faceInd,(kk-1)*szF + (qq-1)*szF + faceInd,pp] +=
          tau[ii,jj,kk,ll] * master.phi1d * jcw1dd * ( master.phi1d * diagm(nL[:,ll] .* nL[:,jj]) )'
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
  # for qq = 1:3
  #   indJ = [i for i=1:3]
  #   deleteat!(indJ,qq)
  #
  #   indFall = abs.(mesh.t2f[pp,:])
  #   deleteat!(indFall,qq)
  #
  #   indF = abs(mesh.t2f[pp,qq])
  #   Rfull[1+(indF-1)*unkUhat:indF*unkUhat] += B_UHATH[1+(qq-1)*unkUhat:qq*unkUhat,1,pp]
  #
  #   # TODO: Figure out how to do sparse matrices in Julia!
  #   indRow[rr:rr-1+unkUhat^2,1] = repmat(1+(indF-1)*unkUhat:indF*unkUhat,1,unkUhat)';
  #   indCol[rr:rr-1+unkUhat^2,1] = reshape(repmat(1+(indF-1)*unkUhat:indF*unkUhat,unkUhat,1),unkUhat^2,1);
  #   unkH(rr:rr-1+unkUhat^2,1)   = reshape(A_UHATH(1+(qq-1)*unkUhat:qq*unkUhat,1+(qq-1)*unkUhat:qq*unkUhat,pp),unkUhat^2,1);
  #
  #   rr = rr + unkUhat^2;
  #
  #   for tt = 1:2
  #     indRow(rr:rr-1+unkUhat^2,1) = repmat(1+(indF-1)*unkUhat:indF*unkUhat,1,unkUhat)';
  #     indCol(rr:rr-1+unkUhat^2,1) = reshape(repmat(1+(indFall(tt)-1)*unkUhat:indFall(tt)*unkUhat,unkUhat,1),unkUhat^2,1);
  #     unkH(rr:rr-1+unkUhat^2,1)   = reshape(A_UHATH(1+(qq-1)*unkUhat:qq*unkUhat,1+(indJ(tt)-1)*unkUhat:indJ(tt)*unkUhat,pp),unkUhat^2,1);
  #
  #     rr = rr + unkUhat^2;
  #   end
  # end
  # -------------------------------------------------------------------------- #
end # end element loop

### Compute approximate trace
# TODO: FIGURE OUT HOW TO DO SPARSE SOLVES IN JULIA

end # end function
