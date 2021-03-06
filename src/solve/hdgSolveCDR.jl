# ---------------------------------------------------------------------------- #
#
#   hdgSolveElas.jl
#
#   Solve convection-diffusion equations (n-dimensional)
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    hdgSolve( master::Master, mesh::Mesh, problem::CDR )

Solves convection-diffusion equations for n-dimensional problems.
"""
function hdgSolve( master::Master, mesh::Mesh, problem::CDR; method=:LU)

dim     = mesh.dim
nelem   = size( mesh.nodes, 3 ) # Mumber of elements in mesh
nnodes  = size( mesh.nodes, 1 ) # Number of nodes in one element
nodfac  = master.nodfac # Number of nodes on a face

unkUhat = nodfac                # Total of uhat unknowns per face
nfaces  = dim + 1               # Number of faces per element

# Stability parameter
τ = 5

### Initialize quantities

# Equations of motion (Newton's second law)
A = fill(0.0, nnodes*dim    , nnodes*dim)       # (v_i, q^h_i)_{T^h}
B = fill(0.0, nnodes*dim    , nnodes    )       # <v_i,i,  >_{∂T^h}
C = fill(0.0, nnodes*dim    , nodfac*nfaces)    # <v_i, func(\hat{u}^h_k) >_{∂T^h}

N = fill(0.0, nnodes        , nnodes*dim)       # <v_i, func(\hat{u}^h_k) >_{∂T^h}
D = fill(0.0, nnodes        , nnodes)           # <v_i, func(\hat{u}^h_k) >_{∂T^h}
E = fill(0.0, nnodes        , nodfac*nfaces)    # <v_i, func(\hat{u}^h_k) >_{∂T^h}

K = fill(0.0, nodfac*nfaces , nnodes*dim)       # <v_i, func(\hat{u}^h_k) >_{∂T^h}
L = fill(0.0, nodfac*nfaces , nnodes)           # <v_i, func(\hat{u}^h_k) >_{∂T^h}
M = fill(0.0, nodfac*nfaces , nodfac*nfaces)    # <v_i, func(\hat{u}^h_k) >_{∂T^h}

G = fill(0.0, nodfac*nfaces , 1)                # (v_{i}, b_{i})_{T^h}

# Matrix with unknowns for uhath
H = fill(0.0, nodfac*nfaces , nodfac*nfaces)
R = fill(0.0, nodfac*nfaces , 1            )    # (v_{i}, b_{i})_{T^h}

zm = fill( 0.0, nnodes*dim, 1)

# Matrices that need to be saved in order to recover u and qh
F    = fill(0.0, nnodes        , 1,              nelem)                # (v_{i}, b_{i})_{T^h}
ABND = fill(0.0, nnodes*(dim+1), nnodes*(dim+1), nelem)
CE   = fill(0.0, nnodes*(dim+1), nodfac*nfaces,  nelem)

# Set up column and row indices for sparse matrix
indRow = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )
indCol = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )
indUnk = fill( 0.0, (dim+1)^2 * nelem * unkUhat^2 )

Rfull = fill( 0.0, size(mesh.f,1)*unkUhat ) # RHS vector

rr = 1 # iterator for unknowns in sparsity matrix

∇ϕc = fill(0.0, size(master.∇ϕ))

κ = problem.κ
c = problem.c

# preallocate
jcwd = fill( 0.0, size(master.∇ϕ,2), size(master.∇ϕ,2) )
∂ξ∂x = fill( 0.0, size(master.∇ϕ,2), dim^2 )
∂x∂ξ = fill( 0.0, size(master.∇ϕ,2), dim^2 )
ϕjϕn = [Array{Float64}(nodfac,nodfac) for idx in 1:dim]


for pp in 1:nelem # Loop over all elements

  # zero all matrices
  A .*= 0.0
  B .*= 0.0
  C .*= 0.0
  N .*= 0.0
  D .*= 0.0
  E .*= 0.0
  K .*= 0.0
  L .*= 0.0
  M .*= 0.0
  G .*= 0.0
  H .*= 0.0
  R .*= 0.0

  # Compute Jacobians
  compJacob!( master, mesh.nodes[:,:,pp], ∂ξ∂x, jcwd, ∂x∂ξ )

  pLoc = master.ϕ' * mesh.nodes[:,:,pp]

  ∇ϕc = getderbfel( master, ∂ξ∂x )

  # mass matrices
  ϕj  = master.ϕ * jcwd
  ϕjϕ = master.ϕ * jcwd * master.ϕ'

  # ------------------------- Volume integrals ------------------------------- #
  ## A
  # (v,qₕ)_{Tₕ}
  for ii in 1:dim
    A[1+(ii-1)*nnodes:ii*nnodes,1+(ii-1)*nnodes:ii*nnodes] = ϕjϕ
  end

  ## B
  # (w_{ij}  , ϵ^h_{ij})_{T^h}
  for ii in 1:dim
    B[1+(ii-1)*nnodes:ii*nnodes,:] = ∇ϕc[:,:,ii] * ϕj'
  end

  ## N
  N += -κ * B'

  ## D (second part)
  # (w_{ij,j}, u^h_i)_{T^h} + (w_{ij,i}, u^h_j)_{T^h}
  for ii in 1:dim
    D -= c[ii] * ∇ϕc[:,:,ii] * ϕj'
  end

  ## F
  # (z_{ij}  , σ^h_{ij})_{T^h}
  src       = problem.source(pLoc)
  F[:,:,pp] = ϕj * src

  # -------------------------------------------------------------------------- #

  # ------------------------- Boundary integrals ----------------------------- #
  for qq in 1:nfaces # loop over number of faces

    (ϕdm, pdm, nod, normal, jcwdm) = compJacobFace( mesh, master, pp, qq )

    jcwddm = diagm( jcwdm )

    faceInd = (qq-1) * nodfac + ( 1:nodfac )

    cn = fill(0.0, size(normal,1))
    for ii in 1:dim
      cn += c[ii] * normal[:,ii]
    end

    # mass matrices
    ϕjf  = ϕdm * jcwddm
    ϕjϕf = ϕjf * ϕdm'
    for ii in 1:dim
        ϕjϕn[ii] = ϕjf * ( ϕdm * diagm(normal[:,ii]) )'
    end

    ## C
    # -<v_i, σ^h_{ij} n_j>_{∂T^h}
    for ii in 1:dim
      C[(ii-1)*nnodes + nod,faceInd] -= ϕjϕn[ii]
    end

    ## D
    # <v_i, τ_{ijkl} u^h_k * n_l * n_j >_{∂T^h}
    D[nod,nod] += τ * ϕjϕf

    ## E
    # <v_i, τ_{ijkl} \hat{u}^h_k * n_l * n_j >_{∂T^h}
    E[nod,faceInd] = - τ * ϕjϕf + ϕdm * ( diagm(cn) * jcwddm ) * ϕdm'

    indF = abs( mesh.t2f[pp,qq] )

    if ( mesh.f[ indF, end ] < 0 )
      # At boundary

      bNo = - mesh.f[ indF, end ]

      if problem.bctype[bNo] == 1
        # Dirichlet boundary condition

        ## M
        # (μ_{i} , \hat{u}^h_{i})_{∂Ω_D}
        M[faceInd,faceInd] = ϕjϕf


        # BC RHS
        ## G
        # (μ_{i} , \bar{u}_i)_{∂Ω_D}
        bcout = problem.bcfunc[bNo](pdm)
        G[faceInd, 1] = ϕjf * bcout

      elseif problem.bctype[bNo] == 2
        # Neumann boundary condition

        ## K
        # <μ_{i} , q_{ij}*n_j>_{∂Ω_N}
        for ii in 1:dim
          K[faceInd, (ii-1)*nnodes + nod] = ϕjϕn[ii]
        end

        # BC RHS
        ## G
        # (μ_{i} , \bar{∂u∂x})_{∂Ω_N})
        bcout = problem.bcfunc[bNo](pdm)
        G[faceInd, 1] = ϕjf * bcout
      else
          error("hdgSolveElas:: BC type not recognized. 1 = Dirichlet, 2 = Neumann")
      end

    else
      # In interior

      ## K
      # <μ_{i} , σ^h_{ij}*n_j>_{∂T^h\∂Ω_D}
      for ii in 1:dim
        K[faceInd, (ii-1)*nnodes + nod] = - ϕjϕn[ii]
      end

      ## L
      # <μ_{i} , -τ_{ijkl} u^h_k n_l n_j ) >_{∂T^h\∂Ω_D}

      L[faceInd,nod] = τ * ϕjϕf

      ## M
      # <μ_{i} , τ_{ijkl} \hat{u}^h_k n_l n_j ) >_{∂T^h\∂Ω_D}
      M[faceInd,faceInd] = - τ * ϕjϕf + ϕdm * ( diagm(cn) * jcwddm ) * ϕdm'

    end # boundary if-statement

  end # end loop over element faces
  # -------------------------------------------------------------------------- #

  # ------------------------- Form elemental quantities ---------------------- #
  ABND[:,:,pp] = [ A B;
                   N D ]

  CE[:,:,pp] = [ C; E ]

  H += M - [K  L] * ( ABND[:,:,pp] \ CE[:,:,pp] )
  R += G - [K  L] * ( ABND[:,:,pp] \ [ zm; F[:,:,pp] ] )

  # -------------------------------------------------------------------------- #

  # ------------------------- Fill up complete H and R matrices -------------- #
  for qq = 1:nfaces
    indJ = [i for i=1:nfaces]
    deleteat!(indJ,qq)

    indFall = abs.(mesh.t2f[pp,1:nfaces])
    deleteat!(indFall,qq)

    indF = abs(mesh.t2f[pp,qq])
    Rfull[1+(indF-1)*unkUhat:indF*unkUhat] += R[1+(qq-1)*unkUhat:qq*unkUhat,1]

    indRow[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, unkUhat )
    indCol[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, 1, unkUhat )'[:]
    indUnk[rr:rr-1+unkUhat^2] = H[ 1+(qq-1)*unkUhat:qq*unkUhat, 1+(qq-1)*unkUhat:qq*unkUhat ][:]

    rr = rr + unkUhat^2

    for tt = 1:(nfaces-1)
      indRow[rr:rr-1+unkUhat^2] = repmat(  1+(indF-1)*unkUhat:indF*unkUhat, unkUhat )
      indCol[rr:rr-1+unkUhat^2] = repmat(  1+(indFall[tt]-1)*unkUhat:indFall[tt]*unkUhat, 1, unkUhat )'[:]
      indUnk[rr:rr-1+unkUhat^2] = H[ 1+(qq-1)*unkUhat:qq*unkUhat, 1+(indJ[tt]-1)*unkUhat:indJ[tt]*unkUhat ][:]

      rr = rr + unkUhat^2
    end
  end
  # -------------------------------------------------------------------------- #
end # end element loop

### Compute approximate trace
Hfull = sparse( indRow, indCol, indUnk, size(mesh.f,1) * unkUhat, size(mesh.f,1) * unkUhat )

#fact = ILU.crout_ilu( Hfull, τ = 0.1 )
uhath = [0.]
if method == :LU
    uhath = Hfull \ Rfull
elseif method == :GMRES
    uhath = IterativeSolvers.gmres( Hfull, Rfull, Pl=fact )
end

# ---------------------------------------------------------------------------- #

## Compute approximate scalar value and flux
uhathTri = fill( 0.0, unkUhat*nfaces, 1 )
uh       = fill( 0.0, nnodes,         1,     nelem )
qh       = fill( 0.0, nnodes,         dim,   nelem )

for pp in 1:nelem
    # ----------- Find uhath corresponding to this element ------------------- #
    indF = abs.( mesh.t2f[pp,1:nfaces] )
    for qq in 1:nfaces
      uhathTri[ 1+(qq-1)*unkUhat:qq*unkUhat, 1 ] = uhath[ 1+(indF[qq]-1)*unkUhat:indF[qq]*unkUhat ]
    end
    # ------------------------------------------------------------------------ #
    # ----------- Compute approximate displacement value --------------------- #

    rhsTemp = [ zm; F[:,:,pp] ] - CE[:,:,pp] * uhathTri

    uTemp   = ABND[:,:,pp] \ rhsTemp

    for ii in 1:dim
      qh[:,ii,pp] = uTemp[ 1+(ii-1)*nnodes:ii*nnodes ]
    end
    uh[:,1,pp] = uTemp[ 1+dim*nnodes:(dim+1)*nnodes ]
end

return ( uhath, uh, qh )

end # end hdgSolveCD
