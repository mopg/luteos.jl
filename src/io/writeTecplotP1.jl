# ---------------------------------------------------------------------------- #
#
#   writeTecplot.jl
#
#   Several functions to write Tecplot output
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    writeTecplot( flname::String, prob::Elas, mesh::Mesh3D,
                  uh::Array{Float64,3}, σh::Array{Float64,3};
                  σmises = Array{Float64,3}(0,0,0) )

Outputs tecplot data file to `flname` for 3D structural elasticity problems.
"""
function writeTecplotP1( flname::String, prob::luteos.Elas, mesh::Mesh3D,
                         uh::Matrix{Float64}, σh::Matrix{Float64}; σmises = Matrix{Float64}(0,0) )


    nelem = size( mesh.nodes, 3 )
    ndat  = 2*mesh.dim + mesh.dim^2
    if size(σmises,1) > 0
        ndat += 1
    end

    # Open file
    fid = open( flname, "w" )

    ### Write header
    @printf( fid, "TITLE = \"%s\"\n", prob.name )
    @printf( fid, "VARIABLES = \"X\", \"Y\", \"Z\"" )
    for jj in 1:mesh.dim
        @printf( fid, ", \"U<sub>%i</sub>\"",jj )
    end
    for ii in 1:mesh.dim, jj in 1:mesh.dim
        @printf( fid, ", \"<greek>s</greek><sub>%i%i</sub>\"", ii, jj )
    end
    if size( σmises, 1 ) > 0
        @printf( fid, ", \"<greek>s</greek><sub>mises</sub>\"" )
    end
    @printf( fid, "\n" )
    @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETETRAHEDRON\n", mesh.n, nelem )

    ### Write values
    for pp in 1:mesh.n
        if size( σmises, 1 ) == 0
            wrtval = [mesh.p[pp,:]; uh[pp,:]; σh[pp,:] ]
        else
            wrtval = [mesh.p[pp,:]; uh[pp,:]; σh[pp,:]; σmises[pp,:] ]
        end

    for ii in 1:ndat
      @printf(fid, "%16.15e\t", wrtval[ii] )
    end
    @printf(fid,"\n")

    end

    ### Write connectivity
    for ee in 1:nelem

        for kk in 1:(mesh.dim+1)
            @printf(fid, "%i\t", mesh.t[ee,kk] )
        end
        @printf(fid,"\n")

    end

    close(fid)

end

"""
    writeTecplot( flname::String, prob::Elas, mesh::Mesh3D,
                  uh::Array{Float64,3}, σh::Array{Float64,3};
                  σmises = Array{Float64,3}(0,0,0) )

Outputs tecplot data file to `flname` for 3D structural elasticity problems.
"""
function writeTecplotP1( flname::String, prob::luteos.Problem, mesh::Mesh3D, solution::Matrix{Float64} )

    nelem   = size( mesh.nodes, 3 )
    nnodes  = size( mesh.nodes, 1 )
    ndat    = size( solution, 2 )

    @assert size(solution,1) == mesh.n

    # Open file
    fid = open( flname, "w" )

    ### Write header
    @printf( fid, "TITLE = \"%s\"\n", prob.name )
    @printf( fid, "VARIABLES = \"X\", \"Y\", \"Z\"" )
    for jj in 1:ndat
        @printf( fid, ", \"U<sub>%i</sub>\"",jj )
    end
    @printf( fid, "\n" )
    @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETETRAHEDRON\n", mesh.n, nelem )

    ### Write values
    for pp in 1:mesh.n

        wrtval = [mesh.p[pp,:] solution[pp,:] ]

        for ii in 1:ndat+3
            @printf(fid, "%16.15e\t", wrtval[ii] )
        end
        @printf(fid,"\n")

    end

    ### Write connectivity
    for ee in 1:nelem

        for kk in 1:(mesh.dim+1)
            @printf(fid, "%i\t", mesh.t[ee,kk] )
        end
        @printf(fid,"\n")

    end

    close(fid)



end

"""
    writeTecplot( flname::String, prob::Elas, mesh::Mesh3D,
                  uh::Array{Float64,3}, σh::Array{Float64,3};
                  σmises = Array{Float64,3}(0,0,0) )

Outputs tecplot data file to `flname` for 3D structural elasticity problems.
"""
function writeTecplotP1Defl( flname::String, prob::luteos.Problem, mesh::Mesh3D, solution::Matrix{Float64} )

  nelem   = size( mesh.nodes, 3 )
  nnodes  = size( mesh.nodes, 1 )
  ndat    = size( solution, 2 )

  @assert size(solution,1) == mesh.n

  # defl = solution ./ maximum(solution) * maximum( mesh.p ) * 0.25
  defl = solution

  # Open file
  fid = open( flname, "w" )

  ### Write header
  @printf( fid, "TITLE = \"%s\"\n", prob.name )
  @printf( fid, "VARIABLES = \"X\", \"Y\", \"Z\"" )
  for jj in 1:ndat
    @printf( fid, ", \"U<sub>%i</sub>\"",jj )
  end
  @printf( fid, "\n" )
  @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETETRAHEDRON\n",
    mesh.n, nelem )

  ### Write values
  for pp in 1:mesh.n

    wrtval = [mesh.p[pp,:]+defl[pp,:] solution[pp,:] ]

    for ii in 1:ndat+3
      @printf(fid, "%16.15e\t", wrtval[ii] )
    end
    @printf(fid,"\n")

  end

  ### Write connectivity
  for ee in 1:nelem

    for kk in 1:(mesh.dim+1)
      @printf(fid, "%i\t", mesh.t[ee,kk] )
    end
    @printf(fid,"\n")

  end

  close(fid)

end

function convertToNodal( mesh::Mesh3D, prob::luteos.Elas, uh::Array{Float64,3}, σh::Array{Float64,3} )

    @assert size(mesh.nodes,1) == 4

    nelem = size( mesh.nodes, 3 )

    ncounts = fill( 0, mesh.n )

    σhmises = luteos.compMises3D( σh )

    uhvec    = fill( 0., mesh.n, mesh.dim )
    σhvec    = fill( 0., mesh.n, mesh.dim^2 )
    σhmisvec = fill( 0., mesh.n, 1 )

    for kk in 1:nelem

        indsp = mesh.t[kk,:]

        ncounts[indsp] .+= 1

        for ii in 1:(mesh.dim+1)
            for jj in 1:mesh.dim
                uhvec[indsp[ii],jj] += uh[ii,jj,kk]
            end
            for jj in 1:mesh.dim^2
                σhvec[indsp[ii],jj] += σh[ii,jj,kk]
            end
            σhmisvec[indsp[ii],1] += σhmises[ii,1,kk]
        end

    end

    for jj in 1:mesh.n
        uhvec[jj,:]   ./= ncounts[jj]
        σhvec[jj,:]   ./= ncounts[jj]
        σhmisvec[jj,:] ./= ncounts[jj]
    end

    return uhvec, σhvec, σhmisvec

end

# TODO: Need to add surface output
