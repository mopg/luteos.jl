# ---------------------------------------------------------------------------- #
#
#   writeTecplotCDR.jl
#
#   Several functions to write Tecplot output for CDR problems
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    writeTecplot( flname::String, prob::CDR, mesh::Mesh2D,
                  uh::Array{Float64,3}, qh::Array{Float64,3} )

Outputs tecplot data file to `flname` for 3D CDR problems.
"""
function writeTecplot( flname::String, prob::CDR, mesh::Mesh2D,
                       uh::Array{Float64,3}, qh::Array{Float64,3} )

  nelem   = size( mesh.nodes, 3 )
  nnodes  = size( mesh.nodes, 1 )
  ntriang = size( mesh.tloc,  1 ) # local number of triangles

  nnodesTot = nnodes * nelem

  nunkTot = mesh.dim + 1 + mesh.dim

  # Open file
  fid = open( flname, "w" )

  ### Write header
  @printf( fid, "TITLE = \"%s\"\n", prob.name )
  @printf( fid, "VARIABLES = \"X\", \"Y\"\n" )
  @printf( fid, ", \"U\"" )
  for jj in 1:mesh.dim
    @printf( fid, ", \"Q<sub>%i</sub>\"",jj )
  end
  @printf( fid, "\n" )
  @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETRIANGLE\n",
    nnodesTot, nelem * ntriang )

  ### Write values
  for pp in 1:nelem

    wrtval = fill( 0.0, nnodes, nunkTot )

    wrtval = [mesh.nodes[:,:,pp] uh[:,:,pp] qh[:,:,pp] ]

    for jj in 1:nnodes
      for ii in 1:nunkTot
        @printf(fid, "%16.15e\t", wrtval[jj,ii] )
      end
      @printf(fid,"\n")
    end
  end

  ### Write connectivity
  tnode = 0
  for pp in 1:nelem

    for jj in 1:ntriang
      for kk in 1:(mesh.dim+1)
        @printf(fid, "%i\t", tnode + mesh.tloc[jj,kk] )
      end
      @printf(fid,"\n")
    end

    tnode = tnode + nnodes

  end

  close(fid)

end

"""
    writeTecplot( prob::CDR, flname::String, mesh::Mesh3D,
                  uh::Array{Float64,3}, σh::Array{Float64,3},
                  ϵh::Array{Float64,3}; σmises = Array{Float64,3}(0,0,0) )

Outputs tecplot data file to `flname` for 3D CDR problems.
"""
function writeTecplotCD( prob::CDR, flname::String, mesh::Mesh3D, uh::Array{Float64,3},
  qh::Array{Float64,3} )

  nelem   = size( mesh.nodes, 3 )
  nnodes  = size( mesh.nodes, 1 )
  ntets   = size( mesh.tloc,  1 ) # local number of tetrahedrons

  nnodesTot = nnodes * nelem

  nunkTot = mesh.dim + 1 + mesh.dim

  # Open file
  fid = open( flname, "w" )

  ### Write header
  @printf( fid, "TITLE = \"%s\"\n", prob.name )
  @printf( fid, "VARIABLES = \"X\", \"Y\", \"Z\"" )
  @printf( fid, ", \"U\"" )
  for jj in 1:mesh.dim
    @printf( fid, ", \"Q<sub>%i</sub>\"",jj )
  end
  @printf( fid, "\n" )
  @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETETRAHEDRON\n",
    nnodesTot, nelem * ntets )

  ### Write values
  for pp in 1:nelem

    wrtval = fill( 0.0, nnodes, nunkTot )

    wrtval = [mesh.nodes[:,:,pp] uh[:,:,pp] qh[:,:,pp] ]

    for jj in 1:nnodes
      for ii in 1:nunkTot
        @printf(fid, "%16.15e\t", wrtval[jj,ii] )
      end
      @printf(fid,"\n")
    end
  end

  ### Write connectivity
  tnode = 0::Int64
  for pp in 1:nelem

    for jj in 1:ntets
      for kk in 1:(mesh.dim+1)
        @printf(fid, "%i\t", tnode + mesh.tloc[jj,kk] )
      end
      @printf(fid,"\n")
    end

    tnode = tnode + nnodes

  end

  close(fid)

end

# TODO: Need to add surface output
