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
    writeTecplot( flname::String, mesh::Mesh2D )

Outputs tecplot data file of mesh to `flname` for 2D mesh.
"""
function writeTecplot( flname::String, mesh::Mesh2D )

  nelem   = size( mesh.nodes, 3 )
  nnodes  = size( mesh.nodes, 1 )
  ntriang = size( mesh.tloc,  1 ) # local number of triangles

  nnodesTot = nnodes * nelem

  nunkTot = mesh.dim

  # Open file
  fid = open( flname, "w" )

  ### Write header
  @printf( fid, "TITLE = \"%s\"\n", "Mesh" )
  @printf( fid, "VARIABLES = \"X\", \"Y\"\n" )
  @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETRIANGLE\n",
    nnodesTot, nelem * ntriang )

  ### Write values
  for pp in 1:nelem

    wrtval = fill( 0.0, nnodes, nunkTot )

    wrtval = mesh.nodes[:,:,pp]

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
    writeTecplot( flname::String, mesh::Mesh3D )

Outputs tecplot data file to `flname` for 3D meshes.
"""
function writeTecplot( flname::String, mesh::Mesh3D )

  nelem   = size( mesh.nodes, 3 )
  nnodes  = size( mesh.nodes, 1 )
  ntets   = size( mesh.tloc,  1 ) # local number of tetrahedrons

  nnodesTot = nnodes * nelem

  nunkTot = mesh.dim

  # Open file
  fid = open( flname, "w" )

  ### Write header
  @printf( fid, "TITLE = \"%s\"\n", "Mesh" )
  @printf( fid, "VARIABLES = \"X\", \"Y\", \"Z\"" )
  @printf( fid, "\n" )
  @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETETRAHEDRON\n",
    nnodesTot, nelem * ntets )

  ### Write values
  for pp in 1:nelem

    wrtval = fill( 0.0, nnodes, nunkTot )

    wrtval = mesh.nodes[:,:,pp]

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
