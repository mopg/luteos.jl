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
    writeTecplot( flname::String, prob::Elas, mesh::Mesh2D,
                  uh::Array{Float64,3}, σh::Array{Float64,3};
                  σmises = Array{Float64,3}(0,0,0) )

Outputs tecplot data file to `flname` for 2D structural elasticity problems.
"""
function writeTecplot( flname::String, prob::Elas, mesh::Mesh2D, uh::Array{Float64,3},
  σh::Array{Float64,3}; σmises = Array{Float64,3}(0,0,0) )

  nelem   = size( mesh.nodes, 3 )
  nnodes  = size( mesh.nodes, 1 )
  ntriang = size( mesh.tloc,  1 ) # local number of triangles

  nnodesTot = nnodes * nelem

  nunkTot = mesh.dim + mesh.dim + mesh.dim^2 + ( size(σmises,1) > 0 )

  # Open file
  fid = open( flname, "w" )

  ### Write header
  @printf( fid, "TITLE = \"%s\"\n", prob.name )
  @printf( fid, "VARIABLES = \"X\", \"Y\"\n" )
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
  @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETRIANGLE\n",
    nnodesTot, nelem * ntriang )

  ### Write values
  for pp in 1:nelem

    wrtval = fill( 0.0, nnodes, nunkTot )

    if size( σmises, 1 ) == 0
        wrtval = [mesh.nodes[:,:,pp] uh[:,:,pp] σh[:,:,pp] ]
    else
        wrtval = [mesh.nodes[:,:,pp] uh[:,:,pp] σh[:,:,pp] σmises[:,:,pp] ]
    end

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
    writeTecplot( flname::String, prob::Elas, mesh::Mesh3D,
                  uh::Array{Float64,3}, σh::Array{Float64,3};
                  σmises = Array{Float64,3}(0,0,0) )

Outputs tecplot data file to `flname` for 3D structural elasticity problems.
"""
function writeTecplot( flname::String, prob::Elas, mesh::Mesh3D, uh::Array{Float64,3},
  σh::Array{Float64,3}; σmises = Array{Float64,3}(0,0,0) )

  nelem   = size( mesh.nodes, 3 )
  nnodes  = size( mesh.nodes, 1 )
  ntets   = size( mesh.tloc,  1 ) # local number of tetrahedrons

  nnodesTot = nnodes * nelem

  nunkTot = mesh.dim + mesh.dim + mesh.dim^2 + ( size(σmises,1) > 0 )

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
    @printf( fid, ", \"<greek>σ</greek><sub>mises</sub>\"" )
  end
  @printf( fid, "\n" )
  @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETETRAHEDRON\n",
    nnodesTot, nelem * ntets )

  ### Write values
  for pp in 1:nelem

    wrtval = fill( 0.0, nnodes, nunkTot )

    if size( σmises, 1 ) == 0
        wrtval = [mesh.nodes[:,:,pp] uh[:,:,pp] σh[:,:,pp] ]
    else
        wrtval = [mesh.nodes[:,:,pp] uh[:,:,pp] σh[:,:,pp] σmises[:,:,pp] ]
    end

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

"""
    writeTecplotP1( flname::String, prob::Elas, mesh::Mesh3D,
                  uh::Array{Float64,3}, σh::Array{Float64,3};
                  σmises = Array{Float64,3}(0,0,0) )

Outputs tecplot data file to `flname` for 3D structural elasticity problems.
"""
function writeTecplotP1( flname::String, prob::Elas, mesh::Mesh3D, uh::Array{Float64,3},
                         σh::Array{Float64,3}, σmises::Array{Float64,3} )

  nelem   = size( mesh.nodes, 3 )
  nnodes  = size( mesh.nodes, 1 )

  nnodesTot = mesh.n

  nunkTot = mesh.dim + mesh.dim + mesh.dim^2 + 1

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
  @printf( fid, ", \"<greek>s</greek><sub>mises</sub>\"" )
  @printf( fid, "\n" )
  @printf( fid, "ZONE, DATAPACKING=POINT, NODES=%i, ELEMENTS=%i, ZONETYPE=FETETRAHEDRON\n",
    nnodesTot, nelem )

  ### Interpolate values on nodes
  nelpoin = fill( 0,  mesh.n )
  values  = fill( 0., mesh.n, nunkTot )
  values[:, 1:mesh.dim] = mesh.p

  for el in 1:nelem, jj in 1:(mesh.dim+1)

      p0 = mesh.t[el,jj]
      values[ p0, mesh.dim+(1:mesh.dim)   ] += uh[jj,:,el]
      values[ p0, 2*mesh.dim+(1:mesh.dim^2) ] += σh[jj,:,el]
      values[ p0, end ] += σmises[jj,1,el]
      nelpoin[p0] += 1

  end

  for pp in 1:mesh.n

    values[pp,(mesh.dim+1):end] ./= nelpoin[pp]

    for ii in 1:nunkTot
      @printf(fid, "%16.15e\t", values[pp,ii] )
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
