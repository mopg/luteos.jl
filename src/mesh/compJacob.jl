function compJacob!( mesh::Mesh2D, master::Master2D )

  xxi = master.dphi[:,:,1]' * mesh.nodes[:,1,:]
  xet = master.dphi[:,:,2]' * mesh.nodes[:,1,:]
  yxi = master.dphi[:,:,1]' * mesh.nodes[:,2,:]
  yet = master.dphi[:,:,2]' * mesh.nodes[:,2,:]

  jac = xxi.*yet - xet.*yxi
  jcw = diagm( master.gwts ) * jac

  xix =  1./jac .* yet
  etx = -1./jac .* yxi
  xiy = -1./jac .* xet
  ety =  1./jac .* xxi

  xix_vec = fill(0.0, size(xxi,1), size(xxi,2), 4)
  xix_vec[:,:,1] = xix
  xix_vec[:,:,2] = xiy
  xix_vec[:,:,3] = etx
  xix_vec[:,:,4] = ety

  mesh.jcw = jcw
  mesh.xix = xix_vec

end
