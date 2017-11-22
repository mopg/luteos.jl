# ---------------------------------------------------------------------------- #
#
#   basisFuncTetLag.jl
#
#   Get Lagrangian basis functions for Tetrahedron
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

function basisFuncTetLag( P::P1, s::Vector{Float64}, t::Vector{Float64}, u::Vector{Float64} )

  dim = 3
  nx  = length( s )

  # Linear

  nphi = 4
  phi  = Array{Float64}(nphi, nx)
  dphi = Array{Float64}(nphi, nx, dim)

  # value
  phi[1,:] = 1 - s - t - u # 1 @ node 0 (s = 0, t = 0, u = 0)
  phi[2,:] =     s         # 1 @ node 1 (s = 1, t = 0, u = 0)
  phi[3,:] =         t     # 1 @ node 2 (s = 0, t = 1, u = 0)
  phi[4,:] =             u # 1 @ node 3 (s = 0, t = 0, u = 1)

  # derivative in x
  dphi[1,:,1] = -1
  dphi[2,:,1] =  1
  dphi[3,:,1] =  0
  dphi[4,:,1] =  0

  # derivative in y
  dphi[1,:,2] = -1
  dphi[2,:,2] =  0
  dphi[3,:,2] =  1
  dphi[4,:,2] =  0

  # derivative in z
  dphi[1,:,3] = -1
  dphi[2,:,3] =  0
  dphi[3,:,3] =  0
  dphi[4,:,3] =  1

  return phi, dphi

end


function basisFuncTetLag( P::P2, s::Vector{Float64}, t::Vector{Float64}, u::Vector{Float64} )

  dim = 3
  nx  = length( s )

  # Quadratric
  nphi = 10
  phi  = Array{Float64}(nphi, nx)
  dphi = Array{Float64}(nphi, nx, dim)

  # precompute
  ss = (s.*s)
  tt = (t.*t)
  uu = (u.*u)

  st = s.*t
  su = s.*u
  tu = t.*u

  # value
  phi[1,:] =  s.*-3.0-t.*3.0-u.*3.0+st.*4.0+su.*4.0+tu.*4.0+ss.*2.0+tt.*2.0+uu.*2.0+1.0
  phi[2,:] =  -s+ss.*2.0
  phi[3,:] =  -t+tt.*2.0
  phi[4,:] =  -u+uu.*2.0
  phi[5,:] =  tu.*4.0
  phi[6,:] =  su.*4.0
  phi[7,:] =  st.*4.0
  phi[8,:] =  t.*4.0-st.*4.0-tu.*4.0-tt.*4.0
  phi[9,:] =  u.*4.0-su.*4.0-tu.*4.0-uu.*4.0
  phi[10,:] =  s.*4.0-st.*4.0-su.*4.0-ss.*4.0

  # derivative in x
  dphi[1,:,1] =  s.*4.0+t.*4.0+u.*4.0-3.0
  dphi[2,:,1] =  s.*4.0-1.0
  dphi[3,:,1] =  0.0
  dphi[4,:,1] =  0.0
  dphi[5,:,1] =  0.0
  dphi[6,:,1] =  u.*4.0
  dphi[7,:,1] =  t.*4.0
  dphi[8,:,1] =  t.*-4.0
  dphi[9,:,1] =  u.*-4.0
  dphi[10,:,1] =  s.*-8.0-t.*4.0-u.*4.0+4.0

  # derivative in y
  dphi[1,:,2] =  s.*4.0+t.*4.0+u.*4.0-3.0
  dphi[2,:,2] =  0.0
  dphi[3,:,2] =  t.*4.0-1.0
  dphi[4,:,2] =  0.0
  dphi[5,:,2] =  u.*4.0
  dphi[6,:,2] =  0.0
  dphi[7,:,2] =  s.*4.0
  dphi[8,:,2] =  s.*-4.0-t.*8.0-u.*4.0+4.0
  dphi[9,:,2] =  u.*-4.0
  dphi[10,:,2] =  s.*-4.0

  # derivative in z
  dphi[1,:,3] =  s.*4.0+t.*4.0+u.*4.0-3.0
  dphi[2,:,3] =  0.0
  dphi[3,:,3] =  0.0
  dphi[4,:,3] =  u.*4.0-1.0
  dphi[5,:,3] =  t.*4.0
  dphi[6,:,3] =  s.*4.0
  dphi[7,:,3] =  0.0
  dphi[8,:,3] =  t.*-4.0
  dphi[9,:,3] =  s.*-4.0-t.*4.0-u.*8.0+4.0
  dphi[10,:,3] =  s.*-4.0

  # transpose

  return phi, dphi

end

function basisFuncTetLag( P::P3, s::Vector{Float64}, t::Vector{Float64}, u::Vector{Float64} )

  dim = 3
  nx  = length( s )

  # Cubic
  nphi = 20
  phi  = Array{Float64}(nphi, nx)
  dphi = Array{Float64}(nphi, nx, dim)

  ss = (s.*s)
  tt = (t.*t)
  uu = (u.*u)

  st = s.*t
  su = s.*u
  tu = t.*u

  ssu = ss.*u
  sst = ss.*t
  tuu =  t.*uu
  ttu = tt.*u
  suu =  s.*uu
  stt =  s.*tt

  stu =  st.*u

  ttt = (t.*t.*t)
  uuu = (u.*u.*u)
  sss = (s.*s.*s)

  # value
  phi[1,:] =  s.*(-1.1E1/2.0)-t.*(1.1E1/2.0)-u.*(1.1E1/2.0)+st.*1.8E1+su.*1.8E1+tu.*1.8E1-stt.*(2.7E1/2.0)-
             sst.*(2.7E1/2.0)-suu.*(2.7E1/2.0)-ssu.*(2.7E1/2.0)-tuu.*(2.7E1/2.0)-ttu.*(2.7E1/2.0)+(s.*
             s).*9.0-sss.*(9.0/2.0)+tt.*9.0-ttt.*(9.0/2.0)+uu.*9.0-uuu.*(9.0/2.0)-stu.*2.7E1+1.0
  phi[2,:] =  s-ss.*(9.0/2.0)+sss.*(9.0/2.0)
  phi[3,:] =  t-tt.*(9.0/2.0)+ttt.*(9.0/2.0)
  phi[4,:] =  u-uu.*(9.0/2.0)+uuu.*(9.0/2.0)
  phi[5,:] =  tu.*(-9.0/2.0)+ttu.*(2.7E1/2.0)
  phi[6,:] =  tu.*(-9.0/2.0)+tuu.*(2.7E1/2.0)
  phi[7,:] =  su.*(-9.0/2.0)+suu.*(2.7E1/2.0)
  phi[8,:] =  su.*(-9.0/2.0)+ssu.*(2.7E1/2.0)
  phi[9,:] =  st.*(-9.0/2.0)+sst.*(2.7E1/2.0)
  phi[10,:] =  st.*(-9.0/2.0)+stt.*(2.7E1/2.0)
  phi[11,:] =  t.*(-9.0/2.0)+st.*(9.0/2.0)+tu.*(9.0/2.0)-stt.*(2.7E1/2.0)-ttu.*(2.7E1/2.0)+tt.*1.8E1-
             ttt.*(2.7E1/2.0)
  phi[12,:] =  t.*9.0-st.*(4.5E1/2.0)-tu.*(4.5E1/2.0)+stt.*2.7E1+sst.*(2.7E1/2.0)+tuu.*(2.7E1/2.0)+
             ttu.*2.7E1-tt.*(4.5E1/2.0)+ttt.*(2.7E1/2.0)+stu.*2.7E1
  phi[13,:] =  u.*9.0-su.*(4.5E1/2.0)-tu.*(4.5E1/2.0)+suu.*2.7E1+ssu.*(2.7E1/2.0)+tuu.*2.7E1+tt.*
             u.*(2.7E1/2.0)-uu.*(4.5E1/2.0)+uuu.*(2.7E1/2.0)+stu.*2.7E1
  phi[14,:] =  u.*(-9.0/2.0)+su.*(9.0/2.0)+tu.*(9.0/2.0)-suu.*(2.7E1/2.0)-tuu.*(2.7E1/2.0)+uu.*1.8E1-
             uuu.*(2.7E1/2.0)
  phi[15,:] =  s.*9.0-st.*(4.5E1/2.0)-su.*(4.5E1/2.0)+stt.*(2.7E1/2.0)+sst.*2.7E1+suu.*(2.7E1/2.0)+
             ssu.*2.7E1-ss.*(4.5E1/2.0)+sss.*(2.7E1/2.0)+stu.*2.7E1
  phi[16,:] =  s.*(-9.0/2.0)+st.*(9.0/2.0)+su.*(9.0/2.0)-sst.*(2.7E1/2.0)-ssu.*(2.7E1/2.0)+ss.*1.8E1-
             sss.*(2.7E1/2.0)
  phi[17,:] =  stu.*2.7E1
  phi[18,:] =  tu.*2.7E1-tuu.*2.7E1-ttu.*2.7E1-stu.*2.7E1
  phi[19,:] =  su.*2.7E1-suu.*2.7E1-ssu.*2.7E1-stu.*2.7E1
  phi[20,:] =  st.*2.7E1-stt.*2.7E1-sst.*2.7E1-stu.*2.7E1

  # derivative in x
  dphi[1,:,1] =  s.*1.8E1+t.*1.8E1+u.*1.8E1-st.*2.7E1-su.*2.7E1-tu.*2.7E1-ss.*(2.7E1/2.0)-tt.*(2.7E1/2.0)-
             uu.*(2.7E1/2.0)-1.1E1/2.0
  dphi[2,:,1] =  s.*-9.0+ss.*(2.7E1/2.0)+1.0
  dphi[3,:,1] =  0.0
  dphi[4,:,1] =  0.0
  dphi[5,:,1] =  0.0
  dphi[6,:,1] =  0.0
  dphi[7,:,1] =  u.*(-9.0/2.0)+uu.*(2.7E1/2.0)
  dphi[8,:,1] =  u.*(-9.0/2.0)+su.*2.7E1
  dphi[9,:,1] =  t.*(-9.0/2.0)+st.*2.7E1
  dphi[10,:,1] =  t.*(-9.0/2.0)+tt.*(2.7E1/2.0)
  dphi[11,:,1] =  t.*(9.0/2.0)-tt.*(2.7E1/2.0)
  dphi[12,:,1] =  t.*(-4.5E1/2.0)+st.*2.7E1+tu.*2.7E1+tt.*2.7E1
  dphi[13,:,1] =  u.*(-4.5E1/2.0)+su.*2.7E1+tu.*2.7E1+uu.*2.7E1
  dphi[14,:,1] =  u.*(9.0/2.0)-uu.*(2.7E1/2.0)
  dphi[15,:,1] =  s.*-4.5E1-t.*(4.5E1/2.0)-u.*(4.5E1/2.0)+st.*5.4E1+su.*5.4E1+tu.*2.7E1+ss.*(8.1E1/2.0)+tt.*
             (2.7E1/2.0)+uu.*(2.7E1/2.0)+9.0
  dphi[16,:,1] =  s.*3.6E1+t.*(9.0/2.0)+u.*(9.0/2.0)-st.*2.7E1-su.*2.7E1-ss.*(8.1E1/2.0)-9.0/2.0
  dphi[17,:,1] =  tu.*2.7E1
  dphi[18,:,1] =  tu.*-2.7E1
  dphi[19,:,1] =  u.*2.7E1-su.*5.4E1-tu.*2.7E1-uu.*2.7E1
  dphi[20,:,1] =  t.*2.7E1-st.*5.4E1-tu.*2.7E1-tt.*2.7E1

  # derivative in y
  dphi[1,:,2] =  s.*1.8E1+t.*1.8E1+u.*1.8E1-st.*2.7E1-su.*2.7E1-tu.*2.7E1-ss.*(2.7E1/2.0)-tt.*(2.7E1/2.0)-
             uu.*(2.7E1/2.0)-1.1E1/2.0
  dphi[2,:,2] =  0.0
  dphi[3,:,2] =  t.*-9.0+tt.*(2.7E1/2.0)+1.0
  dphi[4,:,2] =  0.0
  dphi[5,:,2] =  u.*(-9.0/2.0)+tu.*2.7E1
  dphi[6,:,2] =  u.*(-9.0/2.0)+uu.*(2.7E1/2.0)
  dphi[7,:,2] =  0.0
  dphi[8,:,2] =  0.0
  dphi[9,:,2] =  s.*(-9.0/2.0)+ss.*(2.7E1/2.0)
  dphi[10,:,2] =  s.*(-9.0/2.0)+st.*2.7E1
  dphi[11,:,2] =  s.*(9.0/2.0)+t.*3.6E1+u.*(9.0/2.0)-st.*2.7E1-tu.*2.7E1-tt.*(8.1E1/2.0)-9.0/2.0
  dphi[12,:,2] =  s.*(-4.5E1/2.0)-t.*4.5E1-u.*(4.5E1/2.0)+st.*5.4E1+su.*2.7E1+tu.*5.4E1+ss.*(2.7E1/2.0)+tt.*
             (8.1E1/2.0)+uu.*(2.7E1/2.0)+9.0
  dphi[13,:,2] =  u.*(-4.5E1/2.0)+su.*2.7E1+tu.*2.7E1+uu.*2.7E1
  dphi[14,:,2] =  u.*(9.0/2.0)-uu.*(2.7E1/2.0)
  dphi[15,:,2] =  s.*(-4.5E1/2.0)+st.*2.7E1+su.*2.7E1+ss.*2.7E1
  dphi[16,:,2] =  s.*(9.0/2.0)-ss.*(2.7E1/2.0)
  dphi[17,:,2] =  su.*2.7E1
  dphi[18,:,2] =  u.*2.7E1-su.*2.7E1-tu.*5.4E1-uu.*2.7E1
  dphi[19,:,2] =  su.*-2.7E1
  dphi[20,:,2] =  s.*2.7E1-st.*5.4E1-su.*2.7E1-ss.*2.7E1

  # derivative in z
  dphi[1,:,3] =  s.*1.8E1+t.*1.8E1+u.*1.8E1-st.*2.7E1-su.*2.7E1-tu.*2.7E1-ss.*(2.7E1/2.0)-tt.*(2.7E1/2.0)-
             uu.*(2.7E1/2.0)-1.1E1/2.0
  dphi[2,:,3] =  0.0
  dphi[3,:,3] =  0.0
  dphi[4,:,3] =  u.*-9.0+uu.*(2.7E1/2.0)+1.0
  dphi[5,:,3] =  t.*(-9.0/2.0)+tt.*(2.7E1/2.0)
  dphi[6,:,3] =  t.*(-9.0/2.0)+tu.*2.7E1
  dphi[7,:,3] =  s.*(-9.0/2.0)+su.*2.7E1
  dphi[8,:,3] =  s.*(-9.0/2.0)+ss.*(2.7E1/2.0)
  dphi[9,:,3] =  0.0
  dphi[10,:,3] =  0.0
  dphi[11,:,3] =  t.*(9.0/2.0)-tt.*(2.7E1/2.0)
  dphi[12,:,3] =  t.*(-4.5E1/2.0)+st.*2.7E1+tu.*2.7E1+tt.*2.7E1
  dphi[13,:,3] =  s.*(-4.5E1/2.0)-t.*(4.5E1/2.0)-u.*4.5E1+st.*2.7E1+su.*5.4E1+tu.*5.4E1+ss.*(2.7E1/2.0)+tt.*
             (2.7E1/2.0)+uu.*(8.1E1/2.0)+9.0
  dphi[14,:,3] =  s.*(9.0/2.0)+t.*(9.0/2.0)+u.*3.6E1-su.*2.7E1-tu.*2.7E1-uu.*(8.1E1/2.0)-9.0/2.0
  dphi[15,:,3] =  s.*(-4.5E1/2.0)+st.*2.7E1+su.*2.7E1+ss.*2.7E1
  dphi[16,:,3] =  s.*(9.0/2.0)-ss.*(2.7E1/2.0)
  dphi[17,:,3] =  st.*2.7E1
  dphi[18,:,3] =  t.*2.7E1-st.*2.7E1-tu.*5.4E1-tt.*2.7E1
  dphi[19,:,3] =  s.*2.7E1-st.*2.7E1-su.*5.4E1-ss.*2.7E1
  dphi[20,:,3] =  st.*-2.7E1

  return phi, dphi

end
