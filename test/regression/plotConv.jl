# using PyPlot

flname = "errors_Elas_Dirichlet2D.dat"

Ps = 1:3         # Range of polynomial order
Ns = [9, 17, 33] # Range of grid size

data   = readdlm(flname; header=true)

data = data[1]

Err_uh1 = fill( 0.0, length(Ps), length(Ns) )
Err_uh2 = fill( 0.0, length(Ps), length(Ns) )
Err_σh1 = fill( 0.0, length(Ps), length(Ns) )
Err_σh2 = fill( 0.0, length(Ps), length(Ns) )
Err_σh4 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh1 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh2 = fill( 0.0, length(Ps), length(Ns) )
Err_ϵh4 = fill( 0.0, length(Ps), length(Ns) )

kk = 1
for ii in 1:length(Ps), jj in 1:length(Ns)
  Err_uh1[ii,jj] = data[kk,3]
  Err_uh2[ii,jj] = data[kk,4]
  Err_σh1[ii,jj] = data[kk,5]
  Err_σh2[ii,jj] = data[kk,6]
  Err_σh4[ii,jj] = data[kk,7]
  Err_ϵh1[ii,jj] = data[kk,8]
  Err_ϵh2[ii,jj] = data[kk,9]
  Err_ϵh4[ii,jj] = data[kk,10]
  kk += 1
end

h = 1 ./ ( Ns - 1 )
conv_uh1 = (log.( Err_uh1[:,end-2]) - log.( Err_uh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_uh2 = (log.( Err_uh2[:,end-2]) - log.( Err_uh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

conv_σh1 = (log.( Err_σh1[:,end-2]) - log.( Err_σh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_σh2 = (log.( Err_σh2[:,end-2]) - log.( Err_σh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_σh4 = (log.( Err_σh4[:,end-2]) - log.( Err_σh4[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

conv_ϵh1 = (log.( Err_ϵh1[:,end-2]) - log.( Err_ϵh1[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_ϵh2 = (log.( Err_ϵh2[:,end-2]) - log.( Err_ϵh2[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));
conv_ϵh4 = (log.( Err_ϵh4[:,end-2]) - log.( Err_ϵh4[:,end] ) ) / (log.( h[end-2]) - log.( h[end] ));

clf()
hu1 = loglog(h,Err_uh1)
hu2 = loglog(h,Err_uh2)

hσ1 = loglog(h,Err_σh1)
hσ2 = loglog(h,Err_σh2)
hσ4 = loglog(h,Err_σh4)
