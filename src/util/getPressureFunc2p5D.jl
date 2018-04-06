function pressureFunc2p5D( cpfile::String, Vinf::Float64, ρinf::Float64,
                           Λ::Float64,  λ::Float64, AR::Float64,
                           ea::Float64, cr::Float64 )

    data = readdlm( cpfile )

    s = data[:,1]
    inds = find( s[1:end-1] - s[2:end] .> 0.5 )
    ind = inds[1]

    xup = data[     1:ind,1 ]
    xlo = data[ ind+1:end,1 ]
    cpup = data[     1:ind, 5 ]
    cplo = data[ ind+1:end, 5 ]

    if xup[1] > xup[2]
        inds2 = find( xup .== minimum(xup) )
        ind2 = inds2[1]
        xlo  = vcat( xup[ind2:-1:2],  xlo )
        cplo = vcat( cpup[ind2:-1:2], cplo )
        deleteat!( xup,  1:(ind2-1) )
        deleteat!( cpup, 1:(ind2-1) )
    end

    xup = vcat( -1e-4, xup )
    xlo = vcat( -1e-4, xlo )

    cpup = vcat( cpup[1], cpup )
    cplo = vcat( cplo[1], cplo )

    spl_up = Dierckx.Spline1D( xup, cpup )
    spl_lo = Dierckx.Spline1D( xlo, cplo )

    pu( p::Matrix{Float64} ) = hcat( 0.5 * ρinf * Vinf^2 *
                                     getCpwing2p5D(p[:,1:2], Λ, λ, AR, ea, cr, spl_up), 0.0*p[:,2], 0.0*p[:,2] )
    pl( p::Matrix{Float64} ) = hcat( 0.5 * ρinf * Vinf^2 *
                                     getCpwing2p5D(p[:,1:2], Λ, λ, AR, ea, cr, spl_lo), 0.0*p[:,2], 0.0*p[:,2] )

    return (pu, pl)

end

function getCpwing2p5D( p2d::Matrix{Float64}, Λ::Float64, λ::Float64, AR::Float64,
                        ea::Float64, cr::Float64, spl::Dierckx.Spline1D )

    b2    = AR * ( cr + λ ) / 4
    ct    = cr * λ
    c     = cr - p2d[:,2] / b2 * (cr - ct)
    xnorm = (p2d[:,1] - p2d[:,2]*tand(Λ) ) ./ c + ea

    # if maximum(xnorm) >= 1.0
    #     println( "max(xnorm) ", maximum(xnorm) )
    #     ii = find( xnorm .== maximum(xnorm) )
    #     println( "x org ", p2d[ii,:])
    # end
    # if minimum(xnorm) <= 0.0
    #     println( "min(xnorm) ", minimum(xnorm) )
    #     ii = find( xnorm .== minimum(xnorm) )
    #     println( "x org ", p2d[ii,:])
    # end
    # assert( maximum(xnorm) >= 1.0 )
    # assert( minimum(xnorm) <= 0.0 )

    cpc = spl( xnorm ) .* c .* cosd(Λ)

    return cpc

end
