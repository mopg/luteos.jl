function pressureFunc( cpfile::String, Vinf::Float64, ρinf::Float64,
                       wrat::Float64, c::Float64 )

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

    pu( p::Matrix{Float64} ) = hcat( wrat * 0.5 * ρinf * Vinf^2 * spl_up( p[:,1]/c ), 0.0*p[:,2] )
    pl( p::Matrix{Float64} ) = hcat( wrat * 0.5 * ρinf * Vinf^2 * spl_lo( p[:,1]/c ), 0.0*p[:,2] )

    return (pu, pl)

end
