function compMises3D( σh::Array{Float64,3} )

    σmis = Array{Float64}( size(σh,1), 1, size(σh,3) )
    σmis[:,1,:] = sqrt.( 0.5 * ( (σh[:,1,:] - σh[:,5,:]).^2
                             + (σh[:,5,:] - σh[:,9,:]).^2
                             + (σh[:,9,:] - σh[:,1,:]).^2
                             + 6*(σh[:,2,:].^2 + σh[:,6,:].^2 + σh[:,7,:].^2) ) )

    return σmis

end

function compMises2D( σh::Array{Float64,3} )

    σmis = Array{Float64}( size(σh,1), 1, size(σh,3) )
    σmis[:,1,:] = sqrt.( (σh[:,1,:] - σh[:,4,:]).^2
                         + 3*σh[:,2,:].^2 )

    return σmis

end
