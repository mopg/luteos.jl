type Material

  name::String   # Name of problem
  E::Float64     # Young's modulus
  ν::Float64     # Poisson's ratio
  ρ::Float64     # Material density
  λ::Float64     # Lamé parameter
  μ::Float64     # Shear modulus (or rigidity)

  σt::Float64    # Maximum tensile stress at yield
  σc::Float64    # Maximum compressive stress at yield

  Cstiff::Array{Array{Float64}} # Stiffness matrices (for 2 and 3D)

end

function Material( name::String )

  if name == "PLA"
    E  = 3.5e9  # Pa
    ν  = 0.36
    ρ  = 4.43e3 # kg/m^3
    σt = 50e6   # Pa
    σc = 60e6   # Pa
  elseif name == "Titanium"
    E  = 113.8e9 # Pa
    ν  = 0.33
    ρ  = 1.24e3  # kg/m^3
    σt = 880e6   # Pa
    σc = 970e6   # Pa
  elseif name == "Al6061"
    E  = 68.9e9 # Pa
    ν  = 0.33
    ρ  = 2.7e3  # kg/m^3
    σt = 276e6  # Pa
    σc = 386e6  # Pa
  end

  λ = E * ν / ( (1 + ν) * (1 - 2*ν) )
  μ = E / ( 2 * (1 + ν) )

  Cstiff = Array{Array{Float64}}(3)

  for dd in 1:3
    Cstiff[dd] = compCstiff( λ, μ, dd)
  end

  Material( name, E, ν, ρ, λ, μ, σt, σc, Cstiff )

end

function Material(; name="Unobtainium", E = 68.9e9, ν = 0.33, ρ = 2.7e3, σt = 276e6, σc = 386e6  )

  λ = E * ν / ( (1 + ν) * (1 - 2*ν) )
  μ = E / ( 2 * (1 + ν) )

  Cstiff = Array{Array{Float64}}(3)

  for dd in 1:3
    Cstiff[dd] = compCstiff( λ, μ, dd)
  end

  Material( name, E, ν, ρ, λ, μ, σt, σc, Cstiff )

end

function compCstiff( λ::Float64, μ::Float64, dim::Int64 )

  Cs = Array{Float64}(dim, dim, dim, dim)

  for ii in 1:dim, jj in 1:dim, kk in 1:dim, ll in 1:dim
    Cs[ii,jj,kk,ll] = λ*kdelt(ii,jj)*kdelt(kk,ll) +
        μ * (kdelt(ii,kk)*kdelt(jj,ll) + kdelt(ii,ll)*kdelt(jj,kk))
  end

  return Cs

end

function kdelt( ii::Int64, jj::Int64 )
  # Return kronecker δ
  return Int64( ii == jj )
end
