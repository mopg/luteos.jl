abstract type A

end

# include("temp2.jl")

type B <: A
    x
    y
end

type C <: A
    x
    z
end

prop1(b::B) = b.x
prop2(b::B) = b.y
prop1(c::C) = c.w
prop2(c::C) = c.z  # changed from prop2(c::C)=c.w

mysum(a::A) = prop1(a) + prop2(a)
function getx(a::A)
  b = a.x
  return b
end
