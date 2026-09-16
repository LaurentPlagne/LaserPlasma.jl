# Small fixed-size 3-vector on top of `NTuple{3,Float64}`. Wrapping the tuple
# in a named type (never a type alias) keeps `==` local to this type, and the
# methods keep the C operation order component by component.

"""
    Vec3
    Vec3(x, y, z)

Small fixed-size 3-vector on top of `NTuple{3,Float64}`. Wrapping the tuple in
a named type (never a type alias) keeps `==` local to this type; `+`, `-` and
scalar `*` are defined below in the C operation order, component by component.
"""
struct Vec3
    data::NTuple{3,Float64}
end

Vec3(x::Real, y::Real, z::Real) = Vec3((Float64(x), Float64(y), Float64(z)))

@inline Base.getindex(v::Vec3, i::Int) = v.data[i]
@inline Base.:+(a::Vec3, b::Vec3) = Vec3(a.data .+ b.data)
@inline Base.:-(a::Vec3, b::Vec3) = Vec3(a.data .- b.data)
@inline Base.:*(a::Real, v::Vec3) = Vec3(a .* v.data)
@inline Base.:*(v::Vec3, a::Real) = a * v

"""
    squared_norm(v)

`x*x + y*y + z*z`, in the C accumulation order.
"""
@inline squared_norm(v::Vec3) = v.data[1] * v.data[1] + v.data[2] * v.data[2] +
                                v.data[3] * v.data[3]

"""
    dot(a, b)

Dot product `a · b`, components in the C order.
"""
@inline dot(a::Vec3, b::Vec3) = a.data[1] * b.data[1] + a.data[2] * b.data[2] +
                                a.data[3] * b.data[3]

"""
    cross(a, b)

Cross product `a × b`, with the C component expressions.
"""
@inline cross(a::Vec3, b::Vec3) = Vec3((a.data[2] * b.data[3] - a.data[3] * b.data[2],
                                        a.data[3] * b.data[1] - a.data[1] * b.data[3],
                                        a.data[1] * b.data[2] - a.data[2] * b.data[1]))
