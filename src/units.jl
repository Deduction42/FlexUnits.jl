"""
    map_dimensions(f::F, args::AbstractDimensions...)

Similar to the `map` function, but specifically iterates over dimensions
"""
function map_dimensions(f::F, args::AbstractDimensions...) where {F<:Function}
    D = promote_type(typeof(args).parameters...)
    return  D(
        ( f((getproperty(arg, dim) for arg in args)...) for dim in dimension_names(D) )...
    )
end

#=============================================================================================
 Operations on dimensions
=============================================================================================#
Base.:*(args::AbstractDimensions...) = map_dimensions(+, args...)
Base.:/(args::AbstractDimensions...) = map_dimensions(-, args...)
Base.inv(arg::AbstractDimensions) = map_dimensions(-, arg)
Base.:^(d::AbstractDimensions, p::Integer) = map_dimensions(Base.Fix1(*, p), d)
Base.:^(d::AbstractDimensions{R}, p::Number) where {R} = map_dimensions(Base.Fix1(*, tryrationalize(R, p)), d)
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{p}) where {p} = map_dimensions(Base.Fix1(*, p), d)
Base.sqrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 2))
Base.cbrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 3))
Base.adjoint(d::AbstractDimensions) = d

#=============================================================================================
 Operations on units
=============================================================================================#
function Base.:*(u::U...) where U <: AbstractUnits
    return constructorof(U)(scale=*(map(uscale, u)...), dims=*(map(scalar_dimension, u)...))
end
function Base.:/(u1::U, u2::U) where U <: AbstractUnits
    return constructorof(U)(scale=/(map(uscale, (u1, u2))...), dims=/(map(scalar_dimension, (u1, u2))...))
end
function Base.:inv(arg::U) where U <: AbstractUnits
    return constructorof(U)(scale=inv(uscale(arg)), dims=inv(scalar_dimension(arg)))
end
function Base.:^(u::U, p::Number) where U <:AbstractUnits
    return constructorof(U)(scale=uscale(u)^p, scalar_dimension(u)^p)
end
Base.sqrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 2))
Base.cbrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 3))
Base.adjoint(u::AbstractUnits) = u

#Equality does not compare symbols
Base.:(==)(u1::AbstractUnits, u2::AbstractUnits) = (uscale(u1) == uscale(u2)) & (uoffset(u1) == uoffset(u2)) & (dimension(u1) == dimension(u2))
Base.:(==)(u1::AbstractScalarUnits, u2::AbstractScalarUnits) = (uscale(u1) == uscale(u2)) & (dimension(u1) == dimension(u2))

#=============================================================================================
 Utility functions
=============================================================================================#
assert_scalar(u::AbstractDimensions)  = u
assert_scalar(u::AbstractScalarUnits) = u
assert_scalar(u::AbstractAffineUnits) = iszero(uoffset(u)) ? u : throw(ScalarUnitError(u))
scalar_dimension(u::AbstractUnitLike) = dimension(assert_scalar(u))


#Test code
#=
velocity = Dimensions(length=1, time=-1)
velocity*velocity
velocity^2
velocity^(-0.5)
=#