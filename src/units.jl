
# Operations on core dimensions
Base.:*(args::AbstractDimensions...) = map_dimensions(+, args...)
Base.:/(args::AbstractDimensions...) = map_dimensions(-, args...)
Base.inv(arg::AbstractDimensions) = map_dimensions(-, arg)

# Operations on scalar-like dimensions
function Base.:*(args::U...) where U <: AbstractUnitLike
    return constructorof(U)(
        scale = *(map(uscale, args)...), 
        dims  = *(map(scalar_dimension, args)...)
    )
end

function Base.:/(arg1::U, arg2::U) where U <: AbstractUnitLike
    return constructorof(U)(
        scale = /(map(uscale, (arg1, arg2))...), 
        dims  = /(map(scalar_dimension, (arg1, arg2))...)
    )
end

Base.:(==)(u1::AbstractUnits, u2::AbstractUnits) = all((uscale(u1) == uscale(u2), uoffset(u1) == uoffset(u2), dimension(u1) == dimension(u2)))
Base.:(==)(u1::AbstractScalarUnits, u2::AbstractScalarUnits) = all((uscale(u1) == uscale(u2), dimension(u1) == dimension(u2)))

assert_scalar(u::AbstractDimensions)  = u
assert_scalar(u::AbstractScalarUnits) = u
assert_scalar(u::AbstractAffineUnits) = iszero(uoffset(u)) ? u : throw(ScalarUnitError(u))
scalar_dimension(u::AbstractUnitLike) = dimension(assert_scalar(u))


