#Arithmetic Types in Base that support +-*/ (use instead of Any to avoid ambiguity)
const ARITHMETICS = Union{Number, AbstractArray, Missing}

#=============================================================================================
 Mathematical operations on dimensions
=============================================================================================#
"""
    map_dimensions(f::F, args::AbstractDimensions...)

Similar to the `map` function, but specifically iterates over dimensions.
Useful for defining mathematical operations for dimensions
"""
@inline function map_dimensions(f::F, args::AbstractDimensions...) where {F<:Function}
    D = promote_type(typeof(args).parameters...)
    dnames = dimension_names(D)
    return  D(
        ( f((getproperty(arg, name) for arg in args)...) for name in  dnames)...
    )
end

@inline Base.:+(arg1::AbstractDimensions) = arg1
@inline Base.:+(arg1::AbstractDimensions, arg2::AbstractDimensions) = firstequal(arg1, arg2)
@inline Base.:-(arg1::AbstractDimensions) = arg1
@inline Base.:-(arg1::AbstractDimensions, arg2::AbstractDimensions) = firstequal(arg1, arg2)
Base.:*(arg1::AbstractDimensions, arg2::AbstractDimensions) = map_dimensions(+, arg1, arg2)
Base.:/(arg1::AbstractDimensions, arg2::AbstractDimensions) = map_dimensions(-, arg1, arg2)
Base.inv(arg::AbstractDimensions) = map_dimensions(-, arg)
Base.:^(d::AbstractDimensions, p::Integer) = map_dimensions(Base.Fix1(*, p), d)
Base.:^(d::AbstractDimensions{R}, p::Real) where {R} = map_dimensions(Base.Fix1(*, R(dimensionless(p))), d)
Base.sqrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 2))
Base.cbrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 3))
Base.abs2(d::AbstractDimensions) = d^2
Base.adjoint(d::AbstractDimensions) = inv(d) #Adjoint of units is an inverse (unitful applications of x=A*y)

@inline Base.literal_pow(::typeof(^), d::D, ::Val{0}) where {D <: AbstractDimensions} = D()
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{1}) = d 
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{2}) = d*d 
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{3}) = d*d*d
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{-1}) = inv(d) 
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{-2}) = inv(d*d)

#=============================================================================================
 Mathematical operations on abstract units (mostly for parsing)
=============================================================================================#
Base.:*(x::ARITHMETICS, u::AbstractUnitLike) = quantity(x, u)
Base.:/(x::ARITHMETICS, u::AbstractUnitLike) = quantity(x, inv(u))

Base.:*(q::AbstractQuantity, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)*u)
Base.:*(u::AbstractUnitLike, q::AbstractQuantity) = q*u
Base.:/(q::AbstractQuantity, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)/u)
Base.:/(u::AbstractUnitLike, q::AbstractQuantity) = quantity(inv(ustrip(q)), u/unit(q))

function Base.:*(u1::U, u2::U) where U <: AbstractUnits
    return constructorof(U)(scale=*(map(uscale, (u1, u2))...), dims=*(map(scalar_dimension, (u1, u2))...))
end

function Base.:/(u1::U, u2::U) where U <: AbstractUnits
    return constructorof(U)(scale=/(map(uscale, (u1, u2))...), dims=/(map(scalar_dimension, (u1, u2))...))
end

function Base.:inv(arg::U) where U <: AbstractUnits
    return constructorof(U)(scale=inv(uscale(arg)), dims=inv(scalar_dimension(arg)))
end

function Base.:^(u::U, p::Real) where U <:AbstractUnits
    pn = dimensionless(p)
    return constructorof(U)(scale=uscale(u)^pn, dims=scalar_dimension(u)^pn)
end

Base.:*(u1::AbstractUnitLike, u2::AbstractUnitLike) = *(promote(u1,u2)...)
Base.:/(u1::AbstractUnitLike, u2::AbstractUnitLike) = /(promote(u1,u2)...)

Base.sqrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 2))
Base.cbrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 3))
Base.adjoint(u::AbstractUnits) = inv(u)

#Equality does not compare symbols
Base.:(==)(u1::AbstractAffineUnits, u2::AbstractAffineUnits) = (uscale(u1) == uscale(u2)) & (uoffset(u1) == uoffset(u2)) & (dimension(u1) == dimension(u2))

@inline firstequal(arg1::AbstractUnitLike) = arg1

function firstequal(arg1::AbstractUnitLike, arg2::AbstractUnitLike)
    (arg1, arg2) = promote(arg1, arg2)
    return (arg1 == arg2) ? arg1 : throw(DimensionError((arg1,arg2)))
end

function firstequal(arg1::AbstractUnitLike, arg2::AbstractUnitLike, arg3::AbstractUnitLike, argN::AbstractUnitLike...)
    newargs = promote(arg1, arg2, arg3, argN...)
    return allequal(newargs) ? first(newargs) : throw(DimensionError(newargs))
end


#=============================================================================================
 Mathematical operations on quantities
=============================================================================================#
"""
    with_ubase(f, args::AbstractQuantity...)

Converts all arguments to base units, and applies `f` to values and dimensions
returns the narrowest possible quanity. Useful for defining 
new functions for quantities 
    
Thus, in order to support `f` for quantities, simply define 

f(dims::AbstractDimensions...)
f(args::AbstractQuantity...) = with_ubase(f, args...)
"""
@inline function with_ubase(f, args::AbstractQuantity...)
    baseargs = map(ubase, args)
    basevals = map(ustrip, baseargs)
    basedims = map(unit, baseargs)
    return quantity(f(basevals...), f(basedims...))
end


function Base.:(==)(q1::AbstractQuantity, q2::AbstractQuantity)
    qb1 = ubase(q1)
    qb2 = ubase(q2)
    return (ustrip(qb1) == ustrip(qb2)) && (unit(qb1) == unit(qb2))
end

function Base.:(≈)(q1::AbstractQuantity, q2::AbstractQuantity)
    qb1 = ubase(q1)
    qb2 = ubase(q2)
    return (ustrip(qb1) ≈ ustrip(qb2)) && (unit(qb1) == unit(qb2))
end

Base.:+(q::AbstractQuantity, x::ARITHMETICS) = dimensionless(q) + x 
Base.:+(x::ARITHMETICS, q::AbstractQuantity) = dimensionless(q) + x
Base.:+(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(+, q1, q2)
Base.:+(q1::AbstractQuantity, qN::AbstractQuantity...) = with_ubase(+, q1, qN...)

Base.:-(q::AbstractQuantity, x::ARITHMETICS) = dimensionless(q) - x 
Base.:-(x::ARITHMETICS, q::AbstractQuantity) = x - dimensionless(q)
Base.:-(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(-, q1, q2)
Base.:-(q1::AbstractQuantity) = with_ubase(-, q1)

Base.:*(q0::AbstractQuantity, x::ARITHMETICS) = (q = ubase(q0); quantity(ustrip(q)*x, unit(q)))
Base.:*(x::ARITHMETICS, q0::AbstractQuantity) = (q = ubase(q0); quantity(ustrip(q)*x, unit(q)))
Base.:*(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(*, q1, q2)
Base.:*(q1::AbstractQuantity, qN::AbstractQuantity...) = with_ubase(*, q1, qN...)

Base.:/(q0::AbstractQuantity, x::ARITHMETICS) = (q = ubase(q0); quantity(ustrip(q)/x, unit(q)))
Base.:/(x::ARITHMETICS, q0::AbstractQuantity) = (q = ubase(q0); quantity(x/ustrip(q), inv(unit(q))))
Base.:/(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(/, q1, q2)
Base.:inv(q::AbstractQuantity) = with_ubase(inv, q)
Base.adjoint(q::AbstractQuantity) = with_ubase(adjoint, q)

Base.:^(q::AbstractQuantity, p::Number) = with_ubase(Base.Fix2(^, p), q)
Base.:^(q::Number, p::AbstractQuantity) = q^dimensionless(p)
Base.:^(q::AbstractQuantity, p::AbstractQuantity) = with_ubase(Base.Fix2(^, dimensionless(p)), q)

@inline Base.literal_pow(::typeof(^), q::AbstractQuantity, ::Val{p}) where {p} = with_ubase(x->Base.literal_pow(^, x, Val(dimensionless(p))), q)

Base.sqrt(q::AbstractQuantity) = with_ubase(sqrt, q)
Base.cbrt(q::AbstractQuantity) = with_ubase(cbrt, q)
Base.abs2(q::AbstractQuantity) = with_ubase(abs2, q)

Base.zero(::Type{D}) where D<:AbstractDimensions = D()
Base.zero(::Type{U}) where U<:AffineUnits = U(dims=dimtype(U)())

Base.zero(::Type{<:AbstractQuantity{T, D}}) where {T, D<:AbstractDimensions} = quantity(zero(T), zero(D))
Base.zero(::Type{<:AbstractQuantity{T, <:AbstractUnits{D}}}) where {T, D<:AbstractDimensions} = quantity(zero(T), zero(D))
Base.one(::Type{<:AbstractQuantity{T}}) where T = one(T) #unitless
Base.oneunit(::Type{<:AbstractQuantity{T,D}}) where {T,D} = quantity(one(T), zero(D)) #unitless with type

#Comparison functions
for f in (:<, :<=, :isless)
    @eval function Base.$f(q1::AbstractQuantity, q2::AbstractQuantity)
        (b1, b2) = (ubase(q1), ubase(q2))
        unit(b1) == unit(b2) || throw(DimensionError(unit(b1), unit(b2)))
        return $f(ustrip(b1), ustrip(b2))
    end
end

#Functions that return the same unit
for f in (
        :float, :abs, :real, :imag, :conj, :significand, :zero, :oneunit, :typemax, :transpose
    )
    @eval Base.$f(u::AbstractDimensions) = u
    @eval Base.$f(q::AbstractQuantity) = with_ubase($f, q)
end

#Single-argument funtions that require dimensionless input and produce dimensionless output
#Note that dimensionless input for angles are "radians" so functiosn like "sind" don't apply 
for f in (
        :sin, :cos, :tan, :sinh, :cosh, :tanh, :asin, :acos,
        :asinh, :acosh, :atanh, :sec, :csc, :cot, :asec, :acsc, :acot, :sech, :csch,
        :coth, :asech, :acsch, :acoth, :log, :log2, :log10, :log1p, :exp, :exp2, :exp10, 
        :expm1, :frexp, :exponent,
    )
    @eval Base.$f(u::AbstractDimensions) = dimensionless(u)
    @eval Base.$f(q::AbstractQuantity) = with_ubase($f, q)
end

#Single-argument functions that only operate on values
for f in (:isfinite, :isinf, :isnan, :isreal, :isempty, :one)
    @eval Base.$f(q::AbstractQuantity) = $f(ustrip(q))
end

#Single-argument functions with dimensionless output that only work with scalar units (no offset)
for f in (:iszero, :angle, :signbit, :sign)
    @eval Base.$f(q::AbstractQuantity) = (assert_scalar(unit(q)); $f(ustrip_base(q)))
end

