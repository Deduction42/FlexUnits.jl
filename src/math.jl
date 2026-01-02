#Arithmetic Types in Base that support +-*/ (use instead of Any to avoid ambiguity)
const ARITHMETICS = Union{Number, AbstractArray, Missing}

#=============================================================================================
 Dimension fundamentals
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

#Checks equality of dimensions, returns first non-mirrored dimension
@inline equaldims(arg1::AbstractDimensions) = arg1
@inline equaldims(arg1::AbstractDimensions, arg2::MirrorDims) = arg1
@inline equaldims(arg1::MirrorDims, arg2::AbstractDimensions) = arg2
@inline equaldims(arg1::MirrorDims, arg2::MirrorDims) = arg1

function equaldims(arg1::AbstractDimensions, arg2::AbstractDimensions)
    (d1, d2) = promote(arg1, arg2)
    return (d1 == d2) ? d1 : throw(DimensionError((arg1,arg2)))
end

#=============================================================================================
 Mathematical operations on dimensions
=============================================================================================#
@inline Base.:+(arg1::AbstractDimLike) = arg1
@inline Base.:+(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
@inline Base.:-(arg1::AbstractDimLike) = arg1
@inline Base.:-(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
Base.min(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
Base.max(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
Base.:*(arg1::AbstractDimensions, arg2::AbstractDimensions) = map_dimensions(+, arg1, arg2)
Base.:/(arg1::AbstractDimensions, arg2::AbstractDimensions) = map_dimensions(-, arg1, arg2)
Base.inv(arg::AbstractDimensions) = map_dimensions(-, arg)
Base.:^(d::AbstractDimensions, p::Integer) = map_dimensions(Base.Fix1(*, p), d)
Base.:^(d::AbstractDimensions{R}, p::Real) where {R} = map_dimensions(Base.Fix1(*, R(dimensionless(p))), d)
Base.sqrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 2))
Base.cbrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 3))
Base.abs2(d::AbstractDimensions) = d^2
Base.adjoint(d::AbstractDimensions) = d

@inline Base.literal_pow(::typeof(^), d::D, ::Val{0}) where {D <: AbstractDimensions} = D()
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{1}) = d 
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{2}) = d*d 
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{3}) = d*d*d
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{-1}) = inv(d) 
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{-2}) = inv(d*d)

#Mirror dimensions on two argument operations produces "equaldims"
for op in (:*, :/)
    @eval Base.$op(d1::AbstractDimensions, d2::MirrorDims) = d1
    @eval Base.$op(d1::MirrorDims, d2::AbstractDimensions) = d2
    @eval Base.$op(d1::MirrorDims, d2::MirrorDims) = d1
end

#Mirror dimensions on single argument functions
for op in (:inv, :sqrt, :cbrt, :abs2, :adjoint)
    @eval Base.$op(d::MirrorDims) = d
end
Base.:^(d::MirrorDims, p::Real) = d

#=============================================================================================
 Mathematical operations on abstract units and transforms (mostly for parsing)
=============================================================================================#
const NON_SCALAR_ERROR = ArgumentError("Operation only allowed on scalar transforms")

Base.:+(t::AffineTransform, x::Real) = AffineTransform(offset = t.offset + x, scale = t.scale)
Base.:+(x::Real, t::AffineTransform) = t + x 
Base.:-(t::AffineTransform, x::Real) = AffineTransform(offset = t.offset - x, scale = t.scale)

function Base.:*(t1::AffineTransform, t2::AffineTransform) 
    is_scalar(t1) & is_scalar(t2) || throw(NON_SCALAR_ERROR)
    return AffineTransform(scale = t1.scale*t2.scale, offset = 0) 
end
Base.:*(t::AffineTransform, x::Real) = is_scalar(t) ? AffineTransform(scale=t.scale*x, offset=0) : throw(NON_SCALAR_ERROR)
Base.:*(t::NoTransform, x::Real) = AffineTransform(scale=x, offset=0)

function Base.:/(t1::AffineTransform, t2::AffineTransform) 
    is_scalar(t1) & is_scalar(t2) || throw(NON_SCALAR_ERROR)
    return AffineTransform(scale = t1.scale/t2.scale, offset = 0) 
end
Base.:/(t::AffineTransform, x::Real) = is_scalar(t) ? AffineTransform(scale=t.scale/x, offset=0) : throw(NON_SCALAR_ERROR)
Base.:/(t1::NoTransform, x::Real) = AffineTransform(scale=inv(x), offset=0)

function Base.:^(t::AffineTransform, p::Real) 
    is_scalar(t) || throw(NON_SCALAR_ERROR)
    return AffineTransform(scale = t.scale^p, offset = 0) 
end
Base.:^(t1::NoTransform, p::Real) = t1

Base.:*(x::ARITHMETICS, u::AbstractUnitLike) = Quantity(x, u)
Base.:/(x::ARITHMETICS, u::AbstractUnitLike) = Quantity(x, inv(u))

Base.:*(q::AbstractQuantity, u::AbstractUnitLike) = Quantity(ustrip(q), unit(q)*u)
Base.:*(u::AbstractUnitLike, q::AbstractQuantity) = q*u
Base.:/(q::AbstractQuantity, u::AbstractUnitLike) = Quantity(ustrip(q), unit(q)/u)
Base.:/(u::AbstractUnitLike, q::AbstractQuantity) = Quantity(inv(ustrip(q)), u/unit(q))

function Base.:*(u1::U, u2::U) where U <: AbstractUnits
    return constructorof(U)(scalar_dimension(u1)*scalar_dimension(u2), todims(u1)*todims(u2))
end

function Base.:/(u1::U, u2::U) where U <: AbstractUnits
    return constructorof(U)(scalar_dimension(u1)/scalar_dimension(u2), todims(u1)/todims(u2))
end

function Base.:inv(u::U) where U <: AbstractUnits
    return constructorof(U)(inv(scalar_dimension(u)), inv(todims(u)))
end

function Base.:^(u::U, p::Real) where U <:AbstractUnits
    return constructorof(U)(scalar_dimension(u)^p, todims(u)^p)
end

Base.:*(u1::AbstractUnitLike, u2::AbstractUnitLike) = *(promote(u1,u2)...)
Base.:/(u1::AbstractUnitLike, u2::AbstractUnitLike) = /(promote(u1,u2)...)

Base.sqrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 2))
Base.cbrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 3))
Base.adjoint(u::AbstractUnits) = u

#Equality does not compare symbols
Base.:(==)(u1::AbstractUnits, u2::AbstractUnits) = (todims(u1) == todims(u2)) & (dimension(u1) == dimension(u2))


#=============================================================================================
 Mathematical operations on quantities
=============================================================================================#
"""
    with_ubase(f, args::AbstractQuantity...)

Converts all arguments to base units, and applies `f` to values and dimensions
returns the quanity. Useful for defining new functions for quantities 
    
Thus, in order to support `f` for quantities, simply define 

f(dims::AbstractDimensions...)
f(args::AbstractQuantity...) = with_ubase(f, args...)
"""
function with_ubase(f, args::AbstractQuantity...)
    baseargs = map(ubase, args)
    basevals = map(ustrip, baseargs)
    basedims = map(unit, baseargs)
    scaleval = f(basevals...)
    return Quantity(scaleval, f(basedims...))
end

#=
function with_ubase(f, args::AbstractQuantity{<:Any, Union{U,MirrorDims{P}}}...) where {P, U<:AbstractUnitLike}
    scaleval  = f(map(ustrip_base, args)...)
    resultdim = f(map(dimension, args)...)
    
    return Quantity{typeof(scaleval), Union{dimtype(U),MirrorDims{P}}}(scaleval, resultdim)
end

function with_ubase(f, args::AbstractQuantity{<:Any, MirrorDims{P}}...) where P
    scaleval  = f(map(ustrip_base, args)...)
    return Quantity(scaleval, MirrorDims{P}())
end
=#

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
Base.:(≈)(q1::AbstractQuantity, q2::T) where T<:ARITHMETICS = (convert(T, q1) ≈ q2)
Base.:(≈)(q1::T, q2::AbstractQuantity) where T<:ARITHMETICS = (convert(T, q2) ≈ q1)
Base.:(≈)(q1::Missing, q2::AbstractQuantity) = missing 
Base.:(≈)(q2::AbstractQuantity, q1::Missing) = missing 

Base.:+(q::AbstractQuantity, x::ARITHMETICS) = dimensionless(q) + x 
Base.:+(x::ARITHMETICS, q::AbstractQuantity) = dimensionless(q) + x
Base.:+(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(+, q1, q2)
Base.:+(q1::AbstractQuantity, qN::AbstractQuantity...) = with_ubase(+, q1, qN...)

Base.:-(q::AbstractQuantity, x::ARITHMETICS) = dimensionless(q) - x 
Base.:-(x::ARITHMETICS, q::AbstractQuantity) = x - dimensionless(q)
Base.:-(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(-, q1, q2)
Base.:-(q1::AbstractQuantity) = with_ubase(-, q1)

Base.:*(q0::AbstractQuantity, x::ARITHMETICS) = (q = ubase(q0); Quantity(ustrip(q)*x, unit(q)))
Base.:*(x::ARITHMETICS, q0::AbstractQuantity) = (q = ubase(q0); Quantity(ustrip(q)*x, unit(q)))
Base.:*(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(*, q1, q2)
Base.:*(q1::AbstractQuantity, qN::AbstractQuantity...) = with_ubase(*, q1, qN...)

Base.:/(q0::AbstractQuantity, x::ARITHMETICS) = (q = ubase(q0); Quantity(ustrip(q)/x, unit(q)))
Base.:/(x::ARITHMETICS, q0::AbstractQuantity) = (q = ubase(q0); Quantity(x/ustrip(q), inv(unit(q))))
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
Base.max(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(max, q1, q2)
Base.min(q1::AbstractQuantity, q2::AbstractQuantity) = with_ubase(min, q1, q2)
Base.zero(::Type{D}) where D<:AbstractDimensions = D()
Base.zero(::Type{U}) where {D,T,U<:Units{D,T}} = U(dims=D(), todims=T())

#Common functions for initializers
Base.one(::Type{<:AbstractQuantity{T}}) where T = one(T) #unitless
Base.oneunit(::Type{<:AbstractQuantity{T,D}}) where {T,D} = Quantity(one(T), zero(D)) #unitless with type

for f in (:zero, :typemin, :typemax)
    @eval Base.$f(::Type{<:AbstractQuantity{T, D}}) where {T, D<:AbstractDimensions} = Quantity($f(T), MirrorDims(D))
    @eval Base.$f(::Type{<:AbstractQuantity{T, D}}) where {T, D<:StaticDims} = Quantity($f(T), D())
    @eval Base.$f(::Type{<:AbstractQuantity{T, <:MirrorDims{D}}}) where {T, D<:AbstractDimensions} = $f(Quantity{T,D})
    @eval Base.$f(::Type{<:AbstractQuantity{T, <:MirrorUnion{D}}}) where {T, D<:AbstractDimensions} = $f(Quantity{T,D})
    @eval Base.$f(::Type{<:AbstractQuantity{T, <:AbstractUnits{D}}}) where {T, D<:AbstractDimensions} = $f(Quantity{T,D})
end

#Comparison functions (returns a bool)
for f in (:<, :<=, :isless)
    @eval function Base.$f(q1::AbstractQuantity, q2::AbstractQuantity)
        (b1, b2) = (ubase(q1), ubase(q2))
        equaldims(unit(b1), unit(b2))
        return $f(ustrip(b1), ustrip(b2))
    end
end

#Functions that return the same unit
for f in (
        :float, :abs, :real, :imag, :conj, :significand, :zero, :oneunit, :typemax, :typemin, :transpose
    )
    @eval Base.$f(u::AbstractDimLike) = u
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
    @eval Base.$f(u::AbstractDimLike) = dimensionless(u)
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

