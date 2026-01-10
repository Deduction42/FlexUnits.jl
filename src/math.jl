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

#Mirror dimensions on two argument operations produces "equaldims
Base.:*(d1::AbstractDimensions, d2::MirrorDims) = d1
Base.:*(d1::MirrorDims, d2::AbstractDimensions) = d2
Base.:*(d1::MirrorDims, d2::MirrorDims) = d1
Base.:/(d1::AbstractDimensions, d2::MirrorDims) = d1
Base.:/(d1::MirrorDims, d2::AbstractDimensions) = inv(d2)
Base.:/(d1::MirrorDims, d2::MirrorDims) = d1

#Mirror dimensions on single argument functions
for op in (:inv, :sqrt, :cbrt, :abs2, :adjoint)
    @eval Base.$op(d::MirrorDims) = d
end
Base.:^(d::MirrorDims, p::Real) = d


#============================================================================================================================
Static dimension ops
============================================================================================================================#
Base.:(==)(d1::AbstractDimensions, d2::StaticDims) = (d1 == dimval(d2))
Base.:(==)(d1::StaticDims, d2::AbstractDimensions) = (dimval(d1) == d2)

equaldims(arg1::StaticDims, arg2::AbstractDimensions) = (dimval(arg1) == arg2) ? arg1 : throw(DimensionError((arg1,arg2)))
equaldims(arg1::AbstractDimensions, arg2::StaticDims) = (arg1 == dimval(arg2)) ? arg2 : throw(DimensionError((arg1,arg2)))
equaldims(arg1::StaticDims, arg2::StaticDims) = (dimval(arg1) == dimval(arg2))  ? arg1 : throw(DimensionError((arg1,arg2)))

Base.:*(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1*D2}()
Base.:/(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1/D2}()
Base.inv(arg::StaticDims{D}) where D = StaticDims{inv(D)}()
Base.:^(d::StaticDims{D}, p::Real) where D = StaticDims{D^p}()
Base.sqrt(d::StaticDims{D}) where D = StaticDims{sqrt(D)}()
Base.cbrt(d::StaticDims{D}) where D = StaticDims{cbrt(D)}()
Base.abs2(d::StaticDims{D}) where D = StaticDims{abs2(D)}()
Base.adjoint(d::StaticDims{D}) where D = StaticDims{adjoint(D)}()

@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{0}) = StaticDims{dimtype(d)()}()
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{1}) = d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{2}) = d*d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{3}) = d*d*d
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-1}) = inv(d) 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-2}) = inv(d*d)

#=============================================================================================
 Mathematical operations on abstract units and transforms (mostly for parsing)
=============================================================================================#
const NON_SCALAR_ERROR = ArgumentError("Operation only allowed on scalar transforms")

#Base.:+(t::AffineTransform, x::Real) = AffineTransform(offset = t.offset + x, scale = t.scale)
#Base.:+(x::Real, t::AffineTransform) = t + x 
#Base.:-(t::AffineTransform, x::Real) = AffineTransform(offset = t.offset - x, scale = t.scale)

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

Base.:*(x::MathUnion, u::AbstractUnitLike) = quantity(x, u)
Base.:/(x::MathUnion, u::AbstractUnitLike) = quantity(x, inv(u))
Base.:*(x::Missing, u::AbstractUnitLike) = x
Base.:/(x::Missing, u::AbstractUnitLike) = x

Base.:*(q::QuantUnion, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)*u)
Base.:*(u::AbstractUnitLike, q::QuantUnion) = q*u
Base.:/(q::QuantUnion, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)/u)
Base.:/(u::AbstractUnitLike, q::QuantUnion) = quantity(inv(ustrip(q)), u/unit(q))

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

Base.sqrt(u::AbstractUnits) = u^inv(2)
Base.cbrt(u::AbstractUnits) = u^inv(3)
Base.adjoint(u::AbstractUnits) = u

#Equality does not compare symbols
Base.:(==)(u1::AbstractUnits, u2::AbstractUnits) = (todims(u1) == todims(u2)) & (dimension(u1) == dimension(u2))


#=============================================================================================
 Mathematical operations on quantities
=============================================================================================#
"""
    with_ubase(f, args::QuantUnion...)

Converts all arguments to base units, and applies `f` to values and dimensions
returns the quanity. Useful for defining new functions for quantities 
    
Thus, in order to support `f` for quantities, simply define 
```
f(dims::AbstractDimensions...)
f(args::QuantUnion...) = with_ubase(f, args...)
```
"""
function with_ubase(f, args::QuantUnion...)
    baseargs = map(ubase, args)
    basevals = map(ustrip, baseargs)
    basedims = map(unit, baseargs)
    scaleval = f(basevals...)
    return quantity(scaleval, f(basedims...))
end

function Base.:(==)(q1::QuantUnion, q2::QuantUnion)
    qb1 = ubase(q1)
    qb2 = ubase(q2)
    return (ustrip(qb1) == ustrip(qb2)) && (unit(qb1) == unit(qb2))
end

function Base.:(≈)(q1::QuantUnion, q2::QuantUnion)
    qb1 = ubase(q1)
    qb2 = ubase(q2)
    return (ustrip(qb1) ≈ ustrip(qb2)) && (unit(qb1) == unit(qb2))
end
Base.:(≈)(q1::QuantUnion, q2::T) where T<:NumUnion = (convert(T, q1) ≈ q2)
Base.:(≈)(q1::T, q2::QuantUnion) where T<:NumUnion = (convert(T, q2) ≈ q1)
Base.:(≈)(q1::Missing, q2::QuantUnion) = missing 
Base.:(≈)(q2::QuantUnion, q1::Missing) = missing 

Base.:+(q::QuantUnion, x::NumUnion) = dimensionless(q) + x 
Base.:+(x::NumUnion, q::QuantUnion) = dimensionless(q) + x
Base.:+(q1::QuantUnion, q2::QuantUnion) = with_ubase(+, q1, q2)
Base.:+(q1::QuantUnion, qN::QuantUnion...) = with_ubase(+, q1, qN...)

Base.:-(q::QuantUnion, x::NumUnion) = dimensionless(q) - x 
Base.:-(x::NumUnion, q::QuantUnion) = x - dimensionless(q)
Base.:-(q1::QuantUnion, q2::QuantUnion) = with_ubase(-, q1, q2)
Base.:-(q1::QuantUnion) = with_ubase(-, q1)

Base.:*(q0::QuantUnion, x::NumUnion) = (q = ubase(q0); quantity(ustrip(q)*x, unit(q)))
Base.:*(x::NumUnion, q0::QuantUnion) = (q = ubase(q0); quantity(ustrip(q)*x, unit(q)))
Base.:*(q1::QuantUnion, q2::QuantUnion) = with_ubase(*, q1, q2)
Base.:*(q1::QuantUnion, qN::QuantUnion...) = with_ubase(*, q1, qN...)

Base.:/(q0::QuantUnion, x::NumUnion) = (q = ubase(q0); quantity(ustrip(q)/x, unit(q)))
Base.:/(x::NumUnion, q0::QuantUnion) = (q = ubase(q0); quantity(x/ustrip(q), inv(unit(q))))
Base.:/(q1::QuantUnion, q2::QuantUnion) = with_ubase(/, q1, q2)
Base.:inv(q::QuantUnion) = with_ubase(inv, q)
Base.adjoint(q::QuantUnion) = with_ubase(adjoint, q)

#Operators on explicitly missing values simply return missing
for op in (:+,:-,:*,:/)
    @eval Base.$op(q::QuantUnion, x::Missing) = x
    @eval Base.$op(x::Missing, q::QuantUnion) = x
end

Base.:^(q::QuantUnion, p::Real) = with_ubase(Base.Fix2(^, p), q)
Base.:^(q::QuantUnion, p::Integer) = with_ubase(Base.Fix2(^, p), q)
Base.:^(q::QuantUnion, p::Rational) = with_ubase(Base.Fix2(^, p), q)
#Base.:^(q::QuantUnion, p::Real) = with_ubase(Base.Fix2(^, p), q)
#Base.:^(q::MathUnion, p::Quantity)  = q^dimensionless(p)
#Base.:^(q::QuantUnion, p::Quantity) = with_ubase(Base.Fix2(^, dimensionless(p)), q)

@inline Base.literal_pow(::typeof(^), q::QuantUnion, ::Val{p}) where {p} = with_ubase(x->Base.literal_pow(^, x, Val(dimensionless(p))), q)

Base.sqrt(q::QuantUnion) = with_ubase(sqrt, q)
Base.cbrt(q::QuantUnion) = with_ubase(cbrt, q)
Base.abs2(q::QuantUnion) = with_ubase(abs2, q)
Base.max(q1::QuantUnion, q2::QuantUnion) = with_ubase(max, q1, q2)
Base.min(q1::QuantUnion, q2::QuantUnion) = with_ubase(min, q1, q2)
Base.zero(::Type{D}) where D<:AbstractDimensions = D()
Base.zero(::Type{U}) where {D,T,U<:Units{D,T}} = U(dims=D(), todims=T())

#Common functions for initializers
Base.one(::Type{<:QuantUnion{T}}) where T = one(T) #unitless
Base.oneunit(::Type{<:QuantUnion{T,D}}) where {T,D<:StaticDims} = quantity(one(T), D()) #One with units
Base.oneunit(::Type{<:QuantUnion{T,D}}) where {T,D<:AbstractDimensions} = throw(ArgumentError("Cannot inver value of a dynamic dimension from its type"))
Base.oneunit(::Type{<:QuantUnion{T,D}}) where {T,D<:Units} = throw(ArgumentError("Cannot inver value of a dynamic dimension from its type"))

for f in (:zero, :typemin, :typemax)
    @eval Base.$f(::Type{<:QuantUnion{T, D}}) where {T, D<:AbstractDimensions} = quantity($f(T), MirrorDims(D))
    @eval Base.$f(::Type{<:QuantUnion{T, D}}) where {T, D<:StaticDims} = quantity($f(T), D())
    @eval Base.$f(::Type{<:QuantUnion{T, <:MirrorDims{D}}}) where {T, D<:AbstractDimensions} = $f(quant_type(T){T,D})
    @eval Base.$f(::Type{<:QuantUnion{T, <:MirrorUnion{D}}}) where {T, D<:AbstractDimensions} = $f(quant_type(T){T,D})
    @eval Base.$f(::Type{<:QuantUnion{T, <:AbstractUnits{D}}}) where {T, D<:AbstractDimensions} = $f(quant_type(T){T,D})
end

#Comparison functions (returns a bool)
for f in (:<, :<=, :isless)
    @eval function Base.$f(q1::QuantUnion, q2::QuantUnion)
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
    @eval Base.$f(q::QuantUnion) = with_ubase($f, q)
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
    @eval Base.$f(q::QuantUnion) = with_ubase($f, q)
end

#Single-argument functions that only operate on values
for f in (:isfinite, :isinf, :isnan, :isreal, :isempty, :one)
    @eval Base.$f(q::QuantUnion) = $f(ustrip(q))
end

#Single-argument functions with dimensionless output that only work with scalar units (no offset)
for f in (:iszero, :angle, :signbit, :sign)
    @eval Base.$f(q::QuantUnion) = (assert_scalar(unit(q)); $f(ustrip_base(q)))
end

