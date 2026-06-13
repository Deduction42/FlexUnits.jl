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

unknown(::Type{D}) where {T, D<:AbstractDimensions{T}} = D(typemax(T))
isunknown(d::D, fn::Symbol) where {T, D<:AbstractDimensions{T}} = (getproperty(d, fn) == typemax(T))
isunknown(d::D) where D<:AbstractDimensions = isunknown(d, dimension_names(D)[end]) #Only need to check one dimension for speed
isknown(d::D) where D<:AbstractDimensions = !isunknown(d, dimension_names(D)[end])
isunknown(d::StaticDims{D}) where D = isunknown(D)
isknown(d::StaticDims{D}) where D = isknown(D)

#Checks equality of dimensions, returns first non-mirrored dimension
@inline equaldims(arg1::AbstractDimensions) = arg1
equaldims(arg1::AbstractDimensions, arg2::AbstractDimensions) = equaldims(promote(arg1, arg2)...)
function equaldims(d1::D, d2::D) where D<:AbstractDimensions
    if (d1 == d2)
        return d1
    elseif isunknown(d2)
        return d1
    elseif isunknown(d1)
        return d2
    else
        throw(DimensionError(d1,d2))
    end
end

raw_mul(d1::AbstractDimensions, d2::AbstractDimensions) = map_dimensions(+, d1, d2)
raw_div(d1::AbstractDimensions, d2::AbstractDimensions) = map_dimensions(-, d1, d2)
raw_pow(d::AbstractDimensions, p::Integer) = map_dimensions(Base.Fix1(*, p), d)
raw_pow(d::AbstractDimensions{R}, p::Real) where R = map_dimensions(Base.Fix1(*, convert(R, p)), d)
raw_inv(d::AbstractDimensions) = map_dimensions(-, d)

#=============================================================================================
 Mathematical operations on dimensions
=============================================================================================#
function Base.:*(d1::AbstractDimensions, d2::AbstractDimensions) 
    if isunknown(d1)
        return d1 
    elseif isunknown(d2)
        return d2 
    else
        return raw_mul(d1, d2)
    end
end

function Base.:/(d1::AbstractDimensions, d2::AbstractDimensions)
    if isunknown(d1)
        return d1 
    elseif isunknown(d2)
        return d2 
    else
        return raw_div(d1, d2)
    end
end

function Base.:^(d::AbstractDimensions, p::Real) 
    if iszero(p)
        return typeof(d)()
    elseif isknown(d) 
        return raw_pow(d, p) 
    else
        return d 
    end
end

Base.inv(d::AbstractDimensions) = isknown(d) ? raw_inv(d) : d


@inline Base.:+(arg1::AbstractDimLike) = arg1
@inline Base.:+(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
@inline Base.:-(arg1::AbstractDimLike) = arg1
@inline Base.:-(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
Base.min(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
Base.max(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
Base.rem(arg1::AbstractDimLike, arg2::AbstractDimLike) = equaldims(arg1, arg2)
Base.sqrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 2))
Base.cbrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 3))
Base.abs2(d::AbstractDimensions) = d^2
Base.adjoint(d::AbstractDimensions) = d

@inline zero_pow(d::D) where D<:AbstractDimensions = D()
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{0}) = zero_pow(d)
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{1}) = d 
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{2}) = isknown(d) ? raw_mul(d, d) : d
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{3}) = isknown(d) ? raw_mul(d, raw_mul(d, d)) : d
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{-1}) = inv(d) 
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{-2}) = isknown(d) ? raw_inv(raw_mul(d, d)) : d


#NoDims shortcuts
@inline equaldims(d1::AbstractDimensions, d2::NoDims) = assert_dimensionless(d1)
@inline equaldims(d1::NoDims, a2::AbstractDimensions) = assert_dimensionless(d2)
@inline equaldims(d1::NoDims, d2::NoDims) = d1

Base.:(==)(d1::NoDims, d2::AbstractDimLike) = isdimensionless(d2)
Base.:(==)(d1::AbstractDimLike, d2::NoDims) = isdimensionless(d1)
Base.:(==)(d1::NoDims, d2::NoDims) = true

Base.:*(d1::AbstractDimensions, d2::NoDims) = d1
Base.:*(d1::NoDims, d2::AbstractDimensions) = d2
Base.:*(d1::NoDims, d2::NoDims) = d1

Base.:/(d1::AbstractDimensions, d2::NoDims) = d1
Base.:/(d1::NoDims, d2::AbstractDimensions) = inv(d2)
Base.:/(d1::NoDims, d2::NoDims) = d1

for op in (:inv, :sqrt, :cbrt, :abs2, :adjoint)
    @eval Base.$op(d::NoDims) = d
end
Base.:^(d::NoDims, p::Real) = d


#============================================================================================================================
Static dimension ops
============================================================================================================================#
Base.:(==)(d1::AbstractDimensions, d2::StaticDims) = (d1 == dimval(d2))
Base.:(==)(d1::StaticDims, d2::AbstractDimensions) = (dimval(d1) == d2)
Base.:(==)(d1::NoDims, d2::StaticDims) = (d1 == dimval(d2))
Base.:(==)(d1::StaticDims, d2::NoDims) = (dimval(d1) == d2)

equaldims(arg1::StaticDims, arg2::AbstractDimensions) = (dimval(arg1) == arg2) || isunknown(arg2) ? arg1 : throw(DimensionError(arg1,arg2))
equaldims(arg1::AbstractDimensions, arg2::StaticDims) = (arg1 == dimval(arg2)) || isunknown(arg1) ? arg2 : throw(DimensionError(arg1,arg2))
equaldims(arg1::StaticDims, arg2::StaticDims) = (dimval(arg1) == dimval(arg2)) ? arg1 : throw(DimensionError(arg1,arg2))

Base.:*(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1*D2}()
Base.:/(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1/D2}()
Base.inv(arg::StaticDims{D}) where D = StaticDims{inv(D)}()
Base.:^(d::StaticDims{D}, p::Real) where D = StaticDims{D^p}()
Base.sqrt(d::StaticDims{D}) where D = StaticDims{sqrt(D)}()
Base.cbrt(d::StaticDims{D}) where D = StaticDims{cbrt(D)}()
Base.abs2(d::StaticDims{D}) where D = StaticDims{abs2(D)}()
Base.adjoint(d::StaticDims{D}) where D = StaticDims{adjoint(D)}()

@inline zero_pow(d::StaticDims{D}) where D = StaticDims{zero_pow(D)}()
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{0}) = zero_pow(d)
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{1}) = d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{2}) = d*d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{3}) = d*d*d
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-1}) = inv(d) 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-2}) = inv(d*d)

#Shortcuts to avoid generic fallbacks/promotion
Base.:*(d1::StaticDims, d2::NoDims) = d1
Base.:*(d1::NoDims, d2::StaticDims) = d2
Base.:/(d1::StaticDims, d2::NoDims) = d1
Base.:/(d1::NoDims, d2::StaticDims) = inv(d2)

#=============================================================================================
 Mathematical operations on abstract units and transforms (mostly for parsing)
=============================================================================================#
const NON_SCALAR_ERROR = ArgumentError("Operation only allowed on scalar transforms")

function Base.:*(t1::AffineTransform, t2::AffineTransform) 
    is_scalar(t1) & is_scalar(t2) || throw(NON_SCALAR_ERROR)
    return AffineTransform(scale = t1.scale*t2.scale, offset = 0) 
end
Base.:*(t::AffineTransform{T}, x::Real) where T = is_scalar(t) ? AffineTransform{T}(scale=t.scale*x, offset=0) : throw(NON_SCALAR_ERROR)
Base.:*(t::NoTransform, x::Real) = AffineTransform(scale=x, offset=0)

function Base.:/(t1::AffineTransform, t2::AffineTransform) 
    is_scalar(t1) & is_scalar(t2) || throw(NON_SCALAR_ERROR)
    return AffineTransform(scale = t1.scale/t2.scale, offset = 0) 
end
Base.:/(t::AffineTransform{T}, x::Real) where T = is_scalar(t) ? AffineTransform{T}(scale=t.scale/x, offset=0) : throw(NON_SCALAR_ERROR)
Base.:/(t1::NoTransform, x::Real) = AffineTransform(scale=inv(x), offset=0)

function Base.:^(t::AffineTransform{T}, p::Real) where T
    is_scalar(t) || throw(NON_SCALAR_ERROR)
    return AffineTransform{T}(scale = t.scale^p, offset = 0) 
end
Base.:^(t1::NoTransform, p::Real) = t1

Base.:*(x::MathUnion, u::AbstractUnitLike) = ubase(x, u)
Base.:/(x::MathUnion, u::AbstractUnitLike) = ubase(x, inv(u))
Base.:*(x::Missing, u::AbstractUnitLike) = x
Base.:/(x::Missing, u::AbstractUnitLike) = x

Base.:*(q::QuantUnion, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)*u)
Base.:*(u::AbstractUnitLike, q::QuantUnion) = q*u
Base.:/(q::QuantUnion, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)/u)
Base.:/(u::AbstractUnitLike, q::QuantUnion) = quantity(inv(ustrip(q)), u/unit(q))
Base.:*(m::AbstractArray{<:QuantUnion}, u::AbstractUnitLike) = broadcast(*, m, u)
Base.:*(u::AbstractUnitLike, m::AbstractArray{<:QuantUnion}) = broadcast(*, m, u)

function Base.:*(u1::U, u2::U) where U <: AbstractUnits
    return constructorof(U)(scalar_dimension(u1)*scalar_dimension(u2), tobase(u1)*tobase(u2))
end

function Base.:/(u1::U, u2::U) where U <: AbstractUnits
    return constructorof(U)(scalar_dimension(u1)/scalar_dimension(u2), tobase(u1)/tobase(u2))
end

function Base.:inv(u::U) where U <: AbstractUnits
    return constructorof(U)(inv(scalar_dimension(u)), inv(tobase(u)))
end

function Base.:^(u::U, p::Real) where U <:AbstractUnits
    return constructorof(U)(scalar_dimension(u)^p, tobase(u)^p)
end

Base.:*(u1::AbstractUnitLike, u2::AbstractUnitLike) = *(promote(u1,u2)...)
Base.:/(u1::AbstractUnitLike, u2::AbstractUnitLike) = /(promote(u1,u2)...)

Base.sqrt(u::AbstractUnits) = u^inv(2)
Base.cbrt(u::AbstractUnits) = u^inv(3)
Base.adjoint(u::AbstractUnits) = u

#Equality does not compare symbols
Base.:(==)(u1::AbstractUnits, u2::AbstractUnits) = (tobase(u1) == tobase(u2)) & (dimension(u1) == dimension(u2))


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

Base.:+(q::QuantUnion, x::NumUnion) = scalar(q) + x 
Base.:+(x::NumUnion, q::QuantUnion) = scalar(q) + x
Base.:+(q1::QuantUnion, q2::QuantUnion) = with_ubase(+, q1, q2)
Base.:+(q1::QuantUnion, qN::QuantUnion...) = with_ubase(+, q1, qN...)

Base.:-(q::QuantUnion, x::NumUnion) = scalar(q) - x 
Base.:-(x::NumUnion, q::QuantUnion) = x - scalar(q)
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

Base.muladd(x::QuantUnion, y::QuantUnion, z::QuantUnion) = muladd(dstrip(x), dstrip(y), dstrip(z)) * equaldims(dimension(x)*dimension(y), dimension(z))

#Operators on explicitly missing values simply return missing
for op in (:+,:-,:*,:/)
    @eval Base.$op(q::QuantUnion, x::Missing) = x
    @eval Base.$op(x::Missing, q::QuantUnion) = x
end

Base.:^(q::QuantUnion, p::Real) = with_ubase(Base.Fix2(^, p), q)
Base.:^(q::QuantUnion, p::Integer) = with_ubase(Base.Fix2(^, p), q)
Base.:^(q::QuantUnion, p::Rational) = with_ubase(Base.Fix2(^, p), q)

@inline Base.literal_pow(::typeof(^), q::QuantUnion, ::Val{p}) where {p} = with_ubase(x->Base.literal_pow(^, x, Val(dimensionless(p))), q)

Base.sqrt(q::QuantUnion) = with_ubase(sqrt, q)
Base.cbrt(q::QuantUnion) = with_ubase(cbrt, q)
Base.abs2(q::QuantUnion) = with_ubase(abs2, q)
Base.max(q1::QuantUnion, q2::QuantUnion) = with_ubase(max, q1, q2)
Base.min(q1::QuantUnion, q2::QuantUnion) = with_ubase(min, q1, q2)
Base.rem(q1::QuantUnion, q2::QuantUnion) = with_ubase(rem, q1, q2)
Base.zero(::Type{D}) where D<:AbstractDimensions = D()
Base.zero(::Type{U}) where {D,T,U<:Units{D,T}} = U(dims=D(), tobase=T())

#Common functions for initializers
Base.one(::Type{<:QuantUnion{T}}) where T = one(T) #unitless
Base.rtoldefault(::Type{<:QuantUnion{T}}) where T = Base.rtoldefault(T)

for f in (:zero, :typemin, :typemax, :oneunit, :eps)
    @eval Base.$f(::Type{<:QuantUnion{T, D}}) where {T, D<:AbstractDimensions} = quantity($f(T), unknown(D))
    @eval Base.$f(::Type{<:QuantUnion{T, D}}) where {T, D<:StaticDims} = quantity($f(T), D())
    @eval Base.$f(::Type{<:QuantUnion{T, U}}) where {T, U<:AbstractUnits} = throw(ArgumentError("This operation only supports dimensional quantities"))
end

#Comparison functions (returns a bool)
for f in (:<, :<=, :isless)
    @eval function Base.$f(q1::QuantUnion, q2::QuantUnion)
        (b1, b2) = (ubase(q1), ubase(q2))
        equaldims(unit(b1), unit(b2))
        return $f(ustrip(b1), ustrip(b2))
    end
    @eval Base.$f(q::QuantUnion, x::Real) = $f(dimensionless(q), x)
    @eval Base.$f(x::Real, q::QuantUnion) = $f(x, dimensionless(q))
end

#Functions that return the same unit
for f in (
        :float, :abs, :real, :imag, :conj, :significand, :zero, :eps, :oneunit, :nextfloat,
        :typemax, :typemin, :transpose
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
    @eval Base.$f(q::QuantUnion) = $f(scalar(q))
end

#Single-argument functions that only operate on values
for f in (:isfinite, :isinf, :isnan, :isreal, :isempty, :one)
    @eval Base.$f(q::QuantUnion) = $f(ustrip(q))
end

#Single-argument functions with dimensionless output that only work with scalar units (no offset)
for f in (:iszero, :angle, :signbit, :sign)
    @eval Base.$f(q::QuantUnion) = (assert_scalar(unit(q)); $f(ustrip_base(q)))
end


#================================================================================================================
Logarithmic Quantities, Units and Transforms
================================================================================================================#
"""
    with_logubase(f, args::LogQuant...)

Converts all arguments to log base units, and applies `f` to values and dimensions
returns the quanity. Useful for defining new functions for quantities 
"""
function with_logubase(fv, fd, args::LogQuant...)
    baseargs = map(logubase, args)
    basevals = map(ustrip, baseargs)
    basedims = map(unit, baseargs)
    return logquant(fv(basevals...), fd(basedims...))
end

function Base.:(==)(q1::LogQuant, q2::LogQuant)
    qb1 = logubase(q1)
    qb2 = logubase(q2)
    return (ustrip(qb1) == ustrip(qb2)) && (unit(qb1) == unit(qb2))
end

function Base.:(≈)(q1::LogQuant, q2::LogQuant)
    qb1 = logubase(q1)
    qb2 = logubase(q2)
    return (ustrip(qb1) ≈ ustrip(qb2)) && (unit(qb1) == unit(qb2))
end

Base.:(≈)(q1::Missing, q2::LogQuant) = missing 
Base.:(≈)(q2::LogQuant, q1::Missing) = missing 

#Logarithms/Exponentials of transforms 
Base.log(q::Quantity) = logquant(q)
Base.log(t::ExpAffTransform) = AffineTransform(scale=uscale(t), offset=uoffset(t))

Base.exp(q::LogQuant) = quantity(q)
Base.exp(t::AffineTransform) = ExpAffTransform(scale=uscale(t), offset=uoffset(t))
Base.exp(t::NoTransform) = ExpAffTransform(scale=uscale(t), offset=uoffset(t))

#Mutliplication with Unit{D, ExpAffTransform} produces a LogQuant 
function Base.:*(x::NumUnion, u::Units{<:AbstractDimLike,<:ExpAffTransform})
    logtrans = log(tobase(u))
    return logquant(logtrans(x), dimension(u))
end
Base.:*(x::Quantity, u::Units{<:AbstractDimLike,<:ExpAffTransform}) = scalar(x)*u

#Logarithmic quantity algebra
Base.:+(q::LogQuant) = with_logubase(+, *, q)
Base.:+(q1::LogQuant, q2::LogQuant) = with_logubase(+, *, q1, q2)
Base.:+(x::Real, q0::LogQuant) = (q = logubase(q0); logquant(ustrip(q) + x, unit(q)))
Base.:+(q0::LogQuant, x::Real) = (q = logubase(q0); logquant(ustrip(q) + x, unit(q)))
Base.:+(q1::Quantity, q2::LogQuant) = scalar(q1) + q2
Base.:+(q1::LogQuant, q2::Quantity) = q1 + scalar(q2)

Base.:-(q::LogQuant) = with_logubase(-, inv, q)
Base.:-(q1::LogQuant, q2::LogQuant) = with_logubase(-, /, q1, q2)
Base.:-(x::Real, q0::LogQuant) = (q = logubase(q0); logquant(x - ustrip(q), inv(unit(q))))
Base.:-(q0::LogQuant, x::Real) = (q = logubase(q0); logquant(ustrip(q) - x, unit(q)))
Base.:-(q1::Quantity, q2::LogQuant) = scalar(q1) - q2
Base.:-(q1::LogQuant, q2::Quantity) = q1 - scalar(q2)

Base.:*(q0::LogQuant, x::Real) = (q = logubase(q0); logquant(ustrip(q)*x, unit(q)^x))
Base.:*(x::Real, q0::LogQuant) = (q = logubase(q0); logquant(ustrip(q)*x, unit(q)^x))
Base.:/(q0::LogQuant, x::Real) = (q = logubase(q0); logquant(ustrip(q)/x, unit(q)^inv(x)))
Base.:*(q1::Quantity, q2::LogQuant) = scalar(q1)*q2
Base.:*(q1::LogQuant, q2::Quantity) = q1*scalar(q2)
Base.:/(q1::LogQuant, q2::Quantity) = q1/scalar(q2)

#Addition/subtraction for linear transformations
⊕(q1::LogQuant, q2::LogQuant) = log(ubase(q1) + ubase(q2))
⊕(q1::NumUnion, q2::LogQuant) = log(exp(q1) + ubase(q2))
⊕(q1::LogQuant, q2::NumUnion) = log(ubase(q1) + exp(q2))

⊖(q1::LogQuant, q2::LogQuant) = log(ubase(q1) - ubase(q2))
⊖(q1::NumUnion, q2::LogQuant) = log(exp(q1) - ubase(q2))
⊖(q1::LogQuant, q2::NumUnion) = log(ubase(q1) - exp(q2))

#Comparison functions (returns a bool)
for f in (:<, :<=, :isless)
    @eval function Base.$f(q1::LogQuant, q2::LogQuant)
        (b1, b2) = (logubase(q1), logubase(q2))
        equaldims(unit(b1), unit(b2))
        return $f(ustrip(b1), ustrip(b2))
    end
end

#Unitless remainders 
function Base.rem(lq1::LogQuant{<:Number,D1}, lq2::LogQuant{<:Number,D2}) where {D1<:AbstractDimLike, D2<:AbstractDimLike}
    d = equaldims(dimension(lq1), dimension(lq2))
    isdimensionless(d) || error("Operation only supported for dimensionless LogQuant")
    return LogQuant(rem(dstrip(lq1), dstrip(lq2)), d)
end
Base.rem(lq1::LogQuant, lq2::LogQuant) = rem(logubase(lq1), logubase(lq2))