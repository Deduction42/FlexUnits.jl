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

Base.:+(arg1::AbstractDimensions, args::AbstractDimensions...) = firstequal(arg1, args...)
Base.:-(arg1::AbstractDimensions, args::AbstractDimensions...) = firstequal(arg1, args...)
Base.:*(arg1::AbstractDimensions, args::AbstractDimensions...) = map_dimensions(+, arg1, args...)
Base.:/(arg1::AbstractDimensions, args::AbstractDimensions...) = map_dimensions(-, arg1, args...)
Base.inv(arg::AbstractDimensions) = map_dimensions(-, arg)
Base.:^(d::AbstractDimensions, p::Integer) = map_dimensions(Base.Fix1(*, p), d)
Base.:^(d::AbstractDimensions{R}, p::Real) where {R} = map_dimensions(Base.Fix1(*, R(dimensionless(p))), d)
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{p}) where {p} = map_dimensions(Base.Fix1(*, p), d)
Base.sqrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 2))
Base.cbrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 3))
Base.abs2(d::AbstractDimensions) = d^2

function Base.:(==)(d1::D1, d2::D2) where {D1<:AbstractDimensions, D2<:AbstractDimensions} 
    static_fieldnames(D1) == static_fieldnames(D2)
    return all(fn->d1[fn]==d2[fn], static_fieldnames(D1))
end

#=============================================================================================
 Mathematical operations on abstract units (mostly for parsing)
=============================================================================================#
Base.:*(n::Number, u::AbstractUnitLike) = quantity(n, u)
Base.:/(n::Number, u::AbstractUnitLike) = quantity(n, inv(u))

Base.:*(q::UnionQuantity, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)*u)
Base.:/(q::UnionQuantity, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)/u)

function Base.:*(u1::U, uN::U...) where U <: AbstractUnits
    return constructorof(U)(scale=*(map(uscale, (u1, uN...))...), dims=*(map(scalar_dimension, (u1, uN...))...))
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

Base.:*(args::AbstractUnitLike...) = *(promote(args...)...)
Base.:/(u1::AbstractUnitLike, u2::AbstractUnitLike) = /(promote(u1,u2)...)

Base.sqrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 2))
Base.cbrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 3))

#Equality does not compare symbols
Base.:(==)(u1::AbstractAffineUnits, u2::AbstractAffineUnits) = (uscale(u1) == uscale(u2)) & (uoffset(u1) == uoffset(u2)) & (dimension(u1) == dimension(u2))

@inline firstequal(arg1::AbstractUnitLike) = arg1

@inline function firstequal(arg1::AbstractUnitLike, arg2::AbstractUnitLike)
    (new1, new2) = promote(arg1, arg2)
    (new1 == new2) || throw(DimensionError((arg1, arg2)))
    return new1 
end

@inline function firstequal(arg1::AbstractUnitLike, arg2::AbstractUnitLike, args::AbstractUnitLike...) 
    newargs = promote(arg1, arg2, args...)
    return allequal(newargs) ? first(newargs) : throw(DimensionError(args))
end


#=============================================================================================
 Mathematical operations on quantities
=============================================================================================#
"""
    apply2quantities(f, args::UnionQuantity...)

Converts all arguments to base units, and applies `f` to values and dimensions
returns the narrowest possible quanity (usually RealQuantity). Useful for defining 
new functions for quantities 
    
Thus, in order to support `f` for quantities, simply define 

f(dims::AbstractDimensions...)
f(args::UnionQuantity...) = apply2quantities(f, args...)
"""
@inline function apply2quantities(f, args::UnionQuantity...)
    baseargs = map(ubase, args)
    basevals = map(ustrip, baseargs)
    basedims = map(unit, baseargs)
    return quantity(f(basevals...), f(basedims...))
end


function Base.:(==)(q1::UnionQuantity, q2::UnionQuantity)
    qb1 = ubase(q1)
    qb2 = ubase(q2)
    return (ustrip(qb1) == ustrip(qb2)) && (unit(qb1) == unit(qb2))
end

function Base.:(≈)(q1::UnionQuantity, q2::UnionQuantity)
    qb1 = ubase(q1)
    qb2 = ubase(q2)
    return (ustrip(qb1) ≈ ustrip(qb2)) && (unit(qb1) == unit(qb2))
end

Base.:+(q::UnionQuantity, n::Number) = dimensionless(q) + n 
Base.:+(n::Number, q::UnionQuantity) = dimensionless(q) + n
Base.:+(q1::UnionQuantity, q2::UnionQuantity) = apply2quantities(+, ubase(q1), ubase(q2))
Base.:+(q1::UnionQuantity, qN::UnionQuantity...) = apply2quantities(+, q1, qN...)

Base.:-(q::UnionQuantity, n::Number) = dimensionless(q) - n 
Base.:-(n::Number, q::UnionQuantity) = n - dimensionless(q)
Base.:-(q1::UnionQuantity, q2::UnionQuantity) = apply2quantities(-, ubase(q1), ubase(q2))
Base.:-(q1::UnionQuantity) = apply2quantities(-, ubase(q1))

Base.:*(q0::UnionQuantity, n::Number) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))
Base.:*(n::Number, q0::UnionQuantity) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))
Base.:*(q1::UnionQuantity, q2::UnionQuantity) = apply2quantities(*, q1, q2)
Base.:*(q1::UnionQuantity, qN::UnionQuantity...) = apply2quantities(*, q1, qN...)

Base.:/(q0::UnionQuantity, n::Number) = (q = ubase(q0); quantity(ustrip(q)/n, unit(q)))
Base.:/(n::Number, q0::UnionQuantity) = (q = ubase(q0); quantity(n/ustrip(q), inv(unit(q))))
Base.:/(q1::UnionQuantity, q2::UnionQuantity) = apply2quantities(/, q1, q2)
Base.:inv(q::UnionQuantity) = apply2quantities(inv, q)

for NT in (Number, Complex, Real, Rational, Integer) #Address potential ambiguities
    @eval Base.:^(q::UnionQuantity, p::$NT) = apply2quantities(Base.Fix2(^, p), q)
end

@inline Base.literal_pow(::typeof(^), q::UnionQuantity, ::Val{p}) where {p} = apply2quantities(x->Base.literal_pow(^, x, Val(dimensionless(p))), q)

Base.sqrt(q::UnionQuantity) = apply2quantities(sqrt, q)
Base.cbrt(q::UnionQuantity) = apply2quantities(cbrt, q)


#Functions that return the same unit
for f in (
        :float, :abs, :real, :imag, :conj, :adjoint,
        :transpose, :significand, :zero, :one, :typemax
    )
    @eval Base.$f(u::AbstractDimensions) = u
    @eval Base.$f(q::UnionQuantity) = apply2quantities($f, q)
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
    @eval Base.$f(q::UnionQuantity) = apply2quantities($f, q)
end

#Single-argument functions that only operate on values
for f in (:isfinite, :isinf, :isnan, :isreal, :isempty)
    @eval Base.$f(q::UnionQuantity) = $f(ustrip(q))
end

#Single-argument functions with dimensionless output that don't work with non-scalar units
for f in (:iszero, :signbit, :angle)
    @eval Base.$f(q::UnionQuantity) = (assert_scalar(unit(q)); $f(ustrip_base(q)))
end

#=
#Preliminary test code
K  = Dimensions(temperature=1)
°C = AffineUnits(scale=1, offset=273.15, dims=K, symbol=:°C)
°F = AffineUnits(scale=5/9, offset=(273.15-32*5/9), dims=K, symbol=:°F)


uconvert(°F, Quantity(-40, °C))
uconvert(K, Quantity(0, °F))
uconvert(°C, Quantity(-0, °F))


convert(Quantity{Float64, ScalarUnits}, RealQuantity(-0, °F))
convert(NumberQuantity{Float64, Dimensions}, Quantity(-0, °F))
convert(Quantity{Float64, AffineUnits}, Quantity(-0, °F))
convert(RealQuantity{Float64, AffineUnits}, RealQuantity(-0, °F))
convert(RealQuantity{Float64, AffineUnits}, RealQuantity(-0.0, °F))
=#




#Test code
#=
velocity = Dimensions(length=1, time=-1)
velocity*velocity
velocity^2
velocity^(-0.5)
=#