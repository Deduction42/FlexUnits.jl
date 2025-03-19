#=============================================================================================
 Mathematical operations on dimensions
=============================================================================================#
"""
    map_dimensions(f::F, args::AbstractDimensions...)

Similar to the `map` function, but specifically iterates over dimensions.
Useful for defining mathematical operations for dimensions
"""
function map_dimensions(f::F, args::AbstractDimensions...) where {F<:Function}
    D = promote_type(typeof(args).parameters...)
    return  D(
        ( f((getproperty(arg, name) for arg in args)...) for name in dimension_names(D) )...
    )
end

Base.:+(args::AbstractDimensions...) = firstequal(args...)
Base.:-(args::AbstractDimensions...) = firstequal(args...)
Base.:*(args::AbstractDimensions...) = map_dimensions(+, args...)
Base.:/(args::AbstractDimensions...) = map_dimensions(-, args...)
Base.inv(arg::AbstractDimensions) = map_dimensions(-, arg)
Base.:^(d::AbstractDimensions, p::Integer) = map_dimensions(Base.Fix1(*, p), d)
Base.:^(d::AbstractDimensions{R}, p::Number) where {R} = map_dimensions(Base.Fix1(*, tryrationalize(R, p)), d)
@inline Base.literal_pow(::typeof(^), d::AbstractDimensions, ::Val{p}) where {p} = map_dimensions(Base.Fix1(*, p), d)
Base.sqrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 2))
Base.cbrt(d::AbstractDimensions{R}) where R = d^inv(convert(R, 3))
Base.adjoint(d::AbstractDimensions) = d
Base.abs2(d::AbstractDimensions) = d^2


#=============================================================================================
 Mathematical operations on abstract units (mostly for parsing)
=============================================================================================#
Base.:*(n::Number, u::AbstractUnitLike) = quantity(n, u)
Base.:/(n::Number, u::AbstractUnitLike) = quantity(n, inv(u))

Base.:*(q::UnionQuantity, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)*u)
Base.:/(q::UnionQuantity, u::AbstractUnitLike) = quantity(ustrip(q), unit(q)/u)

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
    return constructorof(U)(scale=uscale(u)^p, dims=scalar_dimension(u)^p)
end

Base.sqrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 2))
Base.cbrt(u::AbstractUnits{D}) where {R, D<:AbstractDimensions{R}} = u^inv(convert(R, 3))

#Equality does not compare symbols
Base.:(==)(u1::AbstractUnits, u2::AbstractUnits) = (uscale(u1) == uscale(u2)) & (uoffset(u1) == uoffset(u2)) & (dimension(u1) == dimension(u2))
Base.:(==)(u1::AbstractScalarUnits, u2::AbstractScalarUnits) = (uscale(u1) == uscale(u2)) & (dimension(u1) == dimension(u2))

function firstequal(args::AbstractUnitLike...) 
    newargs = promote(args...)
    return allequal(newargs) ? first(newargs) : throw(DimensionError(first(args), Base.tail(args)))
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
function apply2quantities(f, args::UnionQuantity...)
    baseargs = map(ubase, args)
    basevals = map(ustrip, baseargs)
    basedims = map(unit, baseargs)
    return quantity(f(basevals...), f(basedims...))
end

Base.:+(q::UnionQuantity, n::Number) = ustrip(assert_dimensionless(ubase(q))) + n 
Base.:+(n::Number, q::UnionQuantity) = ustrip(assert_dimensionless(ubase(q))) + n
Base.:+(q1::UnionQuantity, q2::UnionQuantity) = apply2quantities(+, ubase(q1), ubase(q2))
Base.:+(q::UnionQuantity...) = apply2quantities(+, q...)

Base.:-(q::UnionQuantity, n::Number) = ustrip(assert_dimensionless(ubase(q))) - n 
Base.:-(n::Number, q::UnionQuantity) = n - ustrip(assert_dimensionless(ubase(q)))
Base.:-(q1::UnionQuantity, q2::UnionQuantity) = apply2quantities(-, ubase(q1), ubase(q2))
Base.:-(q1::UnionQuantity) = apply2quantities(-, ubase(q1))

Base.:*(q0::UnionQuantity, n::Number) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))
Base.:*(n::Number, q0::UnionQuantity) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))
Base.:*(q1::UnionQuantity, q2::UnionQuantity) = apply2quantities(*, q1, q2)
Base.:*(q::UnionQuantity...) = apply2quantities(*, q...)

Base.:/(q0::UnionQuantity, n::Number) = (q = ubase(q0); quantity(ustrip(q)/n, unit(q)))
Base.:/(n::Number, q0::UnionQuantity) = (q = ubase(q0); quantity(n/ustrip(q), inv(unit(q))))
Base.:/(q1::UnionQuantity, q2::UnionQuantity) = apply2quantities(/, q1, q2)
Base.:inv(q::UnionQuantity) = apply2quantities(inv, q)

for NT in (Number, Integer, Rational) #Address potential ambiguities
    @eval Base.:^(q::UnionQuantity, p::$NT) = apply2quantities(Base.Fix2(^, p), q)
end
@inline Base.literal_pow(::typeof(^), q::UnionQuantity, ::Val{p}) where {p} = apply2quantities(x->Base.literal_pow(^, x, Val(p)), q)

Base.sqrt(q::UnionQuantity) = apply2quantities(sqrt, q)
Base.cbrt(q::UnionQuantity) = apply2quantities(cbrt, q)


#Functions that return the same unit
for f in (
        :float, :abs, :real, :imag, :conj, :adjoint, :unsigned,
        :nextfloat, :prevfloat, :transpose, :significand, :zero, :one
    )
    @eval Base.$f(q::UnionQuantity) = quantity($f(ustrip(q)), unit(q))
end

#Functions that strip dimensional units
for f in (:angle,)
    @eval Base.$f(q::UnionQuantity) = f(ustrip_base(q))
end

#Single-argument funtions that require dimensionless input
for f in (
    :sin, :cos, :tan, :sinh, :cosh, :tanh, :asin, :acos,
    :asinh, :acosh, :atanh, :sec, :csc, :cot, :asec, :acsc, :acot, :sech, :csch,
    :coth, :asech, :acsch, :acoth, :sinc, :cosc, :cosd, :cotd, :cscd, :secd,
    :sinpi, :cospi, :sind, :tand, :acosd, :acotd, :acscd, :asecd, :asind,
    :log, :log2, :log10, :log1p, :exp, :exp2, :exp10, :expm1, :frexp, :exponent,
)
    @eval Base.$f(q::UnionQuantity) = $f(ustrip_dimensionless(q))
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