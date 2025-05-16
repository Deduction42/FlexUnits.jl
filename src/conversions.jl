#============================================================================================
Unit Transforms are the backbone of the conversion process
    uconvert(u1::AbstractUnit, u2::AbstractUnit) produces a kind of unit transformation
    unit transformations can be applied directly to quantities to do the conversion 
Currently, only conversions between affine units are supported, 
    but other conversions (such as log units) are conceivably possible but they would
    require a separate unit registry (which is why I made separate registries easier)
============================================================================================#
"""
An abstract object representing a unit transformation formula. 
Any object that subtypes this is made callable.

# Callable form 
utrans = uconvert(u"°C", u"°F")
utrans(0.0)
31.999999999999986

# Shorthand callable form (syntactic sugar)
(u"°C" |> u"°F")(0.0)
31.999999999999986
"""
abstract type AbstractUnitTransform end

Base.broadcastable(utrans::AbstractUnitTransform) = Ref(utrans)
(utrans::AbstractUnitTransform)(x) = uconvert(utrans, x)
uconvert(utrans::AbstractUnitTransform, x) = ArgumentError("Conversion formulas not yet implemented for $(utrans)")


"""
    |>(u1::AbstractUnitLike, u2::Union{AbstractUnitLike, AbstractQuantity})

Using `q |> qout` is an alias for `uconvert(u, q)`.
"""
Base.:(|>)(u0::AbstractUnitLike, u::AbstractUnitLike) = uconvert(u, u0)
Base.:(|>)(q::AbstractQuantity, u::AbstractUnitLike,) = uconvert(u, q)


#============================================================================================
Affine transformations
============================================================================================#
"""
    AffineTransform

A type representing an affine transfomration that can be
used to convert values from one unit to another. This is the
output type of `uconvert(u::AbstractUnitLike, u0::AbstractUnitLike)`.
This object is callable.

# Fields
- scale :: Float64
- offset :: Float64

# Constructors
- AffineTransform(scale::Real, offset::Real)
- AffineTransform(; scale, offset)
"""
struct AffineTransform <: AbstractUnitTransform
    scale :: Float64
    offset :: Float64
end
AffineTransform(;scale, offset) = AffineTransform(scale, offset)

"""
    uconvert(utrans::AffineTransform, x)

Apply the unit conversion formula "utrans" to value "x", can be used to apply unit conversions on unitless values

julia> uconvert(u"m/s"|>u"km/hr", 1.0)
3.5999999999999996
"""
uconvert(utrans::AffineTransform, x) = muladd(x, utrans.scale, utrans.offset)
uconvert(utrans::AffineTransform, x::AbstractArray) = utrans.(x)
uconvert(utrans::AffineTransform, x::Tuple) = map(utrans, x)

"""
    uconvert(utarget::AbstractAffineLike, ucurrent::AbstractAffineLike)

Produces an AffineTransform that can convert a quantity with `current`
units to a quantity with `target` units
"""
function uconvert(utarget::AbstractAffineLike, ucurrent::AbstractAffineLike)
    dimension(utarget) == dimension(ucurrent) || throw(ConversionError(utarget, ucurrent))

    return AffineTransform(
        scale  = uscale(ucurrent)/uscale(utarget),
        offset = (uoffset(ucurrent) - uoffset(utarget))/uscale(utarget)
    )
end

"""
    uconvert(u::AbstractUnitLike, q::AbstractQuantity)

Converts quantity `q` to the equivalent quantity having units `u`
"""
function uconvert(u::AbstractUnitLike, q::AbstractQuantity)
    utrans = uconvert(u, unit(q))
    return quantity(utrans(ustrip(q)), u)
end

"""
    ubase(q::AbstractQuantity)

Converts quantity `q` to its raw dimensional equivalent (such as SI units)
"""
ubase(q::AbstractQuantity{<:Any,<:AbstractDimensions}) = q
function ubase(q::AbstractQuantity{<:Any,<:AbstractAffineUnits})
    u = unit(q)
    utrans = AffineTransform(scale=uscale(u), offset=uoffset(u))
    return quantity(utrans(ustrip(q)), dimension(u))
end 

"""
    ustrip(u::AbstractUnitLike, q::AbstractQuantity)

Converts quantity `q` to units `q`` equivalent and removes units
"""
ustrip(u::AbstractUnitLike, q::AbstractQuantity) = uconvert(u, unit(q))(ustrip(q))
ustrip(u::AbstractArray{<:AbstractUnitLike}, q)  = ustrip.(u, q)

"""
    ustrip_base(q::AbstractQuantity)

Converts quantity `q` to its raw dimensional equivalent and removes units
"""
ustrip_base(q::AbstractQuantity) = ustrip(ubase(q))

"""
    ustrip_dimensionless(q::AbstractQuantity)

Converts quantity `q` to its raw dimensional equivalent, asserts 
a dimensionless result, and then removes units
"""
ustrip_dimensionless(q::AbstractQuantity) = ustrip(assert_dimensionless(ubase(q)))

"""
    asunit(q::AbstractQuantity{<:Number})

Convert quantity `q` into a unit of the same magnitude
"""
asunit(q::AbstractQuantity{<:Number, <:Dimensions})  = AffineUnits(scale=ustrip(q), dims=dimension(q))
asunit(q::AbstractQuantity{<:Number, <:AffineUnits}) = ( u = unit(q); AffineUnits(scale=ustrip(q)*uscale(u), offset=uoffset(u), dims=dimension(q)) )
asunit(u::AbstractUnitLike) = u

#=================================================================================================
Conversion and Promotion

The conversion/promotion strategy works according to the engineering workflow:
convert units to SI => do calculations => convert to desired units

The main reason behind this is that Dimensions are the fastest to compute around and easiest to work with.
Moreover, only defining calculations for dimenional units also greatly simplifies the code, thus all operations 
will promote to dimensional units (such as SI). WARNING:: This also means that Affine Units will auto-convert
this can yield potentially unintuitive results like 2°C/1°C = 1.0036476381542951
=================================================================================================#


# Converting quantity types ====================================================
Base.convert(::Type{Q}, q::AbstractQuantity) where Q<:AbstractQuantity = (q isa Q) ? q : Q(ustrip(q), unit(q))
function Base.convert(::Type{Q}, q::AbstractQuantity) where Q<:AbstractQuantity{<:Any, <:AbstractDimensions}
    if (q isa Q) 
        return q 
    else 
        qb = ubase(q)
        return Q(ustrip(qb), unit(qb))
    end
end

# Converting unit types ====================================================
Base.convert(::Type{U}, u::AbstractUnitLike) where U<:AffineUnits = (u isa U) ? u : U(scale=uscale(u), offset=uoffset(u), dims=dimension(u), symbol=usymbol(u))
Base.convert(::Type{D}, u::AbstractUnitLike) where D<:AbstractDimensions = (u isa D) ? u : D(dimension(assert_dimension(u)))
Base.convert(::Type{D}, u::NoDims) where D<:AbstractDimensions = D()

# Converting between units and quantities ====================================================
Base.convert(::Type{U}, q::AbstractQuantity) where U<:AbstractUnits = convert(U, asunit(q))
function Base.convert(::Type{Q}, u::AbstractUnitLike) where {Q<:AbstractQuantity{<:Any, <:AbstractDimensions}}
    assert_scalar(u)
    return Q(uscale(u), dimension(u))
end

# Converting between generic numbers and quantity types 
Base.convert(::Type{T}, q::AbstractQuantity) where {T<:Number} = convert(T, dimensionless(q))

# Promotion rules ======================================================

#We assume dimension types match except for NoDims
function Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:AbstractDimensions, D2<:NoDims}
    return constructorof(D1){promote_type(eltype(D1), eltype(D2))}
end

function Base.promote_rule(::Type{<:Dimensions{P1}}, ::Type{<:Dimensions{P2}}) where {P1, P2}
    return Dimensions{promote_type(P1,P2)}
end

#Unit promotion (favors AffineDimensiosn due to information loss)
function Base.promote_rule(::Type{D1}, ::Type{AffineUnits{D2}}) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return AffineUnits{promote_type(D1, D2)}
end
function Base.promote_rule(::Type{AffineUnits{D1}}, ::Type{AffineUnits{D2}}) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return AffineUnits{promote_type(D1, D2)}
end

#Quantity promotion (favors Dimensions, as converting quantities to SI does not result in information loss)
function Base.promote_rule(::Type{Quantity{T1,U1}}, ::Type{Quantity{T2,U2}}) where {T1, T2, U1, U2}
    D = promote_type(dimtype(U1), dimtype(U2))
    T = promote_type(T1, T2)
    return Quantity{T, D}
end


