#============================================================================================
Unit Converters are the backbone of the conversion process
    uconvert(u1::AbstractUnit, u2::AbstractUnit) produces a kind of unit transformation forumula
    unit transformation formulas can be applied directly to quantities to do the conversion 
Currently, only conversions between affine units are supported, 
    but other conversions (such as log units) are conceivably possible but they would
    require a separate unit registry (which is why I made separate registries easier)
============================================================================================#
"""
An abstract object representing a unit conversion formula. 
Any object that subtypes this is made callable.

# Callable form 
utrans = uconvert(u"°C", u"°F")
utrans(0.0)
31.999999999999986

# Shorthand callable form (syntactic sugar)
(u"°C" |> u"°F")(0.0)
31.999999999999986
"""
abstract type AbstractUnitConverter end

Base.broadcastable(utrans::AbstractUnitConverter) = Ref(utrans)
(utrans::AbstractUnitConverter)(x) = uconvert(utrans, x)
uconvert(utrans::AbstractUnitConverter, x) = ArgumentError("Conversion formulas not yet implemented for $(utrans)")


"""
    uconvert(u::AbstractUnitLike, q::AbstractQuantity)

Converts quantity `q` to the equivalent quantity having units `u`
"""
function uconvert(u::AbstractUnitLike, q::AbstractQuantity)
    utrans = uconvert(u, unit(q))
    return Quantity(utrans(ustrip(q)), u)
end

"""
    dconvert(u::AbstractUnitLike, q::AbstractQuantity)

Converts quantity `q` to the equivalent dimensional quantity having the same dimensions as `u`
"""
dconvert(u::AbstractUnitLike, q::AbstractQuantity) = uconvert(dimension(u), q)


"""
    ubase(q::AbstractQuantity)

Converts quantity `q` to its raw dimensional equivalent (such as SI units)
"""
ubase(q::AbstractQuantity{<:Any,<:AbstractDimLike}) = q


"""
    |>(u1::AbstractUnitLike, u2::Union{AbstractUnitLike, AbstractQuantity})

Using `q |> qout` is an alias for `uconvert(u, q)`.
"""
Base.:(|>)(u0::AbstractUnitLike, u::AbstractUnitLike) = uconvert(u, u0)
Base.:(|>)(q::AbstractQuantity, u::AbstractUnitLike,) = uconvert(u, q)

"""
    ustrip(u::AbstractUnitLike, q::AbstractQuantity)

Converts quantity `q` to units `q`` equivalent and removes units
"""
ustrip(u::AbstractUnitLike, q::AbstractQuantity) = uconvert(u, unit(q))(ustrip(q))
ustrip(u::AbstractArray{<:AbstractUnitLike}, q)  = ustrip.(u, q)

"""
    dstrip(q::AbstractQuantity)

Converts quantity `q` to its raw dimensional equivalent and removes units
"""
dstrip(q::AbstractQuantity) = ustrip(ubase(q))
ustrip_base(q::AbstractQuantity) = dstrip(q)


"""
    ustrip_dimensionless(q::AbstractQuantity)

Converts quantity `q` to its raw dimensional equivalent, asserts 
a dimensionless result, and then removes units
"""
ustrip_dimensionless(q::AbstractQuantity) = ustrip(assert_dimensionless(ubase(q)))


#============================================================================================
Affine transformations
============================================================================================#
"""
    AffineConverter

A type representing an affine transfomration that can be
used to convert values from one unit to another. This is the
output type of `uconvert(u::AbstractUnitLike, u0::AbstractUnitLike)`.
This object is callable.

# Fields
- scale :: Float64
- offset :: Float64

# Constructors
- AffineConverter(scale::Real, offset::Real)
- AffineConverter(; scale, offset)
"""
struct AffineConverter <: AbstractUnitConverter
    scale :: Float64
    offset :: Float64
end
AffineConverter(;scale=1, offset=0) = AffineConverter(scale, offset)

"""
    uconvert(utarget::AbstractAffineLike, ucurrent::AbstractAffineLike)

Produces an AffineConverter that can convert a quantity with `current`
units to a quantity with `target` units
"""
function uconvert(utarget::AbstractAffineLike, ucurrent::AbstractAffineLike)
    dimension(utarget) == dimension(ucurrent) || throw(ConversionError(utarget, ucurrent))

    return AffineConverter(
        scale  = uscale(ucurrent)/uscale(utarget),
        offset = (uoffset(ucurrent) - uoffset(utarget))/uscale(utarget)
    )
end

"""
    uconvert(utrans::AffineConverter, x)

Apply the unit conversion formula "utrans" to value "x", can be used to apply unit conversions on unitless values

julia> uconvert(u"m/s"|>u"km/hr", 1.0)
3.5999999999999996
"""
uconvert(utrans::AffineConverter, x) = muladd(x, utrans.scale, utrans.offset)
uconvert(utrans::AffineConverter, x::AbstractArray) = utrans.(x)
uconvert(utrans::AffineConverter, x::Tuple) = map(utrans, x)

function ubase(q::AbstractQuantity{<:Any,<:AbstractAffineUnits})
    u = unit(q)
    utrans = AffineConverter(scale=uscale(u), offset=uoffset(u))
    return Quantity(utrans(ustrip(q)), dimension(u))
end 



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

# Converting between units and quantities ====================================================
Base.convert(::Type{U}, q::AbstractQuantity) where U<:AbstractUnits = convert(U, asunit(q))
function Base.convert(::Type{Q}, u::AbstractUnitLike) where {Q<:AbstractQuantity{<:Any, <:AbstractDimensions}}
    assert_scalar(u)
    return Q(uscale(u), dimension(u))
end

# Converting between generic numbers and quantity types 
Base.convert(::Type{T}, q::AbstractQuantity) where {T<:Union{Number,AbstractArray}} = convert(T, dimensionless(q))
Base.convert(::Type{Q}, x::Union{Number,AbstractArray}) where {Q<:AbstractQuantity} = Q(x, dimtype(Q)(0))

# Promotion rules ======================================================
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
function Base.promote_rule(::Type{Quantity{T1,U1}}, ::Type{Quantity{T2,U2}}) where {T1, T2, U1<:AbstractUnitLike, U2<:AbstractUnitLike}
    D = promote_type(dimtype(U1), dimtype(U2))
    T = promote_type(T1, T2)
    return Quantity{T, D}
end

#Conversion to dimensions doesn't happen if quantities only differ by value type
function Base.promote_rule(::Type{Quantity{T1,U}}, ::Type{Quantity{T2,U}}) where {T1, T2, U<:AbstractUnitLike}
    return Quantity{promote_type(T1, T2), U}
end

#Cases where values are updated to quantities
Base.promote_rule(::Type{Quantity{T,U}}, ::Type{T}) where {T,U<:AbstractUnitLike} = Quantity{T}
Base.promote_rule(::Type{Quantity{T1,U}}, ::Type{T2}) where {T1<:Number, T2<:Number, U<:AbstractUnitLike} = Quantity{promote_type(T1,T2), U}
Base.promote_rule(::Type{Quantity{T1,U}}, ::Type{T2}) where {T1<:AbstractArray, T2<:AbstractArray, U<:AbstractUnitLike} = Quantity{promote_type(T1,T2), U}
