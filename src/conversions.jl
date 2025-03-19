#============================================================================================
Unit Transforms are the backbone of the conversion process
    uconvert(u1::AbstractUnit, u2::AbstractUnit) produces a kind of unit transform 
    unit transforms can be applied directly to quantities to do the conversion 
Currently, only conversions between affine units are supported, 
    but other conversions (such as log units) are conceivably possible but they would
    require a separate unit registry (which is why I made separate registries easier)
============================================================================================#
abstract type AbstractUnitTransform end

"""
    AffineTransform

A type representing an affine transfomration that can be
used to convert values from one unit to another. This is the
output type of `uconvert(u::AbstractUnitLike, u0::AbstractUnitLike)`

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

Base.broadcastable(trans::AffineTransform) = Ref(trans)
transform(x, trans::AffineTransform) = muladd(x, trans.scale, trans.offset)
transform(x::AbstractArray, trans::AffineTransform) = transform.(x, trans)

"""
    |>(u1::AbstractUnitLike, u2::Union{AbstractUnitLike, UnionQuantity})

Using `q |> qout` is an alias for `uconvert(u, q)`.
"""
Base.:(|>)(u0::AbstractUnitLike, u::AbstractUnitLike) = uconvert(u, u0)
Base.:(|>)(q::UnionQuantity, u::AbstractUnitLike,) = uconvert(u, q)

"""
    uconvert(utarget::AbstractUnitLike, ucurrent::AbstractUnitLike)

Produces an AffineTransform that can convert a quantity with `current`
units to a quantity with `target` units
"""
function uconvert(utarget::AbstractUnitLike, ucurrent::AbstractUnitLike)
    dimension(utarget) == dimension(ucurrent) || throw(DimensionError(first(args), Base.tail(args)))
    return AffineTransform(
        scale  = uscale(ucurrent)/uscale(utarget),
        offset = (uoffset(ucurrent) - uoffset(utarget))/uscale(utarget)
    )
end

"""
    uconvert(u::AbstractUnitLike, q::UnionQuantity)

Converts quantity `q` to the equivalent quantity having units `u`
"""
function uconvert(u::AbstractUnitLike, q::UnionQuantity)
    newval = transform(ustrip(q), uconvert(u, unit(q)))
    return quantity(newval, u)
end

"""
    ubase(q::UnionQuantity)

Converts quantity `q` to its raw dimensional equivalent (such as SI units)
"""
ubase(q::UnionQuantity) = uconvert(dimension(q), q)
ubase(q::UnionQuantity{<:Any,<:AbstractDimensions}) = q
ubase(q::Number) = NoDims()

"""
    ustrip(u::AbstractUnitLike, q::UnionQuantity)

Converts quantity `q` to units `q`` equivalent and removes units
"""
ustrip(u::AbstractUnitLike, q::UnionQuantity) = transform(ustrip(q), uconvert(u, unit(q)))

"""
    basestrip(q::UnionQuantity)

Converts quantity `q` to its raw dimensional equivalent and removes units
"""
ustrip_base(q::UnionQuantity) = ustrip(dimension(q), q)

"""
    ustrip_dimensionless(q::UnionQuantity)

Converts quantity `q` to its raw dimensional equivalent, asserts 
a dimensionless result, and then removes units
"""
ustrip_dimensionless(q::UnionQuantity) = ustrip(assert_dimensionless(ubase(q)))

"""
    asunit(q::UnionQuantity{<:Number})

Convert quantity `q` into a unit of the same magnitude
"""
asunit(q::UnionQuantity{<:Number, <:Dimensions})  = ScalarUnits(scale=ustrip(q), dims=dimension(q))
asunit(q::UnionQuantity{<:Number, <:ScalarUnits}) = ( u = unit(q); ScalarUnits(scale=ustrip(q)*uscale(u), dims=dimension(q)) )
asunit(q::UnionQuantity{<:Number, <:AffineUnits}) = ( u = unit(q); AffineUnits(scale=ustrip(q)*uscale(u), offset=uoffset(u), dims=dimension(q)) )
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


#================================ Conversion between quantity types =======================================#
function Base.convert(::Type{Quantity{T,U}}, q::UnionQuantity) where {T, U<:AbstractUnitLike}
    u = closest_unit(U, unit(q))
    v = transform(ustrip(q), uconvert(u, unit(q)))
    return Quantity{T,U}(v,u)
end

function Base.convert(::Type{Q}, q::UnionQuantity) where {T, U<:AbstractUnitLike, Q<:UnionQuantity{T,U}}
    return constructorof(Q){T,U}(convert(Quantity{T,U}, q))
end

Base.convert(::Type{Quantity{T,U}}, q::Quantity{T,U}) where {T, U<:AbstractUnitLike} = q # Remove potential ambiguities

# NoDims handling ==============================================================================================
Base.convert(::Type{D}, u::NoDims) where {T, D<:AbstractDimensions{T}} = D{T}()
Base.promote_rule(::Type{D}, ::Type{<:NoDims}) where D<:AbstractDimensions = D
Base.promote_rule(::Type{<:NoDims}, ::Type{D}) where D<:AbstractDimensions = D

# Converting between unit types ==============================================================
Base.convert(::Type{U}, u0::AbstractUnitLike) where U<:AffineUnits = U(scale=uscale(u0), offset=uoffset(u0), dims=dimension(u0), symbol=usymbol(u0))
Base.convert(::Type{U}, u0::AbstractUnitLike) where U<:ScalarUnits = (assert_scalar(u0); U(scale=uscale(u0), dims=dimension(u0), symbol=usymbol(u0)))
Base.convert(::Type{U}, u0::AbstractUnitLike) where U<:AbstractDimensions = (assert_dimension(u0); dims=dimension(u0))

# Converting between units and quantities ====================================================
Base.convert(::Type{U}, q::UnionQuantity) where U<:AbstractUnits = convert(U, asunit(q))
function Base.convert(::Type{Q}, u::AbstractUnitLike) where {Q<:UnionQuantity{<:Any, <:AbstractDimensions}}
    assert_scalar(u)
    return Q(uscale(u), dimension(u))
end

# Switching dimension types on affine units (for registries) ======================================================
function Base.convert(::Type{U}, u0::AbstractDimensions{D0}) where {D,D0,U<:AffineUnits{D}}
    return constructorof(U)(dims=convert(D, dimension(u0)))
end

closest_unit(::Type{U}, u::AbstractUnitLike) where U<:AbstractDimensions  = constructorof(U)(dimension(u))
closest_unit(::Type{U}, u::AbstractUnitLike) where U<:AbstractScalarUnits = constructorof(U)(scale=uscale(u), dims=dimension(u))
closest_unit(::Type{U}, u::AbstractUnitLike) where U<:AbstractAffineUnits = constructorof(U)(scale=uscale(u), offset=uoffset(u), dims=dimension(u))


# Promotion rules ======================================================
function Base.promote_rule(::D{P1}, ::D{P2}) where {D<:AbstractDimension, P1, P2}
    return D{promote_type(P1,P2)}
end

#=
function Base.promote_rule(::D1, ::ScalarUnits{D2}) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return ScalarUnits{promote_type(D1, D2)}
end

function Base.promote_rule(::ScalarUnits{D1}, ::D2) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return ScalarUnits{promote_type(D1, D2)}
end

function Base.promote_rule(::ScalarUnits{D1}, ::ScalarUnits{D2}) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return ScalarUnits{promote_type(D1, D2)}
end

function Base.promote_rule(::D1, ::AffineUnits{D2}) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return AffineUnits{promote_type(D1, D2)}
end

function Base.promote_rule(::AffineUnits{D1}, ::D2) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return AffineUnits{promote_type(D1, D2)}
end

function Base.promote_rule(::AffineUnits{D1}, ::ScalarUnits{D2}) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return AffineUnits{promote_type(D1, D2)}
end

function Base.promote_rule(::ScalarUnits{D1}, ::AffineUnits{D2}) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return AffineUnits{promote_type(D1, D2)}
end

function Base.promote_rule(::AffineUnits{D1}, ::AffineUnits{D2}) where {D1<:AbstractDimensions, D2<:AbstractDimensions}
    return AffineUnits{promote_type(D1, D2)}
end

# Generate promotion rules for affine dimensions
let priority_nums = (Dimensions=0, ScalarUnits=1, AffineUnits=2)
    max_priority(D1::Type, D2::Type) = (priority_nums[Symbol(D1)] >= priority_nums[Symbol(D2)]) ? D1 : D2

    for D1 in (:Dimensions, :ScalarUnits, :AffineUnits), D2 in (:AffineDimensions, :Dimensions, :SymbolicDimensions)
        # Skip if both are not affine dimensions
        (D1 != :AffineDimensions && D2 != :AffineDimensions) && continue

        # Determine the output type
        OUT_D = (D1 == :AffineDimensions == D2) ? :AffineDimensions : :Dimensions

        @eval function Base.promote_rule(::Type{$D1{R1}}, ::Type{$D2{R2}}) where {R1,R2}
            return $OUT_D{promote_type(R1,R2)}
        end
    end
end
=#


#=
Preliminary test code
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