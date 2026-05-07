#============================================================================================
uconvert with transformation objects
============================================================================================#
#Generic transformation generator
"""
    uconvert(u_target::AbstractUnitLike, u_current::AbstractUnitLike) 

Produces a conversion rule to convert `u_current` to `u_target`. This result is a callable object that
can be applied directly to numeric objects.

```julia
julia> uconvert(u"K", u"°C")
AffineTransform{Float64}(1.0, 273.15)

julia> uconvert(u"K", u"°C")(0)
273.15
```
"""
function uconvert(u_target::AbstractUnitLike, u_current::AbstractUnitLike) 
    assert_convertable(u_target, u_current)
    return inv(tobase(u_target)) ∘ tobase(u_current)
end
function uconvert(u_target::AbstractUnitLike, u_current::AbstractDimLike)
    assert_convertable(u_target, u_current)
    return Base.Fix1(inv, tobase(u_target))
end
function uconvert(u_target::AbstractDimLike, u_current::AbstractUnitLike)
    assert_convertable(u_target, u_current)
    return tobase(u_current)
end
function uconvert(u_target::AbstractDimLike, u_current::AbstractDimLike)
    assert_convertable(u_target, u_current)
    return NoTransform()
end
uconvert(ft::AbstractUnitTransform, x) = ft(x)
uconvert(u_target::StaticDims{D}, u_current::Units{StaticDims{D}}) where D = tobase(u_current)

assert_convertable(u_target::AbstractUnitLike, u_current::AbstractUnitLike) = compatible_dims(u_target, u_current) ? u_target : throw(ConversionError(u_target, u_current))

#============================================================================================
uconvert with quantities
============================================================================================#
"""
    uconvert(u::AbstractUnitLike, q::Union{LogQuant, QuantUnion})

Converts quantity `q` to the equivalent quantity having units `u`. If unit is a logarithmic unit, a LogQuant is returned

```julia
julia> uconvert(u"K", 25u"°C")
298.15 K
```
"""
function uconvert(u::AbstractUnitLike, q::QuantUnion)
    newval = uconvert(u, unit(q))(ustrip(q))
    return Quantity{typeof(newval), typeof(u)}(newval, u)
end

function uconvert(u::AbstractUnitLike, lq::LogQuant)
    assert_convertable(u, unit(lq))
    ft = inv(tobase(u)) ∘ exp(tobase(unit(lq))) #This produces an affine transform of an ExpAffineTransform
    newval = ft(ustrip(lq))
    return Quantity{typeof(newval), typeof(u)}(newval, u)
end

function uconvert(u::Units{D,<:ExpAffTransform}, q::QuantUnion) where D<:AbstractDimLike
    newval = uconvert(u, unit(q))(ustrip(q))
    newunit = Units{D}(dimension(u), log(tobase(u)), usymbol(u))
    return LogQuant{typeof(newval), typeof(newunit)}(newval, newunit)
end

function uconvert(u::Units{D,<:ExpAffTransform}, lq::LogQuant) where D<:AbstractDimLike
    assert_convertable(u, unit(lq))
    ft = inv(log(tobase(u))) ∘ tobase(unit(lq)) #This produces an AffineTransform
    newval = ft(ustrip(lq))
    newunit = Units{D}(dimension(u), log(tobase(u)), usymbol(u))
    return LogQuant{typeof(newval), typeof(newunit)}(newval, newunit)
end

uconvert(::Type{D}, q::QuantUnion) where D <: StaticDims = uconvert(D(), q)
uconvert(::Type{D}, q::LogQuant) where D <: StaticDims = uconvert(D(), q)


"""
    dconvert(u::AbstractUnitLike, q::QuantUnion)

Converts quantity `q` to the equivalent dimensional quantity having the same dimensions as `u`
```julia
julia> dconvert(u"km/hr", 25u"km/hr")
6.944444444444445 m/s
```
"""
dconvert(u::AbstractUnitLike, q::QuantUnion) = uconvert(dimension(u), q)


"""
    |>(u2::Union{AbstractUnitLike, QuantUnion}, u1::AbstractUnitLike)

Using `q |> qout` is an alias for `uconvert(u, q)`.
"""
Base.:(|>)(q::LogQuantUnion, u::AbstractUnitLike) = uconvert(u, q)
Base.:(|>)(q::LogQuantUnion, ::Type{T}) where T <: StaticDims = uconvert(T(), q)
Base.:(|>)(u0::AbstractUnitLike, u::AbstractUnitLike) = uconvert(u, u0)
Base.:(|>)(u0::AbstractUnitLike, ::Type{T}) where T <: StaticDims = uconvert(T(), u0)

"""
    ustrip(u::AbstractUnitLike, q::QuantUnion)

Converts quantity `q` to units `u` equivalent and removes units
"""
ustrip(u::AbstractUnitLike, q::QuantUnion) = uconvert(u, unit(q))(ustrip(q))
ustrip(u::AbstractUnitLike, q::LogQuant) = uconvert(u, dimension(q))(dstrip(q))
ustrip(u::AbstractArray{<:AbstractUnitLike}, q)  = ustrip.(u, q)

"""
    dstrip(u::AbstractUnitLike, q::QuantUnion)

Converts quantity `q` to its raw dimensional equivalent and verifies if its dimensions are consistent with `u`
Leaving out the unit argument `u` skips the dimensioonal verification process
"""
dstrip(u::AbstractUnitLike, q::Union{QuantUnion,LogQuant}) = ustrip(dimension(u), q)

"""
    ustrip_dimensionless(q::QuantUnion)

Converts quantity `q` to its raw dimensional equivalent, asserts 
a dimensionless result, and then removes units
"""
ustrip_dimensionless(q::QuantUnion) = ustrip(assert_dimensionless(ubase(q)))

"""
    Units(q::QuantUnion{<:Number})

Convert quantity `q` into a unit of the same magnitude
"""
Units(q::QuantUnion{<:Number, <:AbstractUnitLike}) = Units(dims=dimension(q), tobase=tobase(unit(q))*ustrip(q), symbol=DEFAULT_USYMBOL)
Units(u::Units) = u

compatible_dims(d_target::AbstractDimLike, d_current::AbstractDimLike) = (d_target == d_current || isunknown(d_current))
compatible_dims(d_target::StaticDims{d1}, d_current::StaticDims{d2}) where {d1, d2} = (d1==d2)
compatible_dims(target, current) = compatible_dims(dimension(target), dimension(current))

#=================================================================================================
Conversion and Promotion

The conversion/promotion strategy works according to the engineering workflow:
convert units to SI => do calculations => convert to desired units

The main reason behind this is that Dimensions are the fastest to compute around and easiest to work with.
Moreover, only defining calculations for dimenional units also greatly simplifies the code, thus all operations 
will promote to dimensional units (such as SI). WARNING:: This also means that Affine Units will auto-convert
this can yield potentially unintuitive results like 2°C/1°C = 1.0036476381542951
=================================================================================================#
Base.convert(::Type{AffineTransform{T}}, t::AffineTransform) where T = AffineTransform{T}(t.scale, t.offset)
Base.convert(::Type{AffineTransform{T}}, t::NoTransform) where T = AffineTransform{T}(1,0)
Base.convert(::Type{ExpAffTransform{T}}, t::ExpAffTransform) where T = ExpAffTransform{T}(t.scale, t.offset)
Base.convert(::Type{NoTransform}, t::AffineTransform) = is_identity(t) ? NoTransform() : throw(ArgumentError("Cannot convert non-identity AffineTransform to a NoTransform"))
Base.convert(::Type{NoTransform}, t::ExpAffTransform) = throw(ArgumentError("Cannot convert an ExpAffTransform to a NoTransform"))

Base.convert(::Type{Quantity{T,U}}, q::LogQuantUnion) where {T,U} = Quantity{T,U}(q)
Base.convert(::Type{FlexQuant{T,U}}, q::LogQuantUnion) where {T,U} = FlexQuant{T,U}(q)
Base.convert(::Type{LogQuant{T,U}}, q::LogQuantUnion) where {T,U} = LogQuant{T,U}(q)

function Base.convert(::Type{Quantity{T,D}}, u::AbstractUnitLike) where {T,D<:AbstractDimensions}
    assert_scalar(u)
    return Quantity(uscale(u), dimension(u))
end

function Base.convert(::Type{FlexQuant{T,D}}, u::AbstractUnitLike) where {T,D<:AbstractDimensions}
    assert_scalar(u)
    return FlexQuant(uscale(u), dimension(u))
end


# Converting unit types ====================================================
Base.convert(::Type{U}, u::AbstractUnitLike) where {T,D,U<:Units{D,T}} = (u isa Units{D,T}) ? u : Units{D,T}(dims=dimension(u), tobase=tobase(u), symbol=usymbol(u))

Base.convert(::Type{D}, u::AbstractUnitLike) where D<:AbstractDimensions = D(dimension(assert_dimension(u)))
Base.convert(::Type{D}, u::AbstractDimLike) where D<:AbstractDimensions = D(u)
Base.convert(::Type{D}, d::StaticDims) where {D<:AbstractDimensions} = convert(D, dimval(d))

Base.convert(::Type{D}, u::AbstractUnitLike) where {D<:StaticDims} = convert(D, dimension(assert_dimension(u)))
Base.convert(::Type{D}, d::AbstractDimensions) where {D<:StaticDims} = assert_convertable(D(), d)
Base.convert(::Type{D}, d::StaticDims) where {D<:StaticDims} = assert_convertable(D(), d)
Base.convert(::Type{D}, d::NoDims) where {D<:StaticDims} = assert_dimensionless(D())

# Converting transform types ===============================================
Base.convert(::Type{T}, t::NoTransform) where T <: Union{<:AffineTransform, <:NoTransform} = T()

# Converting between units and quantities ====================================================
Base.convert(::Type{U}, q::QuantUnion) where U<:Units = Units(q)

# Converting between generic numbers and quantity types 
Base.convert(::Type{T}, q::QuantUnion) where {T<:MathUnion} = convert(T, dimensionless(q))

# Promotion rules ======================================================
Base.promote_rule(::Type{AffineTransform{T1}}, ::Type{AffineTransform{T2}}) where{T1,T2} = AffineTransform{promote_type{T1,T2}}
Base.promote_rule(::Type{AffineTransform{T}}, ::Type{NoTransform}) where{T} = AffineTransform{T}

function Base.promote_rule(::Type{<:Dimensions{P1}}, ::Type{<:Dimensions{P2}}) where {P1, P2}
    return Dimensions{promote_type(P1,P2)}
end
Base.promote_rule(::Type{T1}, ::Type{T2}) where {T1<:NoTransform, T2<:AbstractUnitTransform} = T2

#Conflicting or uncertain static dimensions get promoted to dynamic version
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:AbstractDimensions, D2<:StaticDims} = promote_type(D1, dimvaltype(D2))
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:StaticDims, D2<:StaticDims} = promote_type(dimvaltype(D1), dimvaltype(D2))
Base.promote_rule(::Type{D}, ::Type{NoDims}) where D<:AbstractDimensions = D
Base.promote_rule(::Type{D}, ::Type{NoDims}) where D<:StaticDims = dimvaltype(D)


#Unit promotion
function Base.promote_rule(::Type{D1}, ::Type{Units{D2,T2}}) where {D1<:AbstractDimLike, D2<:AbstractDimLike, T2<:AbstractUnitTransform}
    return Units{promote_type(D1, D2), T2}
end
function Base.promote_rule(::Type{Units{D1,T1}}, ::Type{Units{D2,T2}}) where {D1<:AbstractDimLike, D2<:AbstractDimLike, T1<:AbstractUnitTransform, T2<:AbstractUnitTransform}
    return Units{promote_type(D1, D2), promote_type(T1,T2)}
end

#Quantity promotion (favors Dimensions, as converting quantities to SI does not result in information loss)
function Base.promote_rule(::Type{Q1}, ::Type{Q2}) where {T1, T2, U1<:AbstractUnitLike, U2<:AbstractUnitLike, Q1<:QuantUnion{T1,U1}, Q2<:QuantUnion{T2,U2}}
    D = promote_type(dimtype(U1), dimtype(U2))
    T = promote_type(T1, T2)
    return quant_type(T){T, D}
end

#Conversion to dimensions doesn't happen if quantities only differ by value type
function Base.promote_rule(::Type{Q1}, ::Type{Q2}) where {T1, T2, U<:AbstractUnitLike, Q1<:QuantUnion{T1,U}, Q2<:QuantUnion{T2,U}}
    T = promote_type(T1, T2)
    return quant_type(T){T, U}
end

#Cases where values are updated to quantities
function Base.promote_rule(::Type{Q}, ::Type{T2}) where {T2<:NumUnion, T1<:NumUnion, U<:AbstractUnitLike, Q<:QuantUnion{T1,U}}
    T = promote_type(T1,T2)
    return quant_type(T){T, nodim_promote(U)}
end
function Base.promote_rule(::Type{Q}, ::Type{T2}) where {T2<:AbstractArray, T1<:AbstractArray, U<:AbstractUnitLike, Q<:QuantUnion{T1,U}}
    T = promote_type(T1,T2)
    return quant_type(T){T, nodim_promote(U)}
end

nodim_promote(::Type{U}) where U<:AbstractUnits = nodim_promote(dimtype(U))
nodim_promote(::Type{D}) where D<:StaticDims = isdimensionless(D) ? D : error("Cannot promote NoDims to $(D) because it's not dimensionless")
nodim_promote(::Type{D}) where D<:AbstractDimensions = D

