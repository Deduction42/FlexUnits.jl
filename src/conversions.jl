#============================================================================================
uconvert with transformation objects
============================================================================================#
#Generic transformation generator
function uconvert(u_target::AbstractUnitLike, u_current::AbstractUnitLike) 
    dimension(u_target) == dimension(u_current) || throw(ConversionError(u_target, u_current))
    return inv(todims(u_target)) ∘ todims(u_current)
end
uconvert(ft::AbstractUnitTransform, x) = ft(x)
uconvert(utarget::StaticDims{D}, ucurrent::StaticUnits{D}) where D = ucurrent.todims

#============================================================================================
uconvert with quantities
============================================================================================#
"""
    uconvert(u::AbstractUnitLike, q::AbstractQuantity)

Converts quantity `q` to the equivalent quantity having units `u`
"""
function uconvert(u::AbstractUnitLike, q::AbstractQuantity)
    dimension(u) == dimension(q) || throw(ConversionError(u, unit(q)))
    ft = uconvert(u, unit(q))
    newval = ft(ustrip(q))
    return Quantity{typeof(newval), typeof(u)}(newval, u)
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
function ubase(q::AbstractQuantity{<:Any,<:AbstractUnitLike})
    u  = unit(q)
    ft = todims(u)
    return Quantity(ft(ustrip(q)), dimension(u))
end 
function ubase(q::AbstractQuantity{T, <:StaticUnits{D}}) where {T,D}
    x = unit(q).todims(ustrip(q))
    return Quantity{typeof(x), StaticDims{D}}(x, StaticDims{D}())
end
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

"""
    Units(q::AbstractQuantity{<:Number})

Convert quantity `q` into a unit of the same magnitude
"""
Units(q::AbstractQuantity{<:Number, <:AbstractUnitLike}) = Units(dims=dimension(q), todims=todims(unit(q))*ustrip(q), symbol=DEFAULT_USYMBOL)
Units(u::Units) = u

#=================================================================================================
Conversion and Promotion

The conversion/promotion strategy works according to the engineering workflow:
convert units to SI => do calculations => convert to desired units

The main reason behind this is that Dimensions are the fastest to compute around and easiest to work with.
Moreover, only defining calculations for dimenional units also greatly simplifies the code, thus all operations 
will promote to dimensional units (such as SI). WARNING:: This also means that Affine Units will auto-convert
this can yield potentially unintuitive results like 2°C/1°C = 1.0036476381542951
=================================================================================================#

#Converting dynamic quantity types ====================================================
Base.convert(::Type{Quantity{T,D}}, q::AbstractQuantity) where {T,D<:AbstractDimensions} = Quantity{T,D}(dstrip(q), dimension(q))
Base.convert(::Type{Quantity{T,U}}, q::AbstractQuantity) where {T,U<:Units} = Quantity{T,U}(ustrip(q), unit(q))
Base.convert(::Type{Quantity{T,M}}, q::AbstractQuantity) where {T, D<:AbstractDimensions, M<:MirrorUnion{D}} = Quantity{T,M}(T(dstrip(q)), dimension(q))
Base.convert(::Type{Quantity{T,D}}, x::Union{Number,AbstractArray}) where {T,D<:AbstractDimensions} = Quantity{T,D}(x, D(0))
Base.convert(::Type{Quantity{T,U}}, x::Union{Number,AbstractArray}) where {T,U<:Units} = Quantity{T,U}(x, dimtype(U)(0))

#Converting static quantity types ====================================================
Base.convert(::Type{Quantity{T,D}}, q::AbstractQuantity) where {T, D<:StaticDims} = Quantity{T,D}(ustrip(D(), q), D())
Base.convert(::Type{Quantity{T,U}}, q::AbstractQuantity) where {T, D, C, U<:StaticUnits{D,C}} = Quantity{T,U}(ustrip(D, q), StaticUnits{D,C}(C()))
function Base.convert(::Type{Quantity{T,D}}, u::AbstractUnitLike) where {T,D<:AbstractDimensions}
    assert_scalar(u)
    return Q(uscale(u), dimension(u))
end

# Converting unit types ====================================================
Base.convert(::Type{U}, u::AbstractUnitLike) where {T,D,U<:Units{D,T}} = (u isa Units{D,T}) ? u : Units{D,T}(dims=dimension(u), todims=todims(u), symbol=usymbol(u))
Base.convert(::Type{U}, u::AbstractUnitLike) where {T,D,U<:StaticUnits{D,T}} = (u isa StaticUnits{D,T}) ? u : StaticUnits{D,T}(todims=todims(u), symbol=usymbol(u))

Base.convert(::Type{D}, u::AbstractUnitLike) where D<:AbstractDimensions = D(dimension(assert_dimension(u)))
Base.convert(::Type{D}, u::AbstractDimLike) where D<:AbstractDimensions = D(u)
Base.convert(::Type{D}, d::StaticDims) where {D<:AbstractDimensions} = convert(D, dimval(d))

# Converting transform types ===============================================
Base.convert(::Type{T}, t::NoTransform) where T <: AbstractUnitTransform = T()

# Converting between units and quantities ====================================================
Base.convert(::Type{U}, q::AbstractQuantity) where U<:Units = Units(q)

# Converting between generic numbers and quantity types 
Base.convert(::Type{T}, q::AbstractQuantity) where {T<:Union{Number,AbstractArray}} = convert(T, dimensionless(q))

# Promotion rules ======================================================
function Base.promote_rule(::Type{<:Dimensions{P1}}, ::Type{<:Dimensions{P2}}) where {P1, P2}
    return Dimensions{promote_type(P1,P2)}
end
Base.promote_rule(::Type{T1}, ::Type{T2}) where {T1<:NoTransform, T2<:AbstractUnitTransform} = T2

#Conflicting or uncertain static dimensions get promoted to dynamic version
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:AbstractDimensions, D2<:StaticDims} = promote_type(D1, dimtype(D2))
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:StaticDims, D2<:StaticDims} = promote_type(dimtype(D1), dimtype(D2))
Base.promote_rule(::Type{U1}, ::Type{U2}) where {T1, d1, U1<:StaticUnits{d1,T1}, D2, T2, U2<:Units{D2,T2}} = Units{promote_type(typeof(d1),D2), promote_type(T1,T2)}
Base.promote_rule(::Type{U1}, ::Type{U2}) where {T1, d1, U1<:StaticUnits{d1,T1}, d2, T2, U2<:StaticUnits{d2,T2}} = Units{promote_type(typeof(d1),typeof(d2)), promote_type(T1,T2)}

#Promote static unit quantities to static dimenion quantities, double definition needed for specificity
function Base.promote_rule(::Type{Quantity{T1,U1}}, ::Type{Quantity{T2,U2}}) where {T1, T2, U1<:StaticUnits, U2<:StaticDims}
    D = equaldims(dimval(U1), dimval(U2))
    T = promote_type(T1, T2)
    return Quantity{T, StaticDims{D}}
end
function Base.promote_rule(::Type{Quantity{T1,U1}}, ::Type{Quantity{T2,U2}}) where {T1, T2, U1<:StaticDims, U2<:StaticUnits}
    D = equaldims(dimval(U1), dimval(U2))
    T = promote_type(T1, T2)
    return Quantity{T, StaticDims{D}}
end

#Unit promotion
function Base.promote_rule(::Type{D1}, ::Type{Units{D2,T2}}) where {D1<:AbstractDimensions, D2<:AbstractDimensions, T2<:AbstractUnitTransform}
    return Units{promote_type(D1, D2), T2}
end
function Base.promote_rule(::Type{Units{D1,T1}}, ::Type{Units{D2,T2}}) where {D1<:AbstractDimensions, D2<:AbstractDimensions, T1<:AbstractUnitTransform, T2<:AbstractUnitTransform}
    return Units{promote_type(D1, D2), promote_type(T1,T2)}
end
function Base.promote_rule(::Type{StaticDims{D}}, ::Type{StaticUnits{D,T}}) where {D, T<:AbstractUnitTransform}
    return StaticUnits{D,T}
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
Base.promote_rule(::Type{Quantity{T,U}}, ::Type{T}) where {T,U<:AbstractUnitLike} = Quantity{T,U}
Base.promote_rule(::Type{Quantity{T1,U}}, ::Type{T2}) where {T1<:Number, T2<:Number, U<:AbstractUnitLike} = Quantity{promote_type(T1,T2), U}
Base.promote_rule(::Type{Quantity{T1,U}}, ::Type{T2}) where {T1<:AbstractArray, T2<:AbstractArray, U<:AbstractUnitLike} = Quantity{promote_type(T1,T2), U}
