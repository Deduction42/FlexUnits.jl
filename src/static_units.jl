#=
include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")
=#

#============================================================================================================================
Type definitions
============================================================================================================================#
struct StaticDims{D} <: AbstractDimLike
    function StaticDims{D}() where D
        return (D isa AbstractDimensions) ? new{D}() : error("Type parameter must be a dimension")
    end
end 
StaticDims(D::AbstractDimensions) = StaticDims{D}()
StaticDims{D}(d::AbstractDimensions) where D = (D == d) ? StaticDims{D} : throw(ArgumentError("Dimesion $(d) must be equal to $(D)"))
dimtype(::Type{<:StaticDims{D}}) where D = typeof(D)
dimtype(d::StaticDims) = dimtype(typeof(d))
dimval(::Type{<:StaticDims{D}}) where D = D
dimval(d::StaticDims) = dimval(typeof(d))
udynamic(u::StaticDims{D}) where D = D

struct StaticUnits{D, C<:AbstractUnitTransform} <: AbstractUnitLike
    todims :: C
    symbol :: Symbol
    function StaticUnits{D}(conv::C, symb=DEFAULT_USYMBOL::Symbol) where {D, C<:AbstractUnitTransform}
        return (D isa AbstractDimensions) ? new{D,C}(conv, symb) : error("Type parameter must be a dimension")
    end
end
StaticUnits(u::AffineUnits) = StaticUnits{dimension(u)}(uconvert(dimension(u), u), u.symbol)
AffineUnits(u::StaticUnits) = AffineUnits{dimtype(u)}(scale=u.todims.scale, offset=u.todims.offset, dims=dimval(u), symbol=u.symbol)
todims(u::StaticUnits) = u.todims
udynamic(u::StaticUnits{D, U}) where {D, U<:AffineTransform} = AffineUnits(u)
dimtype(::Type{StaticUnits{D,C}}) where {D,C} = typeof(D)
dimtype(d::StaticUnits) = dimtype(typeof(d))
dimval(::Type{StaticUnits{D,C}}) where {D,C} = D
dimval(d::StaticUnits) = dimval(typeof(d))
dimension(::Type{StaticUnits{D,T}}) where {D,T} = StaticDims{D}()
dimension(d::StaticUnits) = dimension(typeof(d))



Base.:(==)(d1::AbstractDimensions, d2::StaticDims) = (d1 == dimval(d2))
Base.:(==)(d1::StaticDims, d2::AbstractDimensions) = (dimval(d1) == d2)

#============================================================================================================================
Conversions
============================================================================================================================#
uconvert(utarget::StaticDims{D}, ucurrent::StaticUnits{D}) where D = ucurrent.todims

function ubase(q::AbstractQuantity{T, <:StaticUnits{D}}) where {T,D}
    x = unit(q).todims(ustrip(q))
    return Quantity{typeof(x), StaticDims{D}}(x, StaticDims{D}())
end

Quantity(x::T, u::StaticUnits{D}) where {T,D} = Quantity(u.todims(x), StaticDims{D}())
Quantity{T}(x, u::StaticUnits{D}) where {T,D} = Quantity(convert(T, u.todims(x)), StaticDims{D}())

Base.convert(::Type{Q}, q::AbstractQuantity{<:Any, <:StaticUnits{D}}) where {T,D,Q<:Quantity{T,StaticDims{D}}} = Q(dstrip(q), StaticDims{D}())
Base.convert(::Type{Q}, q::AbstractQuantity{<:Any, <:StaticDims{D}}) where {T,D,Q<:Quantity{T,StaticDims{D}}}  = Q(dstrip(q), StaticDims{D}())
Base.convert(::Type{Q}, q::AbstractQuantity{<:Any, <:StaticDims{D}}) where {T,D,C,Q<:Quantity{T,StaticUnits{D,C}}} = Q(ustrip(q), StaticUnits{D}(C()))

#Conflicting or uncertain dimensions get promoted to "Dimension"
Base.convert(::Type{D}, d::StaticDims) where {D<:AbstractDimensions} = convert(D, dimval(d))
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:AbstractDimensions, D2<:StaticDims} = promote_type(D1, dimtype(D2))
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:StaticDims, D2<:StaticDims} = promote_type(dimtype(D1), dimtype(D2))

#Promote static units to static dimenions, double definition needed for specificity
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

#============================================================================================================================
Utilities
============================================================================================================================#
function Base.show(io::IO, u::StaticDims; pretty=PRETTY_DIM_OUTPUT[])
    return Base.show(io, udynamic(u), pretty=pretty)
end

function Base.show(io::IO, u::StaticUnits; pretty=PRETTY_DIM_OUTPUT[])
    return Base.show(io, udynamic(u), pretty=pretty)
end

#============================================================================================================================
Static dimension ops
============================================================================================================================#
equaldims(arg1::StaticDims, arg2::AbstractDimensions) = (dimval(arg1) == arg2) ? arg1 : throw(DimensionError((arg1,arg2)))
equaldims(arg1::AbstractDimensions, arg2::StaticDims) = (dimval(arg1) == arg2) ? arg2 : throw(DimensionError((arg1,arg2)))
equaldims(arg1::StaticDims, arg2::StaticDims) = (dimval(arg1) == dimval(arg2))  ? arg1 : throw(DimensionError((arg1,arg2)))

Base.:*(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1*D2}()
Base.:/(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1/D2}()
Base.inv(arg::StaticDims{D}) where D = StaticDims{inv(D)}()
Base.:^(d::StaticDims{D}, p::Real) where D = StaticDims{D^p}()
Base.sqrt(d::StaticDims{D}) where D = StaticDims{sqrt(D)}()
Base.cbrt(d::StaticDims{D}) where D = StaticDims{sqrt(D)}()
Base.abs2(d::StaticDims{D}) where D = StaticDims{abs2(D)}()
Base.adjoint(d::StaticDims{D}) where D = StaticDims{adjoint(D)}()

@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{0}) = StaticDims{dimtype(d)()}
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{1}) = d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{2}) = d*d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{3}) = d*d*d
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-1}) = inv(d) 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-2}) = inv(d*d)


#=
#Tests
using Test
using BenchmarkTools
import .UnitRegistry.@u_str


@testset "Static Dimensions" begin
    @test 1*StaticUnits(u"km/hr") isa Quantity{<:Any, <:StaticDims}

    #Test promotion rules
    q1 = Quantity{Float64, StaticUnits{dimension(u"km/hr"), AffineTransform}}(5, StaticUnits(u"km/hr"))
    q2 = 20.0*StaticUnits(u"m/s")
    q3 = 10*StaticUnits(u"kg/s")

    @test eltype([q1,q1]) <: Quantity{Float64, <:StaticUnits{dimension(u"m/s")}} #StaticUnits preserved
    @test eltype([q1,q2]) <: Quantity{Float64, <:StaticDims{dimension(u"m/s")}} #StaticUnits and StaticDims promote to StaticDims if dimension is the same
    @test eltype([q1,q2,q3]) <: Quantity{Float64, <:Dimensions} #Different dimensions promote to "Dimensions"

end



#=
# This is an example where static FlexUnits VASTLY outperforms Unitful due to the design choice of sticking with dimensions
import Unitful
import DynamicQuantities

v1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s"]
v1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s"]
v1flex = [1.0*StaticUnits(u"m/s"), 1.0*StaticUnits(u"m/s"), 1.0*StaticUnits(u"m/s")]

@btime sum(x->x^2.0, v1uni)
@btime sum(x->x^2.0, $v1dyn)
@btime sum(x->x^2.0, $v1flex)
=#
=#