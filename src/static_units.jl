include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")

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

struct StaticUnits{D, C<:AbstractUnitConverter} <: AbstractUnitLike
    todims :: C
    symbol :: Symbol
    function StaticUnits{D}(conv::C, symb=DEFAULT_USYMBOL::Symbol) where {D, C<:AbstractUnitConverter}
        return (D isa AbstractDimensions) ? new{D,C}(conv, symb) : error("Type parameter must be a dimension")
    end
end
StaticUnits(u::AffineUnits) = StaticUnits{dimension(u)}(uconvert(dimension(u), u), u.symbol)
AffineUnits(u::StaticUnits) = AffineUnits{dimtype(u)}(scale=u.todims.scale, offset=u.todims.offset, dims=dimval(u), symbol=u.symbol)
udynamic(u::StaticUnits{D, U}) where {D, U<:AffineConverter} = AffineUnits(u)

dimtype(::Type{<:StaticUnits{D,C}}) where {D,C} = typeof(D)
dimtype(d::StaticUnits) = dimtype(typeof(d))
dimval(::Type{<:StaticUnits{D,C}}) where {D,C} = D
dimval(d::StaticUnits) = dimval(typeof(d))

dimtype(q::Quantity) = dimtype(unit(q))

#============================================================================================================================
Conversions
============================================================================================================================#
uconvert(utarget::StaticDims{D}, ucurrent::StaticUnits{D}) where D = ucurrent.todims

function ubase(q::AbstractQuantity{T, <:StaticUnits{D}}) where {T,D}
    x = unit(q).todims(ustrip(q))
    return Quantity{typeof(x), StaticDims{D}}(x, StaticDims{D}())
end

Quantity(x::T, u::StaticUnits{D}) where {T,D} = Quantity(u.todims(x), StaticDims{D}())

Base.convert(::Type{Q}, q::AbstractQuantity{<:Any, <:StaticUnits{D}}) where {T,D,Q<:Quantity{T,StaticDims{D}}} = Q(dstrip(q), StaticDims{D}())
Base.convert(::Type{Q}, q::AbstractQuantity{<:Any, <:StaticDims{D}}) where {T,D,Q<:Quantity{T,StaticDims{D}}}  = Q(dstrip(q), StaticDims{D}())
Base.convert(::Type{Q}, q::AbstractQuantity{<:Any, <:StaticDims{D}}) where {T,D,C,Q<:Quantity{T,StaticUnits{D,C}}} = Q(ustrip(q), StaticUnits{D}(C()))

#Conflicting or uncertain dimensions get promoted to "Dimension"
Base.convert(::Type{D}, d::StaticDims) where {D<:AbstractDimensions} = convert(D, dimval(d))
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:AbstractDimensions, D2<:StaticDims} = promote_type(D1, dimtype(D2))
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:StaticDims, D2<:StaticDims} = promote_type(dimtype(D1), dimtype(D2))

#Promote static units to static dimenions
function Base.promote_rule(::Type{Quantity{T1,StaticUnits{D,C}}}, ::Type{Quantity{T2,StaticDims{D}}}) where {T1, T2, D, C}
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
equaldims(arg1::StaticDims, arg2::StaticDims) = (dimval(arg1) == dimal(arg2))  ? arg1 : throw(DimensionError((arg1,arg2)))

Base.:*(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1*D2}()
Base.:/(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1/D2}()
Base.inv(arg::StaticDims{D}) where D = StaticDims{inv(D)}()
Base.:^(d::StaticDims{D}, p::Real) where D = StaticDims{D^p}()
Base.sqrt(d::StaticDims{D}) where D = StaticDims{sqrt(D)}()
Base.cbrt(d::StaticDims{D}) where D = StaticDims{sqrt(D)}()
Base.abs2(d::StaticDims{D}) where D = StaticDims{abs2(D)}()
Base.adjoint(d::StaticDims{D}) where D = StaticDims{adjoint(D)}()

@inline Base.literal_pow(::typeof(^), d::D, ::Val{0}) where {D <: AbstractDimensions} = D()
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{1}) = d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{2}) = d*d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{3}) = d*d*d
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-1}) = inv(d) 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-2}) = inv(d*d)



#Tests
using Test
import .UnitRegistry.@u_str

@testset "Static Dimensions" begin
    @test 1*StaticUnits(u"km/hr") isa Quantity{<:Any, <:StaticDims}

    q1 = Quantity{Float64, StaticUnits{dimension(u"km/hr"), AffineConverter}}(5, StaticUnits(u"km/hr"))
    q2 = 20.0*StaticUnits(u"m/s")
    [q1,q2]

end