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
dimtype(::Type{StaticDims{D}}) where D = typeof(D)
dimval(::Type{StaticDims{D}}) where D = D
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

dimtype(::Type{StaticUnits{D,C}}) where {D,C} = typeof(D)
dimval(::Type{StaticUnits{D,C}}) where {D,C} = D

#============================================================================================================================
Conversions
============================================================================================================================#
uconvert(utarget::StaticDims{D}, ucurrent::StaticUnits{D}) where D = ucurrent.todims

function ubase(q::AbstractQuantity{T, StaticUnits{D}}) where {T,D}
    x = unit(q).todims(ustrip(q))
    return Quantity{typeof(x), StaticDims{D}}(x, StaticDims{D}())
end

Quantity(x::T, u::StaticUnits{D}) where {T,D} = Quantity(u.todims(x), StaticDims{D}())


Base.convert(::Type{D}, d::StaticDims) where {D<:AbstractDimensions} = convert(D, dimval(d))
Base.convert(::Type{Quantity{T,D}}, q::Quantity{<:Any, StaticDims{d}}) where {T,d,D<:AbstractDimensions} = Quantity{T,D}(ustrip(q), d)

Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:AbstractDimensions, D2<:StaticDims} = promote_type(D1, dimtype(D2))
Base.promote_rule(::Type{D1}, ::Type{D2}) where {D1<:StaticDims, D2<:StaticDims} = promote_type(dimtype(D1), dimtype(D2))

#============================================================================================================================
Utilities
============================================================================================================================#
function Base.show(io::IO, u::StaticDims; pretty=PRETTY_DIM_OUTPUT[])
    return Base.show(io, udynamic(u), pretty=pretty)
end

function Base.show(io::IO, u::StaticUnits; pretty=PRETTY_DIM_OUTPUT[])
    return Base.show(ui, udynamic(u), pretty=pretty)
end

#============================================================================================================================
Promotions
============================================================================================================================#




#Tests
import .UnitRegistry.@u_str
q = 1*StaticUnits(u"km/hr")