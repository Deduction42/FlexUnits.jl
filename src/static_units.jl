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
dimtype(::Type{StaticDims{D}}) where D= typeof(D)

struct StaticUnit{D, C<:AbstractUnitConverter} <: AbstractUnitLike
    todims :: C
    symbol :: Symbol
    function StaticUnit{D}(conv::C, symb=DEFAULT_USYMBOL::Symbol) where {D, C<:AbstractUnitConverter}
        return (D isa AbstractDimensions) ? new{D,C}(conv, symb) : error("Type parameter must be a dimension")
    end
end
StaticUnit(u::AffineUnits) = StaticUnit{dimension(u)}(uconvert(dimension(u), u), u.symbol)
dimtype(::Type{<:StaticUnit{D}}) where D = typeof(D)

#============================================================================================================================
Conversions
============================================================================================================================#
uconvert(utarget::StaticDims{D}, ucurrent::StaticUnit{D}) where D = ucurrent.todims

function ubase(q::AbstractQuantity{T, StaticUnit{D}}) where {T,D}
    x = unit(q).todims(ustrip(q))
    return Quantity{typeof(x), StaticDims{D}}(x, StaticDims{D}())
end

Quantity(x::T, u::StaticUnit{D}) where {T,D} = Quantity(u.todims(x), StaticDims{D}())


convert(::Type{D}, d::StaticDims{D0}) where {D0, D<:AbstractDimensions} = convert(D, D0)
convert(::Type{Quantity{T,D}}, q::Quantity{<:Any, StaticDims{d}}) where {T,d,D<:AbstractDimensions} = Quantity{T,D}(ustrip(q), d)

promote_rule(::D1, ::D2) where {D1<:AbstractDimensions, D2<:StaticDims} = D1
promote_rule(::D1, ::D2) where {D1<:StaticDims, D2<:StaticDims} = promote_type(dimtype(D1), dimtype(D2))

#============================================================================================================================
Utilities
============================================================================================================================#
function Base.show(io::IO, u::StaticDims{D}; pretty=PRETTY_DIM_OUTPUT[]) where D
    return Base.show(io, D, pretty=pretty)
end

function Base.show(io::IO, u::StaticUnit{D}; pretty=PRETTY_DIM_OUTPUT[]) where D
    usymb = u.symbol
    if usymb == DEFAULT_USYMBOL
        return Base.show(io, D, pretty=pretty)
    else
        return Base.print(io, usymb)
    end
end

#============================================================================================================================
Promotions
============================================================================================================================#




#Tests
import .UnitRegistry.@u_str
q = 1*StaticUnit(u"km/hr")