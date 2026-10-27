module DatesExt 

import Dates: Nanosecond, Microsecond, Millisecond, Second, Minute, Hour, Day, Week, FixedPeriod, CompoundPeriod
import FlexUnits: seconds_unit, ustrip, tobase, usymbol, quantity, ubase, uscale, dimtype, dimension, uconvert
import FlexUnits: AbstractDimensions, AffineTransform, AbstractDimLike, StaticDims, Quantity, Units, AbstractUnitLike

#==============================================================================================================================
Constructors and conversions
==============================================================================================================================#
quantity(::Type{D}, p::FixedPeriod) where {D<:AbstractDimLike} = quantity(p.value, Units{D}(typeof(p)))
quantity(::Type{D}, p::CompoundPeriod) where {D<:AbstractDimLike} = ubase(D, p)
ubase(::Type{D}, p::Union{FixedPeriod,CompoundPeriod}) where {D<:AbstractDimLike} = ubase(toseconds(p), seconds_unit(D))

uconvert(u::Units{D}, p::FixedPeriod) where {D<:AbstractDimLike} = uconvert(u, ubase(D, p))
uconvert(u::D, p::FixedPeriod) where {D<:AbstractDimLike} = uconvert(u, ubase(D, p))
uconvert(::Type{P}, q::Quantity) where {P<:FixedPeriod} = uconvert(Units{dimtype(q)}(p), q)

Quantity{T,U}(p::FixedPeriod) where {T,U} = Quantity{T,U}(p.value, Units{dimtype(U)}(p))
Base.convert(::Type{Quantity{T,U}}, p::Union{FixedPeriod,CompoundPeriod}) where {T,U} = convert(Quantity{T,U}, quantity(dimtype(U), p))
Base.convert(::Type{P}, q::Quantity) where {P<:FixedPeriod} = P(ustrip(Units{dimtype(q)}(P), q))

function Units{D}(::Type{P}) where {D<:AbstractDimLike, P<:FixedPeriod}
    u = seconds_unit(D)
    return Units{D}(dimension(u), AffineTransform(uscale(u)*toseconds(P), 0), usymbol(P))
end
Units{D}(p::P) where {D<:AbstractDimLike, P<:FixedPeriod} = Units{D}(P)

Base.promote_rule(::Type{Quantity{T,U}}, ::Type{<:FixedPeriod}) where {T,U} = Quantity{promote_type(T,Int),U}
Base.promote_rule(::Type{Quantity{T,U}}, ::Type{<:CompoundPeriod}) where {T,U} = Quantity{promote_type(T,Int),U}

#==============================================================================================================================
Common operators
==============================================================================================================================#
for op = (:+, :-, :*, :/, ://, :(==), :isequal, :<, :isless, :≤)
    @eval Base.$op(p::Union{FixedPeriod,CompoundPeriod}, q::Quantity) = $op(ubase(dimtype(q), p), q)
    @eval Base.$op(q::Quantity, p::Union{FixedPeriod,CompoundPeriod}) = $op(q, ubase(dimtype(q), p))
end

for op = (:*, :/, ://)
    @eval Base.$op(p::Union{FixedPeriod,CompoundPeriod}, u::AbstractUnitLike) = $op(ubase(dimtype(u), p), u)
    @eval Base.$op(u::AbstractUnitLike, p::Union{FixedPeriod,CompoundPeriod}) = $op(u, ubase(dimtype(u), p))
end

Base.isapprox(p::Union{FixedPeriod,CompoundPeriod}, q::Quantity; kwargs...) = isapprox(ubase(dimtype(q), p), q, kwargs...)
Base.isapprox(q::Quantity, p::Union{FixedPeriod,CompoundPeriod}; kwargs...) = isapprox(q, ubase(dimtype(q), p), kwargs...)

#==============================================================================================================================
Rounding
==============================================================================================================================#
Base.round(::Type{P}, q::Quantity, r::RoundingMode=RoundNearest; kwargs...) where {P<:FixedPeriod} = P(round(Units{dimtype(q)}(P), q, r), kwargs...)
Base.round(u::AbstractUnitLike, p::FixedPeriod, r::RoundingMode=RoundNearest; kwargs...) = round(u, ubase(dimtype(u), p), r, kwargs...)

Base.trunc(u::AbstractUnitLike, p::FixedPeriod; kwargs...) = round(u, p, RoundToZero; kwargs...)
Base.floor(u::AbstractUnitLike, p::FixedPeriod; kwargs...) = round(u, p, RoundDown; kwargs...)
Base.ceil(u::AbstractUnitLike, p::FixedPeriod; kwargs...)  = round(u, p, RoundUp; kwargs...)

Base.trunc(::Type{P}, q::Quantity; kwargs...) where {P<:FixedPeriod} = round(P, q, RoundToZero; kwargs...)
Base.floor(::Type{P}, q::Quantity; kwargs...) where {P<:FixedPeriod} = round(P, q, RoundDown; kwargs...)
Base.ceil(::Type{P}, q::Quantity; kwargs...) where {P<:FixedPeriod} = round(P, q, RoundUp; kwargs...)

#==============================================================================================================================
Utility functions
==============================================================================================================================#
# Convert all periods to seconds
toseconds(::Type{Nanosecond}) = 1e-9
toseconds(::Type{Microsecond}) = 1e-6
toseconds(::Type{Millisecond}) = 1e-3
toseconds(::Type{Second}) = true
toseconds(::Type{Minute}) = 60
toseconds(::Type{Hour}) = 3600
toseconds(::Type{Day}) = 86400
toseconds(::Type{Week}) = 604800
toseconds(p::FixedPeriod) = toseconds(typeof(p))*p.value
toseconds(p::CompoundPeriod) = sum(toseconds, p.periods, init=0.0)
toseconds(p::Second) = p.value

# Extract units
usymbol(::Type{Nanosecond}) = :ns
usymbol(::Type{Microsecond}) = :μs
usymbol(::Type{Millisecond}) = :ms 
usymbol(::Type{Second}) = :s 
usymbol(::Type{Minute}) = :min
usymbol(::Type{Hour}) = :hr
usymbol(::Type{Day}) = :day 
usymbol(::Type{Week}) = :wk 



end