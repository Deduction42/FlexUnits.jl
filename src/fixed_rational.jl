#===============================================================================================
 FixedRational was heavily inspired by DynamicQuantities.FixedRational
 https://github.com/SymbolicML/DynamicQuantities.jl?tab=Apache-2.0-1-ov-file
===============================================================================================#
"""
    Numerator{T<:Signed}

Internal struct that indicates the value should be interpreated as a numerator, mainly intended 
for dispatching FixedRational constructors to directly assign a numerator.
"""
struct Numerator{T<:Signed}
    val :: T
end


"""
    FixedRational{B,T<:Signed} <: Real

Rational number type, with fixed base (or denominator) of B. FixedRational is faster and simpler 
than Julia's base Rational type, as operations on FixedRational usually result in simple integer, 
operations; however, precision of FixedRational limited to 1/B.

When constructing, FixedRational{B,T}(x::Real) uses integer division to assign a numerator that yields 
the closest value to "x". To assign a numerator directly, use FixedRational{B,T}(Numerator(x)).
"""
struct FixedRational{B,T<:Signed} <: Real
    num :: T
    FixedRational{B}(x::Numerator{T}) where {B,T} = new{B,T}(x.val)
    FixedRational{B,T}(x::Numerator) where {B,T} = new{B,T}(x.val)
    FixedRational{B}(x::T) where {B,T<:Signed} = new{B,T}(x*B)
    FixedRational{B,T}(x::Signed) where {B,T} = new{B,T}(x*B)
    FixedRational{B,T}(x::Real) where {B,T} = new{B,T}(round(T, x*B))
end

#Default Rational Base, common factors of prime numbers multiplying to less than sqrt(typemax(Int32))
const DEFAULT_NUMERATOR_TYPE = Int32
const DEFAULT_DENOM = 2^4 * 3^2 * 5^2 * 7
const DEFAULT_DENOM_64 = 2^7 * 3^4 * 5^4 * 7^2 * 11
const FixRat32 = FixedRational{DEFAULT_DENOM, Int32}
const FixRat64 = FixedRational{DEFAULT_DENOM_64, Int64}

FixedRational(x::FixedRational) = x
FixedRational(x::Real) = FixedRational{DEFAULT_DENOM, DEFAULT_NUMERATOR_TYPE}(x)
FixedRational(x::Signed) = FixedRational{DEFAULT_DENOM}(x)
Numerator(x::FixedRational) = Numerator(x.num)

#Julia's Rational API
Base.numerator(x::FixedRational) = x.num 
Base.denominator(x::FixedRational{B,T}) where {B,T} = convert(T, B)

#Conversion to float and rational
function (::Type{T})(x::FixedRational) where {T<:Union{AbstractFloat,Signed}}
    return convert(T, numerator(x)/denominator(x))
end
Base.Bool(x::FixedRational) = iszero(x) ? false : isone(x) ? true : throw(InexactError(:Bool, Bool, x))
Base.Rational{T}(x::FixedRational) where T = Rational{T}(numerator(x)//denominator(x))
Base.Rational(x::FixedRational) = numerator(x)//denominator(x)

#Mathematical operations on similar objects
Base.:+(x1::FixedRational{B}, x2::FixedRational{B}) where {B} = FixedRational{B}(Numerator(x1.num + x2.num))
Base.:-(x::FixedRational{B}) where {B} = FixedRational{B}(Numerator(-x.num))
Base.:-(x1::FixedRational{B}, x2::FixedRational{B}) where {B} = FixedRational{B}(Numerator(x1.num - x2.num))
Base.:*(x1::FixedRational{B,T1}, x2::FixedRational{B,T2}) where {B,T1,T2} = FixedRational{B,promote_type(T1,T2)}(Numerator(widemul(x1.num, x2.num) ÷ B))
Base.:/(x1::FixedRational{B,T1}, x2::FixedRational{B,T2}) where {B,T1,T2} = FixedRational{B,promote_type(T1,T2)}(Numerator(widemul(B, x1.num) ÷ x2.num))
Base.:inv(x::FixedRational{B,T}) where {B,T} = FixedRational{B,T}(Numerator(widemul(B, B) ÷ x.num))

#Specialized mathematical operations on different number types that don't require promotion
Base.:*(xi::Signed, xr::FixedRational{B,T}) where {B,T} = FixedRational{B,T}(Numerator(xi*xr.num))
Base.:*(xr::FixedRational{B,T}, xi::Signed) where {B,T} = FixedRational{B,T}(Numerator(xi*xr.num))

#Rounding
Base.round(::Type{T}, x::FixedRational, r::RoundingMode=RoundNearest) where {T} = div(convert(T, numerator(x)), convert(T, denominator(x)), r)
Base.round(::Type{>:Missing}, x::FixedRational, r::RoundingMode=RoundNearest) = missing

#Comparisons
for comp in (:(==), :isequal, :<, :(isless), :<=)
    @eval Base.$comp(x::FixedRational{B}, y::FixedRational{B}) where {B} = $comp(x.num, y.num)
end

#Check values
Base.iszero(x::FixedRational) = iszero(numerator(x))
Base.isone(x::FixedRational)  = (numerator(x) == denominator(x))
Base.isinteger(x::FixedRational) = iszero(numerator(x) % denominator(x))

#Fixed rational conversions
(::Type{F})(x::FixedRational{B}) where {B,F<:FixedRational{B}} = F(Numerator(x)) #Same-base shortcut
(::Type{F})(x::FixedRational) where {F<:FixedRational} = F(numerator(x)/denominator(x))

#Promotion rules:
function Base.promote_rule(::Type{FixedRational{B1,T1}}, ::Type{FixedRational{B2,T2}}) where {B1,B2,T1,T2}
    B1 == B2 || error("Refusing to promote `FixedRational` types with mixed denominators. Use `Rational` instead.")
return FixedRational{B1, promote_type(T1,T2)}
end
function Base.promote_rule(::Type{FixedRational{B,I}}, ::Type{T}) where {B,I,T<:Real}
    return promote_type(Rational{I}, T)
end
function Base.promote_rule(::Type{FixedRational{B,I}}, ::Type{Rational{T}}) where {B,I,T}
    return Rational{promote_type(I,T)}
end
function Base.promote_rule(::Type{F}, ::Type{<:Signed}) where {F<:FixedRational}
    return F
end

#Printing/showing
function Base.show(io::IO, x::FixedRational{B,T}) where {B,T}
    if isinteger(x)
        return print(io, convert(T, x))
    end
    g = gcd(x.num, B)
    print(io, div(x.num, g))
    print(io, "//")
    return print(io, div(B, g))
end

function Base.show(io::IO, ::Type{R}) where {B,T,R<:FixedRational{B,T}}
    return print(io, replace("FixRat$(T)", "Int"=>""))
end