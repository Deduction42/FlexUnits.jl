#===============================================================================================
 FixedRational was heavily inspired by DynamicQuantities.FixedRational
 https://github.com/SymbolicML/DynamicQuantities.jl?tab=Apache-2.0-1-ov-file
===============================================================================================#
"""
    Numerator{T<:Integer}

Internal struct that indicates result should be interpreated as a numerator
"""
struct Numerator{T<:Integer}
    val :: T
end


"""
    FixedRational{B,T<:Integer} <: Real

Rational number type, with fixed base (or denominator) of B
Faster and simpler than the base rational type, as operations
can usually convert to integer ops, but precision is limited to 1/B.
"""
struct FixedRational{B,T<:Integer} <: Real
    num :: T
    FixedRational{B}(x::Numerator{T}) where {B,T} = new{B,T}(x.val)
    FixedRational{B,T}(x::Numerator) where {B,T} = new{B,T}(x.val)
    FixedRational{B}(x::T) where {B,T<:Integer} = new{B,T}(x*B)
    FixedRational{B,T}(x::Integer) where {B,T} = new{B,T}(x*B)
    FixedRational{B,T}(x::Real) where {B,T} = new{B,T}(round(T, x*B))
end

#Default Rational Base, common factors of prime numbers multiplying to less than sqrt(typemax(Int32))
const DEFAULT_NUMERATOR_TYPE = Int32
const DEFAULT_DENOM = 2^4 * 3^2 * 5^2 * 7
FixedRational(x::Number)  = FixedRational{DEFAULT_DENOM, DEFAULT_NUMERATOR_TYPE}(x)
FixedRational(x::Integer) = FixedRational{DEFAULT_DENOM}(x)

#Julia's Rational API
Base.numerator(x::FixedRational) = x.num 
Base.denominator(x::FixedRational{B,T}) where {B,T} = convert(T, B)
numtype(::Type{<:FixedRational{B,T}}) where {B,T} = T
Numerator(x::FixedRational) = Numerator(x.num)

#Conversion to float and rational
function (::Type{T})(x::FixedRational) where {T<:Union{AbstractFloat,Integer}}
    return convert(T, numerator(x)/denominator(x))
end
Base.Bool(x::FixedRational) = iszero(x) ? false : isone(x) ? true : throw(InexactError(:Bool, Bool, x))
Base.Rational(x::FixedRational) = numerator(x)//denominator(x)

#Promotion rule:
promote_rule(::Type{FixedRational{B,T}}, ::Type{S}) where {B, T<:Integer, S<:AbstractFloat} = promote_type(T,S)

#Mathematical operations on similar objects
Base.:+(x1::FixedRational{B}, x2::FixedRational{B}) where {B} = FixedRational{B}(Numerator(x1.num + x2.num))
Base.:-(x::FixedRational{B}) where {B} = FixedRational{B}(Numerator(-x.num))
Base.:-(x1::FixedRational{B}, x2::FixedRational{B}) where {B} = FixedRational{B}(Numerator(x1.num - x2.num))
Base.:*(x1::FixedRational{B,T1}, x2::FixedRational{B,T2}) where {B,T1,T2} = FixedRational{B,promote_type(T1,T2)}(Numerator(widemul(x1.num, x2.num) ÷ B))
Base.:/(x1::FixedRational{B,T1}, x2::FixedRational{B,T2}) where {B,T1,T2} = FixedRational{B,promote_type(T1,T2)}(Numerator(widemul(B, x1.num) ÷ x2.num))
Base.:inv(x::FixedRational{B,T}) where {B,T} = FixedRational{B,T}(Numerator(widemul(B, B) ÷ x.num))

#Specialized mathematical operations on different number types that don't require promotion
Base.:*(xi::Integer, xr::FixedRational{B,T}) where {B,T} = FixedRational{B,T}(Numerator(xi*xr.num))
Base.:*(xr::FixedRational{B,T}, xi::Integer) where {B,T} = FixedRational{B,T}(Numerator(xi*xr.num))

#Rounding
Base.round(::Type{T}, x::FixedRational, r::RoundingMode=RoundNearest) where {T} = div(convert(T, numerator(x)), convert(T, denominator(x)), r)

#Comparisons
for comp in (:(==), :isequal, :<, :(isless), :<=)
    @eval Base.$comp(x::FixedRational{B}, y::FixedRational{B}) where {B} = $comp(x.num, y.num)
end

#Check values
Base.iszero(x::FixedRational) = iszero(numerator(x))
Base.isone(x::FixedRational)  = (numerator(x) == denominator(x))
Base.isinteger(x::FixedRational) = iszero(numerator(x) % denominator(x))

function Base.convert(::Type{F}, x::FixedRational{B}) where {B, T, F<:FixedRational{B,T}}
    return FixedRational{B,T}(Numerator(x.num))
end

# Promotion rules with self and other types
function Base.promote_rule(::Type{FixedRational{B1,T1}}, ::Type{FixedRational{B2,T2}}) where {B1,B2,T1,T2}
        B1 == B2 || error("Refusing to promote `FixedRational` types with mixed denominators. Use `Rational` instead.")
    return FixedRational{B1, promote_type(T1,T2)}
end
function Base.promote_rule(::Type{FixedRational{B,T1}}, ::Type{Rational{T2}}) where {B,T1,T2}
    return Rational{promote_type(T1, T2)}
end
function Base.promote_rule(::Type{F}, ::Type{<:Integer}) where {F<:FixedRational}
    return F
end
function Base.promote_rule(::Type{FixedRational{B,T1}}, ::Type{T2}) where {B,T1,T2<:Real}
    return promote_type(Rational{T1}, T2)
end

#Ambiguities
function Base.promote_rule(::Type{F}, ::Type{Bool}) where {F<:FixedRational}
    return F
end
function Base.promote_rule(::Type{FixedRational{B,T}}, ::Type{BigFloat}) where {B,T}
    return promote_type(Rational{T}, BigFloat)
end
function Base.promote_rule(::Type{FixedRational{B,T}}, ::Type{T2}) where {B,T,T2<:AbstractIrrational}
    return promote_type(Rational{T}, T2)
end

#Printing/showing
function Base.string(x::FixedRational{B,T}) where {B,T}
    if isinteger(x)
        return string(convert(T, x))
    end
    g = gcd(x.num, B)
    return string(div(x.num, g)) * "//" * string(div(B, g))
end
Base.show(io::IO, x::FixedRational) = print(io, string(x))