using Test
using FlexUnits

#===============================================================================================================================
Define your custom dimensions
===============================================================================================================================#
@kwdef struct MoneyDimensions{P} <: AbstractDimensions{P}
    m   ::P = zero(FixRat32)
    kg  ::P = zero(FixRat32)
    s   ::P = zero(FixRat32)
    A   ::P = zero(FixRat32)
    K   ::P = zero(FixRat32)
    cd  ::P = zero(FixRat32)
    mol ::P = zero(FixRat32)
    €   ::P = zero(FixRat32)
end
MoneyDimensions(args::Real...) = MoneyDimensions{FixRat32}(args...)

#===============================================================================================================================
Create your unit registry as a module using RegistryTools
===============================================================================================================================#
module CurrencyUnits

using FlexUnits
using FlexUnits.RegistryTools
import ..MoneyDimensions

const UNITS = PermanentDict{Symbol,Units{MoneyDimensions{FixRat32},AffineTransform{Float64}}}()

registry_defaults!(UNITS)
register_unit!(UNITS, "€" => MoneyDimensions{FixRat32}(€=1))
register_unit!(UNITS, "EUR" => UNITS[:€])

uparse(str::String) = RegistryTools.uparse(str, UNITS)
qparse(str::String) = RegistryTools.qparse(str, UNITS)

macro u_str(str); return suparse_expr(str, UNITS); end
macro ud_str(str); return uparse_expr(str, UNITS); end
macro q_str(str); return qparse_expr(str, UNITS); end
macro U_str(str); suexpr = suparse_expr(str, UNITS); return :($typeof($suexpr)); end
macro D_str(str); suexpr = suparse_expr(str, UNITS); return :($dimtype($suexpr)); end

utype() = RegistryTools.regunittype(UNITS)
dtype() = RegistryTools.regdimtype(UNITS)

export @u_str, @ud_str, @q_str, @U_str, @D_str, uparse, qparse, utype, dtype

end


#===============================================================================================================================
If you want simplify(...) to work with your units, 
    => define a constant vector of preferred units that are compatible with your new type
    => overload FlexUnits.preferred_units(::Type{T}) to refer to that list of preferred units
===============================================================================================================================#
const CURRENCY_UNITS = [ CurrencyUnits.UNITS[k] for k in [:F, :H, :T, :Ω, :V, :W, :J, :Pa, :N, :C, :L, :€] ]
FlexUnits.preferred_units(::Type{<:MoneyDimensions}) = CURRENCY_UNITS


#===============================================================================================================================
If you want to combine operations with default unit registry
    => Add promotion rules and a conversion function 
===============================================================================================================================#
Base.promote_rule(::Type{MoneyDimensions{T1}}, ::Type{Dimensions{T2}}) where {T1, T2} = MoneyDimensions{promote_type(T1,T2)}
Base.promote_rule(::Type{Dimensions{T1}}, ::Type{MoneyDimensions{T2}}) where {T1, T2} = MoneyDimensions{promote_type(T1,T2)}


@testset "Custom Dimensions" begin

@test string(1*CurrencyUnits.u"Pa") == "1.0 kg/(m s²)"
@test string(1*CurrencyUnits.u"Pa" |> simplify) == "1.0 Pa"

@test (6*UnitRegistry.u"kg") + (5*CurrencyUnits.u"kg") ≈ 11*CurrencyUnits.u"kg"
@test typeof((6*UnitRegistry.u"kg") + (5*CurrencyUnits.u"kg")) == Quantity{Float64, StaticDims{MoneyDimensions(kg=1)}}

@test (5*CurrencyUnits.u"EUR") / (1*UnitRegistry.u"hr") ≈ 5*CurrencyUnits.u"EUR/hr"
@test typeof((5*CurrencyUnits.u"EUR") / (1*UnitRegistry.u"hr")) == Quantity{Float64, StaticDims{MoneyDimensions(€=1, s=-1)}}


end
