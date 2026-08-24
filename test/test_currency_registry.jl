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

using FlexUnits.RegistryTools
import ..MoneyDimensions

const UNITS = PermanentDict{Symbol,Units{MoneyDimensions{FixRat32},AffineTransform{Float64}}}()
registry_defaults!(UNITS)
register_unit!(UNITS, "EUR" => UNITS[:€])

const PREFERRED_UNITS = [UNITS[u] for u in [:F, :H, :T, :Ω, :V, :W, :J, :Pa, :N, :C, :L]]

@generate_unit_simplifier(PREFERRED_UNITS)
@generate_registry_exports(UNITS)
end


@testset "Custom Dimensions" begin

@test string(1*CurrencyUnits.u"Pa") == "1.0 kg/(m s²)"
@test string(1*CurrencyUnits.u"Pa" |> simplify) == "1.0 Pa"

@test (6*UnitRegistry.u"kg") + (5*CurrencyUnits.u"kg") ≈ 11*CurrencyUnits.u"kg"
@test typeof((6*UnitRegistry.u"kg") + (5*CurrencyUnits.u"kg")) == Quantity{Float64, StaticDims{MoneyDimensions(kg=1)}}

@test (5*CurrencyUnits.u"EUR") / (1*UnitRegistry.u"hr") ≈ 5*CurrencyUnits.u"EUR/hr"
@test typeof((5*CurrencyUnits.u"EUR") / (1*UnitRegistry.u"hr")) == Quantity{Float64, StaticDims{MoneyDimensions(€=1, s=-1)}}


end
