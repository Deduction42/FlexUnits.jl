using Test
using FlexUnits

module RationalRegistry
    using FlexUnits.RegistryTools #RegistryTools contains all you need to build a registry in one simple import

    const UNITS = PermanentDict{Symbol, Units{Dimensions{FixRat32}, AffineTransform{Rational{Int64}}}}() #Just change the AffineTransform type
    registry_defaults!(UNITS) #Auto-populate the new registry

    @generate_registry_exports(UNITS) #Use macros to generate the boilerplate code for registry exports
end

@testset "Rational Registry Exact Conversions" begin
    tc = 1*RationalRegistry.u"°C" 
    tf = tc |> RationalRegistry.u"°F"

    @test ustrip(tc) isa Rational 
    @test ustrip(tf) isa Rational
end