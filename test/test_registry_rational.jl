using Test
using FlexUnits
using Random

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

#x*s correctly rounded: the exact product, rounded once
exact(x::T, s::Rational) where T<:AbstractFloat = setprecision(BigFloat, 4096) do
    T(big(x)*numerator(s)/denominator(s))
end

@testset "A rational scale applied to a float is rounded once" begin
    m  = RationalRegistry.u"m"
    mm = RationalRegistry.u"mm"
    inch = RationalRegistry.u"inch"
    @test ustrip(9.0*mm) === 0.009
    @test ustrip(142.0*mm) === 0.142            #142.0*0.001 was 0.14200000000000002
    @test ustrip(mm, 0.142*m) === 142.0         #0.142/0.001 was 141.99999999999997
    @test ustrip(inch, 0.0508*m) === 2.0

    rng = Xoshiro(42)
    samples = vcat(
        [parse(Float64, string(k, "e", e)) for k in 1:999 for e in (-3, 0, 3)], #short decimals
        [ldexp(rand(rng), rand(rng, -60:60)) for _ in 1:20_000],                #random magnitudes
        [-1.5, 1.0, 3.0, 7.0, 1e300, 1e-300, 5e-324, floatmax(Float64)],
    )
    units = (RationalRegistry.u"km", RationalRegistry.u"cm", mm, RationalRegistry.u"μm",
             RationalRegistry.u"nm", inch, RationalRegistry.u"ft")
    for u in units
        s = uscale(u)
        @test all(x -> ustrip(x*u) === exact(x, s), samples)                 #to base units
        @test all(x -> ustrip(u, x*m) === exact(x, inv(s)), samples)         #from base units
        @test all(x -> uconvert(mm, u)(x) === exact(x, s*1000), samples)     #between units: one composed scale
    end

    #Signs, zeros and non-finite values
    @test ustrip(-142.0*mm) === -0.142
    @test ustrip(-0.0*mm) === -0.0 && ustrip(0.0*mm) === 0.0
    @test isnan(ustrip(NaN*mm)) && ustrip(Inf*mm) === Inf && ustrip(-Inf*inch) === -Inf

    #Float32 rounds once too
    @test ustrip(9f0*mm) === 0.009f0
    @test all(x -> ustrip(x*inch) === exact(x, 127//5000), Float32.(samples[1:3000]))

    #An integer or a rational stays exact, an offset (°C) is added after
    @test ustrip(142*mm) === 71//500
    tc = AffineTransform{Rational{Int64}}(scale=1//1, offset=5463//20)
    @test tc(25.0) === 25.0 + 273.15
end
