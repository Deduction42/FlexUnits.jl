using Revise
using Test
using DimensionalUnits
using DimensionalUnits: FAST_RATIONAL
using .UnitRegistry

@testset "basic utilities" begin
    tmp_io = IOBuffer()

    for Q in [Quantity, NumberQuantity, RealQuantity], T in [Float16, Float32, Float64], R in [FAST_RATIONAL, Rational{Int16}, Rational{Int32}]
        
        D = Dimensions{R}
        x = Q(T(0.2), D(length=1, mass=2.5, time=-1))

        @test typeof(x).parameters[1] == T
        @test typeof(x).parameters[2] == D
        @test ustrip(x) ≈ T(0.2)
        @test dimension(x) == D(length=1, mass=5//2, time=-1)
        if R == FAST_RATIONAL
            @test dimension(x) == Dimensions(length=1, mass=5//2, time=-1)
        end

        y = x^2

        @test typeof(y) <: RealQuantity
        @test typeof(x).parameters[1] == T
        @test typeof(x).parameters[2] == D
        @test ustrip(y) ≈ T(0.04)

        if R <: Rational
            DimensionalUnits.PRETTY_DIM_OUTPUT[] = true    
            @test string(x) == "0.2 (m kg⁵ᐟ²)/s"
            @test string(inv(x)) == "5.0 s/(m kg⁵ᐟ²)"
        end
    end

    #Test parsing of non-pretty unit output
    x = quantity(0.2, Dimensions(length=1, mass=2.5, time=-1))
    show(tmp_io, unit(x), pretty=false)
    xp = ustrip(x)*uparse(String(take!(tmp_io)))
    @test ubase(xp) ≈ ubase(x)

end

