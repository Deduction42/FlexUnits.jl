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

        y = x - x
        @test iszero(x) == false
        @test iszero(y) == true
        @test iszero(dimension(y)) == false
    
        y = -x
        @test ustrip(y) == -ustrip(x)
        @test dimension(y) == dimension(x)
    
        y = x / x
        @test iszero(dimension(x)) == false
        @test iszero(dimension(y)) == true
    
        y = x * Inf32
        @test typeof(y).parameters[1] == promote_type(T, Float32)
        @test typeof(y).parameters[2] == D
        @test isfinite(x)
        @test !isfinite(y)
       
        u = Dimensions(length=2//5)
        x = RealQuantity(-1.2, u)

        @test typemax(x) == RealQuantity(typemax(-1.2), u)
    
        @test abs(x) == RealQuantity(1.2, u)
        @test abs(x) == abs(RealQuantity(1.2, u))
        @test abs2(x) == RealQuantity(abs2(-1.2), u^2)
    
        @test copy(x) == x
    
        @test iszero(x) == false
        @test iszero(x * 0) == true
        @test isfinite(x) == true
        @test isfinite(x * Inf) == false
        @test isfinite(x * NaN) == false
        @test isinf(x * Inf) == true
        @test isnan(x) == false
        @test isnan(x * NaN) == true
        @test isreal(x) == true
        @test isreal(x * (1 + 2im)) == false
        @test signbit(x) == true
        @test signbit(-x) == false
        @test isempty(x) == false
        @test isempty(Quantity([0.0, 1.0], u)) == false
        @test isempty(Quantity(Float64[], u)) == true 

    end


    #Test parsing of non-pretty unit output
    x = quantity(0.2, Dimensions(length=1, mass=2.5, time=-1))
    show(tmp_io, unit(x), pretty=false)
    xp = ustrip(x)*uparse(String(take!(tmp_io)))
    @test ubase(xp) ≈ ubase(x)
end

