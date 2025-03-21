using Revise
using Test
using DimensionalUnits
using DimensionalUnits: FAST_RATIONAL, FixedRational
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


@testset "basic unit parsing" begin
    x = 1.3u"km/s^2"
    @test ustrip(x) == 1.3
    @test ustrip_base(x) == 1300  # SI base units

    y = 0.9u"sqrt(mΩ)"
    @test typeof(y) == RealQuantity{Float64, AffineUnits{Dimensions{FAST_RATIONAL}}}
    @test typeof(ubase(y)) == RealQuantity{Float64, Dimensions{FAST_RATIONAL}}
    @test ustrip_base(y) ≈ 0.02846049894151541


    y = BigFloat(0.3) * u"mΩ"
    @test typeof(y) == RealQuantity{BigFloat, AffineUnits{Dimensions{FAST_RATIONAL}}}
    @test ustrip_base(y) ≈ 0.0003

    y32 = convert(RealQuantity{Float32, AffineUnits{Dimensions{FAST_RATIONAL}}}, y)
    @test typeof(y32) == RealQuantity{Float32, AffineUnits{Dimensions{FAST_RATIONAL}}}

    z = 1*u"yr"
    @test ustrip_base(z) ≈ 60 * 60 * 24 * 365.25
    @test z == 1*uparse("yr")

    # Test type stability of extreme range of units
    @test typeof(u"s"^2) == AffineUnits{Dimensions{FAST_RATIONAL}}
    @test typeof(u"Ω") == AffineUnits{Dimensions{FAST_RATIONAL}}

    @test_throws ArgumentError uparse(":x")
    @test_throws "Symbol x not found" uparse("x")
    @test_throws "Unexpected expression" uparse("import ..Units")
    @test_throws "Unexpected expression" uparse("(m, m)")
    @test_throws LoadError eval(:(us"x"))
end

@testset "Additional tests of FixedRational" begin
    @test convert(Int64, FixedRational{Int64,1000}(2 // 1)) == 2
    @test convert(Int32, FixedRational{Int64,1000}(3 // 1)) == 3
    @test convert(Bool, FixedRational{Int8,6}(1//1)) === true
    @test convert(Bool, FixedRational{Int8,6}(0//1)) === false

    @test_throws InexactError convert(Int32, FixedRational{Int8,6}(2//3))
    @test_throws InexactError convert(Bool, FixedRational{Int8,6}(2//1))

    @test_throws "Refusing to" promote(FixedRational{Int,10}(2), FixedRational{Int,4}(2))

    f64 = FixedRational{Int,10}(2)
    f8 = FixedRational{Int8,10}(2)
    @test promote(f64, f8) == (2, 2)
    @test typeof(promote(f64, f8)) == typeof((f64, f64))
    @test typeof(promote(FixedRational{Int8,10}(2), FixedRational{Int8,10}(2))) == typeof((f8, f8))
    @test promote_type(Float64, typeof(f64)) == Float64

    # Required to hit integer branch (otherwise will go to `literal_pow`)
    f(i::Int) = Dimensions(length=1, mass=-1)^i
    @test f(2) == Dimensions(length=2, mass=-2)

    # Null conversion
    @test typeof(FixedRational{Int,10}(FixedRational{Int,10}(2))) == FixedRational{Int,10}

    # Conversion to Rational without specifying type
    @test convert(Rational, FixedRational{UInt8,6}(2)) === Rational{UInt8}(2)

    # Promotion rules
    @test promote_type(FixedRational{Int64,10},FixedRational{BigInt,10}) == FixedRational{BigInt,10}
    @test promote_type(Rational{Int8}, FixedRational{Int,12345}) == Rational{Int}
    @test promote_type(Int8, FixedRational{Int,12345}) == FixedRational{Int,12345}

    # Bug where user would create a FixedRational{::Type{Int32}, ::Int64} and get stack overflow,
    # because the stored type was FixedRational{::Type{Int32}, ::Int32}
    x = 10u"m"
    user_quantity = Quantity(10.0, Dimensions{FixedRational{Int32,25200}}(1, 0, 0, 0, 0, 0, 0))
    @test x == user_quantity
end