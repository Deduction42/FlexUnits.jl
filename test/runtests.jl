using FlexUnits
using TestItems: @testitem
using TestItemRunner

#============================================================================================================================
Run these commands at startup to see coverage
julia --startup-file=no --depwarn=yes --threads=auto -e 'using Coverage; clean_folder(\"src\"); clean_folder(\"test\")'
julia --startup-file=no --depwarn=yes --threads=auto --code-coverage=user --project=. -e 'using Pkg; Pkg.test(coverage=true)'
julia --startup-file=no --depwarn=yes --threads=auto coverage.jl
============================================================================================================================#

#To see the actual coverage in VSCode, install the Coverage Gutters extension
#https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters

using Revise
using Test
using BenchmarkTools
using FlexUnits, .UnitRegistry
using FlexUnits: DEFAULT_RATIONAL, FixedRational, map_dimensions
using Aqua

const DEFAULT_UNIT_TYPE = typeof(first(values(UnitRegistry.UNITS)))
const DEFAULT_DIM_TYPE  = FlexUnits.dimtype(DEFAULT_UNIT_TYPE)

@testset "Basic utilities" begin

    for Q in [Quantity, NumberQuantity, RealQuantity], T in [Float16, Float32, Float64], R in [DEFAULT_RATIONAL, Rational{Int16}, Rational{Int32}]
        
        D = Dimensions{R}
        x = Q(T(0.2), D(length=1, mass=2.5, time=-1))

        @test typeof(x).parameters[1] == T
        @test typeof(x).parameters[2] == D
        @test ustrip(x) ≈ T(0.2)
        @test dimension(x) == D(length=1, mass=5//2, time=-1)
        if R == DEFAULT_RATIONAL
            @test dimension(x) == Dimensions(length=1, mass=5//2, time=-1)
        end

        y = x^2

        @test typeof(y) <: RealQuantity
        @test typeof(x).parameters[1] == T
        @test typeof(x).parameters[2] == D
        @test ustrip(y) ≈ T(0.04)

        if R <: Rational
            FlexUnits.PRETTY_DIM_OUTPUT[] = true    
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
        @test sin(5u"rad") == sin(5) 
        @test cos(60u"deg") ≈ 0.5

        #Cannot check iszero on non-affine units 
        @test_throws NotScalarError iszero(5u"°C")

    end

end


@testset "Basic unit functionality" begin
    x = quantity(0.2, Dimensions(length=1, mass=2.5, time=-1))
    u = ustrip(x)

    #Test round-trip parsing
    @test ubase(x) ≈ ustrip(x)*uparse(string(unit(x), pretty=false))

    #Test toggling pretty-print
    pretty_print_units(false)
    xp = ustrip(x)*uparse(string(unit(x)))
    pretty_print_units(true)
    @test ubase(xp) ≈ ubase(x)


    x = 1.3u"km/s^2"
    @test ustrip(x) == 1.3
    @test ustrip_base(x) == 1300  # SI base units
    @test x == q"1.3km/s^2"
    @test x !== q"1.3km/s^2"
    @test abs(x) === q"1.3km/s^2"

    y = 0.9u"sqrt(mΩ)"
    @test typeof(y) == RealQuantity{Float64, AffineUnits{Dimensions{DEFAULT_RATIONAL}}}
    @test typeof(ubase(y)) == RealQuantity{Float64, Dimensions{DEFAULT_RATIONAL}}
    @test ustrip_base(y) ≈ 0.02846049894151541
    @test y ≈ q"(0.9*sqrt(mΩ))"

    y = BigFloat(0.3) * u"mΩ"
    @test typeof(y) == RealQuantity{BigFloat, AffineUnits{Dimensions{DEFAULT_RATIONAL}}}
    @test ustrip_base(y) ≈ 0.0003

    y32 = convert(RealQuantity{Float32, AffineUnits{Dimensions{DEFAULT_RATIONAL}}}, y)
    @test typeof(y32) == RealQuantity{Float32, AffineUnits{Dimensions{DEFAULT_RATIONAL}}}

    z = 1*u"yr"
    @test ustrip_base(z) ≈ 60 * 60 * 24 * 365.25
    @test z == 1*uparse("yr")
    @test z == qparse("1yr")
    @test 1/z == qparse("1/yr")
    @test_throws MethodError qparse("yr")

    # Test type stability of extreme range of units
    U = typeof(first(values(UnitRegistry.UNITS)))
    @test typeof(u"s"^2) == U
    @test typeof(u"Ω") == U

    @test_throws ArgumentError uparse(":x")
    @test_throws "Symbol x not found" uparse("x")
    @test_throws "Unexpected expression" uparse("import ..Units")
    @test_throws "Unexpected expression" uparse("(m, m)")
    @test_throws LoadError eval(:(us"x"))


    #Basic mathematical operations
    xp = 1u"percent"
    xv = 1u"m/s"

    @test xp + 1 == 1.01
    @test 1 + xp == 1.01
    @test xv + xv == 2u"m/s"
    @test (xv + xv) isa RealQuantity{Float64, <:Dimensions}

    @test 1 - xp == 0.99 
    @test xp - 1 == -0.99
    @test xv - xv == 0u"m/s"
    @test (xv - xv) isa RealQuantity{Float64, <:Dimensions}

    @test_throws DimensionError xv + 1
    @test_throws DimensionError 1 + xv

    @test xv*1 == 1u"m/s"
    @test xv*1 !== 1u"m/s"
    @test xv*1 === ubase(1u"m/s")
    @test xp*1 == 0.01*u""
    @test xv*xv == 1u"m^2/s^2"

    @test xv/1 == 1u"m/s"
    @test xp/1 == 0.01*u""
    @test xv/xv == 1u""

    @test sqrt(4*xv) == 2u"m^0.5/s^0.5"
    @test (4*xv)^0.5 == 2u"m^0.5/s^0.5"
    @test (2*xv)^2 == 4u"m^2/s^2"

    @test_throws DimensionError 2.0^(1u"m/s")
    @test_throws DimensionError (1u"m/s")^(1u"m/s")

    @test_throws ArgumentError uparse("s[1]")
    @test_throws ArgumentError uparse("pounds_per_hour")


    #Tests on truly affine units
    °C  = u"°C"
    °F  = u"°F"
    K   = u"K"
    @test dimension(°C).temperature == 1
    @test dimension(°C).length == 0
    @test °C == u"degC"
    @test °F == u"degF"
    @test 5°C - 4°C == 1K
    @test 5°C + 4°C == 282.15°C
    @test 5°C + 4°C - 0°C == 9°C
    @test qparse("0°C") == 0°C
    @test 0*uparse("°C") == 0°C


    # Constructors
    kelvin  = AffineUnits(dims=u"K")
    @test 1*kelvin == 1K

    rankine = AffineUnits(scale=5/9, offset=0.0, dims=K)
    @test 1*rankine == (5/9)K

    fahrenheit = AffineUnits(scale=5/9, offset=(273.15-32*5/9), dims=K)
    @test 1*fahrenheit ≈ 1°F

    celsius = AffineUnits(offset=273.15, dims=K)
    @test 1*celsius ≈ 1°C

    # Round-trip sanity checks
    @test -40°C ≈ -40°F
    @test -40.0*celsius ≈ -40.0*fahrenheit
    
    # Test promotion explicitly for coverage:
    @test promote_type(AffineUnits{Dimensions{Int16}}, AffineUnits{Dimensions{Int32}}) === AffineUnits{Dimensions{Int32}}
    @test promote_type(Dimensions{Int16}, AffineUnits{Dimensions{Int32}}) === AffineUnits{Dimensions{Int32}}
    
    # Test conversions
    @test 1°C |> K isa RealQuantity{<:Real, <:AffineUnits}
    @test 1°C |> unit(ubase(1K)) isa RealQuantity{<:Real, <:Dimensions}
    @test  °C |> unit(ubase(1K)) isa AffineTransform

    @test 0°C |> u"K" == 273.15u"K"
    @test 1u"K" |> °C isa RealQuantity{<:Real, <:AffineUnits}
    @test 0u"K" |> °C  == -273.15°C
    @test °C |> °F isa AffineTransform
    @test 0°C |> °F == 32°F


    # Test display against errors
    celsius = AffineUnits(offset=273.15, dims=u"K")
    psi = FlexUnits.asunit(6.89476u"kPa")
    io = IOBuffer()
    @test isnothing(show(io, (dimension(°F), dimension(u"K"), psi, celsius, fahrenheit)))

    # Test updating affine units
    @test register_unit("°C" => °C) isa AbstractDict # Updating the same value does nothing for unit
    @test register_unit("K" => Dimensions(temperature=1)) isa AbstractDict # same value yields nothing for dimension
    @test_throws MethodError register_unit(u"K") # cannot register only a unit

    # Cannot re-register a unit if its value changes
    @test_throws RegistryTools.PermanentDictError register_unit("°C"=>u"°F")

    # Cannot register non-parsable unit
    @test_throws ArgumentError register_unit("m/s"=>u"m/s")

    # Test map_dimensions
    @test map_dimensions(+, dimension(u"m/s"), dimension(u"m/s")) == Dimensions(length=2, time=-2)
    @test map_dimensions(-, dimension(u"m"), dimension(u"s"))     == Dimensions(length=1, time=-1)
    @test map_dimensions(Base.Fix1(*,2), dimension(u"m/s"))       == Dimensions(length=2, time=-2)
end


@testset "Type conversions" begin
    d = Dimensions{Rational{Int16}}(mass=2)
    d32 = convert(Dimensions{Rational{Int32}}, d)
    @test typeof(d) == Dimensions{Rational{Int16}}
    @test typeof(d32) == Dimensions{Rational{Int32}}

    # Should not change:
    @test convert(Dimensions{Rational{Int16}}, d) === d

    q = Quantity(0.5, d)
    q32_32 = convert(Quantity{Float32,Dimensions{Rational{Int32}}}, q)
    @test typeof(q) == Quantity{Float64,Dimensions{Rational{Int16}}}
    @test typeof(q32_32) == Quantity{Float32,Dimensions{Rational{Int32}}}
    @test ustrip(q) == 0.5
    @test ustrip(q32_32) == 0.5
    @test typeof(ustrip(q)) == Float64
    @test typeof(ustrip(q32_32)) == Float32
    @test dimension(q32_32) == dimension(q)
    @test typeof(convert(Quantity{Float16}, q)) == Quantity{Float16,Dimensions{Rational{Int16}}}
    @test convert(Quantity, q) === q
    @test convert(Quantity{Float64, AffineUnits{DEFAULT_DIM_TYPE}}, ubase(1u"kg")) isa Quantity{Float64, AffineUnits{DEFAULT_DIM_TYPE}}

    # Test that regular type promotion applies:
    q = Quantity(2, d)
    @test typeof(q) == Quantity{Int64,typeof(d)}
    @test typeof(q ^ 2) == RealQuantity{Int64,typeof(d)}
    @test typeof(0.5 * q) == RealQuantity{Float64,typeof(d)}
    @test typeof(inv(q)) == RealQuantity{Float64,typeof(d)}

    # Automatic conversions via constructor:
    for T in [Float16, Float32, Float64, BigFloat], R in [DEFAULT_RATIONAL, Rational{Int16}, Rational{Int32}]
        D = Dimensions{R}
        q = Quantity{T,D}(2, D(length=1.5))
        @test typeof(q) == Quantity{T,D}
        @test typeof(ustrip(q)) == T

        # Now, without R, the default will be DEFAULT_DIM_BASE_TYPE:
        q = Quantity{T}(2, D(length=1.5))
        @test typeof(q) == Quantity{T,D}
        @test typeof(ustrip(q)) == T

        # Just dimensions:
        d = D(length=1.5)
        @test typeof(d) == D
    end
end


@testset "Unit conversions (uconvert)" begin
    U = typeof(first(values(UnitRegistry.UNITS)))
    @test uconvert(u"nm", 5e-9u"m") ≈ (5e-9u"m" |> u"nm") ≈ 5u"nm"
    @test_throws ConversionError uconvert(u"nm*J", 5e-9u"m")

    # Types:
    @test typeof(uconvert(u"nm", 5e-9u"m")) <: RealQuantity{Float64, U} 
    @test typeof(uconvert(u"nm", Quantity(5e-9, u"m"))) <: RealQuantity{Float64, U}
    @test uconvert(u"nm", Quantity(5e-9, u"m")) ≈ 5u"nm"

    @test dimension(1u"m" |> u"nm")[:length] == 1

    for Q in (RealQuantity, NumberQuantity, Quantity)
        # Different types require converting both arguments:
        q = convert(Q{Float16}, 1.5u"g")

        # Broadcasting conversions over Arrays
        x = [1.0, 2.0, 3.0] .* Q(1, u"kg")
        x2 = x .|> u"g"
        @test typeof(x2) <: Vector{<:RealQuantity{Float64,<:AffineUnits{<:Any}}}
        @test x2[2] ≈ Q(2000u"g")
    end
end

@testset "Additional tests of FixedRational" begin
    #Tests were basically copied from DynamicQuantities, since FixedRational isn't its own package
    @test convert(Int64, FixedRational{1000,Int64}(2 // 1)) == 2
    @test convert(Int32, FixedRational{1000,Int64}(3 // 1)) == 3
    @test convert(Bool, FixedRational{6,Int8}(1//1)) === true
    @test convert(Bool, FixedRational{6,Int8}(0//1)) === false

    @test_throws InexactError convert(Int32, FixedRational{6,Int8}(2//3))
    @test_throws InexactError convert(Bool, FixedRational{6,Int8}(2//1))

    @test_throws "Refusing to" promote(FixedRational{10,Int}(2), FixedRational{4,Int}(2))

    f64 = FixedRational{10,Int}(2)
    f8  = FixedRational{10,Int8}(2)
    @test promote(f64, f8) == (2, 2)
    @test typeof(promote(f64, f8)) == typeof((f64, f64))
    @test typeof(promote(FixedRational{10,Int8}(2), FixedRational{10,Int8}(2))) == typeof((f8, f8))
    @test promote_type(Float64, typeof(f64)) == Float64

    # Required to hit integer branch (otherwise will go to `literal_pow`)
    f(i::Int) = Dimensions(length=1, mass=-1)^i
    @test f(2) == Dimensions(length=2, mass=-2)

    # Null conversion
    @test typeof(FixedRational{10,Int}(FixedRational{10,Int}(2))) == FixedRational{10,Int}

    # Conversion to Rational without specifying type
    @test convert(Rational, FixedRational{6,UInt8}(2)) === Rational{UInt8}(2)

    # Promotion rules
    @test promote_type(FixedRational{10,Int64},FixedRational{10,BigInt}) == FixedRational{10,BigInt}
    @test promote_type(Rational{Int8}, FixedRational{12345,Int}) == Rational{Int}
    @test promote_type(Int8, FixedRational{12345,Int}) == FixedRational{12345,Int}

    # Bug where user would create a FixedRational{::Type{Int32}, ::Int64} and get stack overflow,
    # because the stored type was FixedRational{::Type{Int32}, ::Int32}
    x = 10u"m"
    user_quantity = Quantity(10.0, Dimensions{FixedRational{25200,Int32}}(1, 0, 0, 0, 0, 0, 0))
    @test x == user_quantity
end

@testset "Aqua.jl" begin
    Aqua.test_all(FlexUnits)
end

R = Ref(8.314 * u"J/(mol*K)")
@benchmark let R=R[]
    v_satp = R * (25u"°C") / (101.3u"kPa")
end


