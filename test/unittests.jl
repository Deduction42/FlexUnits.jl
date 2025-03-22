using Revise
using Test
using DimensionalUnits
using DimensionalUnits: FAST_RATIONAL, FixedRational
using .UnitRegistry

@testset "Basic utilities" begin

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
    U = typeof(first(values(UnitRegistry.UNITS)))
    @test typeof(u"s"^2) == U
    @test typeof(u"Ω") == U

    @test_throws ArgumentError uparse(":x")
    @test_throws "Symbol x not found" uparse("x")
    @test_throws "Unexpected expression" uparse("import ..Units")
    @test_throws "Unexpected expression" uparse("(m, m)")
    @test_throws LoadError eval(:(us"x"))
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

    # Test that regular type promotion applies:
    q = Quantity(2, d)
    @test typeof(q) == Quantity{Int64,typeof(d)}
    @test typeof(q ^ 2) == RealQuantity{Int64,typeof(d)}
    @test typeof(0.5 * q) == RealQuantity{Float64,typeof(d)}
    @test typeof(inv(q)) == RealQuantity{Float64,typeof(d)}

    # Automatic conversions via constructor:
    for T in [Float16, Float32, Float64, BigFloat], R in [FAST_RATIONAL, Rational{Int16}, Rational{Int32}]
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

#=
@testset "AffineUnits" begin
    °C  = ua"°C"
    °F  = ua"°F"
    mps = ua"m/s"
    
    import DynamicQuantities.AffineUnits
    import DynamicQuantities.affine_base_dim

    @test °C == AffineUnits.°C
    @test °C == AffineUnits.°C

    @test aff_uparse("m/(s^2.5)") == ua"m/(s^2.5)"
    @test_throws ArgumentError aff_uparse("s[1]")
    @test_throws ArgumentError aff_uparse("pounds_per_hour")
    @test °C isa Quantity{T,AffineDimensions{R}} where {T,R}
    @test dimension(°C) isa AffineDimensions
    @test dimension(°C) isa AbstractAffineDimensions

    @test affine_base_dim(dimension(°C)).temperature == 1
    @test affine_base_dim(dimension(°C)).length == 0

    @test inv(mps) == us"s/m"
    @test inv(mps) == u"s/m"
    @test mps^2 == u"m^2/s^2"

    @test °C == ua"degC"
    @test °F == ua"degF"
    @test dimension(°C) == dimension(ua"degC")
    @test_throws AffineOffsetError 5.0ua"°C" - 4.0ua"°C"

    # Constructors
    @test with_type_parameters(AffineDimensions, Float64) == AffineDimensions{Float64}
    @test constructorof(AffineDimensions) == AffineDimensions
    @test constructorof(AffineDimensions{Float64}) == AffineDimensions
    @test Quantity(1.0, AffineDimensions(dimension(u"K"))) == u"K"
    @test AffineDimensions(scale=1, offset=0, basedim=dimension(u"K")) == AffineDimensions(basedim=dimension(u"K"))
    @test AffineDimensions(scale=1, offset=0, basedim=u"K") == AffineDimensions(basedim=dimension(ua"K"))
    @test AffineDimensions(scale=1.0, offset=273.15u"K", basedim=dimension(u"K")) == AffineDimensions(basedim=ua"°C")

    kelvin  = AffineDimensions(basedim=u"K")
    @test Quantity(1.0, kelvin) == u"K"

    rankine = AffineDimensions(scale=5/9, offset=0.0, basedim=dimension(u"K"))
    @test Quantity(1.0, rankine) == (5/9)u"K"

    fahrenheit = AffineDimensions(scale=1.0, offset=Quantity(459.67, rankine), basedim=rankine)
    @test Quantity(1.0, fahrenheit) ≈ °F

    celsius = AffineDimensions(scale=9/5, offset=Quantity(32.0, rankine), basedim=°F)
    @test Quantity(1.0, celsius) ≈ °C

    # Round-trip sanity checks
    @test -40°C ≈ -40°F
    @test Quantity(-40.0, celsius) ≈ Quantity(-40.0, fahrenheit)
    
    # Test promotion explicitly for coverage:
    @test promote_type(AffineDimensions{Int16}, AffineDimensions{Int32}) === AffineDimensions{Int32}
    @test promote_type(Dimensions{Int16}, AffineDimensions{Int32}) === Dimensions{Int32}
    @test promote_type(AffineDimensions{Int16}, Dimensions{Int32}) === Dimensions{Int32}
    @test promote_type(SymbolicDimensions{Int16}, AffineDimensions{Int32}) === Dimensions{Int32}
    @test promote_type(AffineDimensions{Int16}, SymbolicDimensions{Int32}) === Dimensions{Int32}

    # Type conversions
    @test convert(Quantity{Float64, AffineDimensions}, u"kg") isa Quantity{Float64, AffineDimensions{DynamicQuantities.DEFAULT_DIM_BASE_TYPE}}
    @test convert(Quantity{Float64, AffineDimensions{Float64}}, u"kg") isa Quantity{Float64, AffineDimensions{Float64}}
    @test convert(Quantity{Float64, Dimensions}, ua"kg") isa Quantity{Float64, Dimensions{DynamicQuantities.DEFAULT_DIM_BASE_TYPE}}
    
    # Test uncovered operations
    @test (2.0ua"m")^2 == (2.0u"m")^2
    @test dimension(ua"m")^Int32(2) == dimension(ua"m^2")
    @test 2.0u"m" + 2.0ua"m" === 4.0u"m"
    @test 2.0ua"m" + 2.0ua"m" === 4.0u"m"
    @test 2.0u"m" - 2.0ua"m" === 0.0u"m"
    @test 2.0ua"m" - 2.0ua"cm" === 1.98u"m"
    @test mod(2.0u"m", 2.0ua"m") === 0.0u"m"
    @test mod(2.0ua"m", 2.0ua"m") === 0.0u"m"
    @test_throws AffineOffsetError mod(2.0ua"°C", 2.0ua"°C")

    @test 2.0u"K" ≈ 2.0ua"K"
    @test 2.0ua"K" ≈ 2.0ua"K"
    @test 2.0ua"K" ≈ 2.0u"K"
    @test_throws AffineOffsetError (2ua"°C")^2
    @test uexpand(2ua"°C") == 275.15u"K"

    # Test conversions
    @test °C |> us"K" isa Quantity{<:Real, <:SymbolicDimensions}
    @test 0°C |> us"K" == 273.15us"K"
    @test us"K" |> °C isa Quantity{<:Real, <:AffineDimensions}
    @test 0us"K" |> °C  == -273.15°C
    @test °C |> °F isa Quantity{<:Real, <:AffineDimensions}
    @test 0°C |> °F == 32°F

    @test QuantityArray([0,1]°C) |> °F isa QuantityArray{T, <:Any, AffineDimensions{R}} where {T,R}

    # Test display against errors
    celsius = AffineDimensions(offset=273.15, basedim=u"K")
    psi = AffineDimensions(basedim=6.89476us"kPa")
    io = IOBuffer()
    @test isnothing(show(io, (dimension(°F), dimension(ua"K"), psi, celsius, fahrenheit)))

    # Test updating affine units
    @test isnothing(DynamicQuantities.update_external_affine_unit(:°C, °C)) # same value yields nothing for quantity
    @test isnothing(DynamicQuantities.update_external_affine_unit(:°C, dimension(°C))) # same value yields nothing for dimension with symbol
    @test isnothing(DynamicQuantities.update_external_affine_unit(dimension(°C))) # same value yields for dimension
    @test_throws "Cannot register an affine dimension without symbol declared" DynamicQuantities.update_external_affine_unit(celsius) # cannot register :nothing

    # Cannot re-register a unit if its value changes
    @test_throws "Unit `°C` already exists as `(scale = 1.0, offset = 273.15, basedim = K)`, its value cannot be changed to `(scale = 2.0, offset = 273.15, basedim = K)`" DynamicQuantities.update_external_affine_unit(:°C, 2*°C)

    # Test map_dimensions
    @test map_dimensions(+, dimension(ua"m/s"), dimension(ua"m/s")) == AffineDimensions(Dimensions(length=2, time=-2))
    @test map_dimensions(-, dimension(ua"m"), dimension(ua"s"))     == AffineDimensions(Dimensions(length=1, time=-1))
    @test map_dimensions(Base.Fix1(*,2), dimension(ua"m/s"))        == AffineDimensions(Dimensions(length=2, time=-2))

end

=#





















@testset "Additional tests of FixedRational" begin
    #Tests were basically copied from DynamicQuantities, since FixedRational isn't its own package
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