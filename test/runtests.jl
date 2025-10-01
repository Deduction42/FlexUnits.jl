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

#using Revise
using Test
using BenchmarkTools
using FlexUnits, .UnitRegistry
using FlexUnits: DEFAULT_RATIONAL, FixedRational, map_dimensions
using Aqua

const DEFAULT_UNIT_TYPE = typeof(first(values(UnitRegistry.UNITS)))
const DEFAULT_DIM_TYPE  = FlexUnits.dimtype(DEFAULT_UNIT_TYPE)

@testset "Basic utilities" begin
    d = Dimensions(length=2, mass=1, time=-2)
    @test d.length == 2
    @test d.mass == 1
    @test d.time == -2
    @test FlexUnits.dimtype(typeof(d)) == DEFAULT_DIM_TYPE
    @test FlexUnits.dimtype(typeof(u"m/s")) == DEFAULT_DIM_TYPE
    @test uscale(d) == 1
    @test uoffset(d) == 0
    @test FlexUnits.usymbol(d) == FlexUnits.DEFAULT_USYMBOL

    d = NoDims()
    @test d.length == 0
    @test d.mass == 0
    @test FlexUnits.unit_symbols(NoDims) ==  NoDims{Symbol}()

    @test FlexUnits.remove_offset(u"°C") == u"K"

    @test FlexUnits.constructorof(typeof(Dimensions())) == Dimensions
    @test FlexUnits.constructorof(typeof(u"m")) == AffineUnits
    @test FlexUnits.constructorof(typeof(1.0*u"m")) == Quantity
    @test FlexUnits.constructorof(typeof((1.0+im)*u"m")) == Quantity
    @test FlexUnits.constructorof(typeof(quantity("this", u"m"))) == Quantity
    @test FlexUnits.constructorof(Array{Float64}) == Array

    @test string(AffineUnits(scale=1, offset=0, dims=dimension(u"m"), symbol=:_)) == "AffineUnits(scale=1.0, offset=0.0, dims=m)"
    @test string(ubase(1.0u"kg*m^2/s^2"), pretty=true)  == "1.0 (m² kg)/s²"
    @test string(ubase(1.0u"kg*m^2/s^2"), pretty=false) == "(1.0)(m^2*kg)/s^2"
    @test string(1.0u"kg*m^2/s^2", pretty=true)  == "1.0 kg*m^2/s^2"
    @test string(1.0u"kg*m^2/s^2", pretty=false) == "(1.0)kg*m^2/s^2"
    @test string(ubase(1.0u"1/s^2"), pretty=true)  == "1.0 1/s²"
    @test string(ubase(1.0u"1/s^2"), pretty=false) == "(1.0)1/s^2"
    @test string(1.0*AffineUnits(dims=dimension(u"m/s^2")), pretty=true) == "1.0 m/s²"
    @test string(1.0*AffineUnits(dims=dimension(u"m/s^2")), pretty=false) == "(1.0)m/s^2"

    #Vector operations
    vq = quantity([1,2], u"m/s")
    @test [1,2].*u"m/s" == [1u"m/s", 2u"m/s"]
    @test [1,2].*(5u"m/s") == [5.0u"m/s", 10.0u"m/s"]
    @test vq[:] == quantity([1,2], u"m/s")
    @test vq[1] == 1*u"m/s"
    @test vq[CartesianIndex(1)] == 1*u"m/s"
    @test all([q for q in vq] .== vq)
    @test vq[begin] == 1*u"m/s"
    @test vq[end] == 2*u"m/s"

    #Size indicators for quantities
    tx = [1, 3.6, 501.3]
    tq = tx*u"kg/hr"
    @test size(tq) == size(tx)
    @test length(tq) == length(tx)
    @test axes(tq) == axes(tx)
    @test ndims(tq) == ndims(tx)
    @test ndims(typeof(tq)) == ndims(typeof(tx))
    @test Base.broadcastable(tq) == tq

    #Size indicators for units
    @test size(u"kW") == size(1)
    @test length(u"kW") == length(1)
    @test axes(u"kW") == axes(1)
    @test ndims(u"kW") == ndims(1)
    @test ndims(typeof(u"kW")) == ndims(Float64)
    @test iterate(u"kW") == (u"kW", nothing)
    @test (5u"kg/hr")[] == 5u"kg/hr"

end

@testset "Math on various types" begin
    FlexUnits.PRETTY_DIM_OUTPUT[] = true  

    for Q in [Quantity], T in [Float16, Float32, Float64], R in [DEFAULT_RATIONAL, Rational{Int16}, Rational{Int32}]
        
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

        @test typeof(y) <: Quantity
        @test typeof(x).parameters[1] == T
        @test typeof(x).parameters[2] == D
        @test ustrip(y) ≈ T(0.04)

        if R <: Rational  
            @test string(x, pretty=true) == "0.2 (m kg⁵ᐟ²)/s"
            @test string(inv(x), pretty=true) == "5.0 s/(m kg⁵ᐟ²)"
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
       
        u = Dimensions{R}(length=2//5)
        x = Quantity(-1.2, u)

        @test typemax(x) == Quantity(typemax(-1.2), u)
    
        @test abs(x) == Quantity(1.2, u)
        @test abs(x) == abs(Quantity(1.2, u))
        @test abs2(x) == Quantity(abs2(-1.2), u^2)
    
        @test deepcopy(x) == x
    
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
        @test zero(Dimensions{R}) === Dimensions{R}()
        @test zero(AffineUnits{Dimensions{R}}) === AffineUnits{Dimensions{R}}(dims=zero(Dimensions{R}))
        @test zero(Quantity{T, Dimensions{R}}) === Quantity(zero(T), zero(Dimensions{R}))
        @test zero(Quantity{T, AffineUnits{Dimensions{R}}}) == Quantity(zero(T), zero(AffineUnits{Dimensions{R}}))
        @test one(Quantity{T, AffineUnits{Dimensions{R}}}) === one(T)
        @test oneunit(Quantity{T, AffineUnits{Dimensions{R}}}) == Quantity(one(T), zero(AffineUnits{Dimensions{R}}))
        @test oneunit(Quantity{T, Dimensions{R}}) == Quantity(one(T), zero(Dimensions{R}))

        @test min(0.0u"kg/hr", 1.2u"lb/s") == 0u"kg/s"
        @test max(0.0u"kg/hr", 1.2u"lb/s") == 1.2u"lb/s"
        @test max(0.0u"°C", 0.0u"°F") == ubase(0.0u"°C")
        @test maximum([1.0u"lb/hr", 1.0u"kg/hr", 1.0u"kg/s"]) == 1.0u"kg/s"
        @test minimum([1.0u"lb/hr", 1.0u"kg/hr", 1.0u"kg/s"]) == 1.0u"lb/hr"
        @test sort([1.0u"kg/s", 1.0u"lb/hr", 1.0u"kg/hr"]) == [1.0u"lb/hr", 1.0u"kg/hr", 1.0u"kg/s"]

        #Cannot check iszero on non-affine units 
        @test_throws NotScalarError iszero(5u"°C")

    end

    #Other mathematical operators/functions
    @test sin(5u"rad") == sin(5) 
    @test cos(60u"deg") ≈ 0.5
    @test 5u"km/hr" < 6u"km/hr"
    @test 5u"m/s" > 5u"km/hr"
    @test 5u"km/hr" <= 6u"km/hr"
    @test 5u"m/s" >= 5u"km/hr"
    @test 5u"m/s" >= 5u"m/s"
    @test 5u"m/s" <= 5u"m/s"
    @test_throws DimensionError 5u"m/s" <= 5u"kg"

    #Math on arrays of number quantities 
    mq = [5*u"m/s" 2u"m/s^2"; 1*u"kg/s" 4*u"kg/s^2"]
    vq = [1u"s", 2u"s^2"]
    @test mq*vq == [9.0u"m", 9.0u"kg"]
    @test mq + 2*mq == 3*mq
    @test 3.0.*mq .- 2.0.*mq == mq
    @test vq.*vq == [1u"s^2", 4u"s^4"]

    #Math on arrays of quantities
    qm = [1 2; 3 4]*u"m/s"
    qi = inv(qm)
    @test qm*qi ≈ [1 0; 0 1]*u""
    @test sum(qm) ≈ 10u"m/s"
    @test qm[1] == 1u"m/s"

end

@testset "UnitfulCallable" begin
    #Test callable application
    angle_coords(θ::Real, r::Real) = r.*(cos(θ), sin(θ))
    unitful_angle_coords = UnitfulCallable(angle_coords, (u"", u"m") => (Dimensions(length=1), Dimensions(length=1)))
    c = unitful_angle_coords(30u"deg", 6u"cm")
    @test all(c .≈ (cosd(30), sind(30)).*(ubase(0.06u"m")))

    #Test unit_call application
    function angle_coords(θ::Quantity{<:Real}, r::Quantity{<:Real})
        u = UnitfulCallable( (u"", unit(r)) => unit(r) )
        return unitful_call(angle_coords, u, θ, r)
    end
    c = angle_coords(30u"deg", 100.0u"cm")
    @test all(c .≈ 1u"m".*(cosd(30), sind(30)))
    @test unit(c) == u"cm"

    #Test rotating arm application
    struct RotatingArm
        len :: Float64
    end
    angle_coords(arm::RotatingArm, θ::Real) = arm.len.*(cos(θ), sin(θ))
    angle_coords(arm::Quantity{<:RotatingArm}, θ::Quantity{<:Real}) = UnitfulCallable(Base.Fix1(angle_coords, ustrip(arm)), u""=>unit(arm))(θ)
    c = angle_coords(Quantity(RotatingArm(1.0), u"m"), 30u"deg")
    @test all(c .≈ 1u"m".*(cosd(30), sind(30)))


    #Test linear transfomration where x is assumed to be [u"kg/s", u"m^3/s"]
    heat_rate(x::AbstractVector) = [4.136,0.235]'*x
    heat_rate(x::AbstractVector{<:Quantity}) = UnitfulCallable(heat_rate, [u"kg/s", u"m^3/s"]=>u"kW")(x)
    @test heat_rate([1000u"g/s", 1000u"L/s"]) ≈ sum([4.136,0.235])*u"kW"

    #Test no-argument call
    unitless_const() = 5.0
    unitful_const() = unitful_call(unitless_const, UnitfulCallable(()=>u"m"))
    unitful_const()

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
    @test typeof(x) == typeof(q"1.3km/s^2")
    @test abs(x) === ubase(q"1.3km/s^2")

    y = 0.9u"sqrt(mΩ)"
    @test typeof(y) == Quantity{Float64, AffineUnits{Dimensions{DEFAULT_RATIONAL}}}
    @test typeof(ubase(y)) == Quantity{Float64, Dimensions{DEFAULT_RATIONAL}}
    @test ustrip_base(y) ≈ 0.02846049894151541
    @test y ≈ q"(0.9*sqrt(mΩ))"

    y = BigFloat(0.3) * u"mΩ"
    @test typeof(y) == Quantity{BigFloat, AffineUnits{Dimensions{DEFAULT_RATIONAL}}}
    @test ustrip_base(y) ≈ 0.0003

    y32 = convert(Quantity{Float32, AffineUnits{Dimensions{DEFAULT_RATIONAL}}}, y)
    @test typeof(y32) == Quantity{Float32, AffineUnits{Dimensions{DEFAULT_RATIONAL}}}

    z = 1.0*u"yr"
    @test ustrip_base(z) ≈ 60 * 60 * 24 * 365.25
    @test z === 1.0*uparse("yr")
    @test z === qparse("1yr")
    @test 1/z === ubase(qparse("1/yr"))
    @test_throws MethodError qparse("yr")

    # Test type stability of extreme range of units
    U = typeof(first(values(UnitRegistry.UNITS)))
    @test typeof(u"s"^2) == U
    @test typeof(u"Ω") == U

    @test_throws ArgumentError uparse(":x")
    @test_throws ArgumentError qparse(":x")
    @test_throws "Symbol x not found" uparse("x")
    @test_throws "Unexpected expression" uparse("import ..Units")
    @test_throws "Unexpected expression" uparse("(m, m)")
    @test_throws LoadError eval(:(us"x"))
    @test_throws NotScalarError u"J/°C"


    #Basic mathematical operations
    xp = 1u"percent"
    xv = 1u"m/s"

    @test +(xp, xp, xp) + 1 == 1.03
    @test xp + 1 == 1.01
    @test 1 + xp == 1.01
    @test xv + xv == 2u"m/s"
    @test (xv + xv) isa Quantity{Float64, <:Dimensions}

    @test 1 - xp == 0.99 
    @test xp - 1 == -0.99
    @test xv - xv == 0u"m/s"
    @test (xv - xv) isa Quantity{Float64, <:Dimensions}

    @test_throws DimensionError xv + 1
    @test_throws DimensionError 1 + xv

    @test xv*NoDims() == xv
    @test xv*1 == 1u"m/s"
    @test xv*1 !== 1u"m/s"
    @test xv*1 === ubase(1u"m/s")
    @test xp*1 == 0.01*u""
    @test xv*xv == 1u"m^2/s^2"
    @test *(xv,xv,xv) == 1u"m^3/s^3"
    @test cbrt(*(xv,xv,xv)) == 1u"m/s"
    @test u"m/s" / unit(ubase(1u"m/s")) == u""
    @test cbrt(u"m^3/s^3") == u"m/s"
    @test abs2(dimension(u"m/s")) == dimension(u"m^2/s^2")

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

    #Ambiguity tests 
    cplx  = 1.0 + im 
    cplxb = Complex{Bool}(true + im)
    @test xv*cplx  === ubase(Quantity(cplx, u"m/s"))
    @test xv*cplxb === ubase(Quantity(cplxb, u"m/s"))
    @test cplx*xv  === ubase(Quantity(cplx, u"m/s"))
    @test cplxb*xv === ubase(Quantity(cplxb, u"m/s"))

    @test xv/cplx  === ubase(Quantity(inv(cplx), u"m/s"))
    @test xv/cplxb === ubase(Quantity(inv(cplxb), u"m/s"))
    @test cplx/xv  === ubase(Quantity(cplx, inv(u"m/s")))
    @test cplxb/xv === ubase(Quantity(cplxb, inv(u"m/s")))

    @test cplx + xp === cplx + 0.01
    @test xp + cplx === cplx + 0.01
    @test cplxb + xp === 0.01 + cplxb
    @test xp + cplxb === 0.01 + cplxb

    @test cplx - xp === cplx - 0.01
    @test xp - cplx === 0.01 - cplx
    @test cplxb - xp === cplxb - 0.01
    @test xp - cplxb === 0.01 - cplxb

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
    @test 1°C |> K isa Quantity{<:Real, <:AffineUnits}
    @test 1°C |> unit(ubase(1K)) isa Quantity{<:Real, <:Dimensions}
    @test  °C |> unit(ubase(1K)) isa AffineTransform

    @test 0°C |> u"K" == 273.15u"K"
    @test 1u"K" |> °C isa Quantity{<:Real, <:AffineUnits}
    @test 0u"K" |> °C  == -273.15°C
    @test °C |> °F isa AffineTransform
    @test 0°C |> °F == 32°F
    @test (°C |> °F)(0) ≈ 32

    @test AffineUnits(dims=u"Pa") == u"Pa"
    @test_throws NotDimensionError AffineUnits(dims=u"kPa")

    # Test display against errors
    celsius = AffineUnits(offset=273.15, dims=u"K")
    psi = FlexUnits.asunit(6.89476u"kPa")
    io = IOBuffer()
    @test isnothing(show(io, (dimension(°F), dimension(u"K"), psi, celsius, fahrenheit)))

    # Test updating affine units
    @test register_unit("°C" => °C) isa AbstractDict # Updating the same value does nothing for unit
    @test register_unit("K" => Dimensions(temperature=1)) isa AbstractDict # same value yields nothing for dimension
    @test_throws MethodError register_unit(u"K") # cannot register only a unit

    # Cannot re-register a unit if its value changes nor can we delete a unit
    @test_throws RegistryTools.PermanentDictError register_unit("°C"=>u"°F")
    @test_throws RegistryTools.PermanentDictError delete!(UnitRegistry.UNITS, :m)

    # Cannot register non-parsable unit
    @test_throws ArgumentError register_unit("m/s"=>u"m/s")

    # Test map_dimensions
    @test map_dimensions(+, dimension(u"m/s"), dimension(u"m/s")) == Dimensions(length=2, time=-2)
    @test map_dimensions(-, dimension(u"m"), dimension(u"s"))     == Dimensions(length=1, time=-1)
    @test map_dimensions(Base.Fix1(*,2), dimension(u"m/s"))       == Dimensions(length=2, time=-2)

    # Parsing tests 
    @test u"1.0" == u""
    @test q"5 kg/s" == 5u"kg/s"
    @test uparse("1.0") == u""
    @test qparse("10 km/hr") == 10*u"km/hr"
    @test u"%" == u"percent"
    @test q"1%" == 0.01u""
    @test q"1.0" == 1.0u""
    @test qparse("1%") == 0.01u""
    @test qparse("1.0") == 1.0u""
    @test qparse("5  5kg") == 25u"kg"

    #Test showerror 
    testio = IOBuffer()

    showerror(testio, ConversionError(u"m/s", u"m"))
    @test String(take!(testio)) == "ConversionError: Cannot convert unit 'm' to target unit 'm/s'. Consider multiplying 'm/s' by 's' or similar."

    showerror(testio, DimensionError(dimension(u"m/s")))
    @test String(take!(testio)) == "DimensionError: m/s is not dimensionless"

    showerror(testio, DimensionError(1u"m/s"))
    @test String(take!(testio)) == "DimensionError: 1 m/s is not dimensionless"

    showerror(testio, DimensionError((dimension(u"m/s"), dimension(u"s/m"))))
    @test String(take!(testio)) == "DimensionError: (m/s, s/m) have incompatible dimensions"

    showerror(testio, NotScalarError(u"°C"))
    @test String(take!(testio)) == "NotScalarError: °C cannot be treated as scalar, operation only valid for scalar units"

    showerror(testio, NotDimensionError(u"kPa"))
    @test String(take!(testio)) == "NotDimensionError: kPa cannot be treated as dimension, operation only valid for dimension units"
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
    @test dimension(q32_32) == convert(typeof(dimension(q32_32)), dimension(q))
    @test typeof(convert(Quantity{Float16}, q)) == Quantity{Float16,Dimensions{Rational{Int16}}}
    @test convert(Quantity, q) === q
    @test convert(Quantity{Float64, AffineUnits{DEFAULT_DIM_TYPE}}, ubase(1u"kg")) isa Quantity{Float64, AffineUnits{DEFAULT_DIM_TYPE}}

    # Test that regular type promotion applies:
    q = Quantity(2, d)
    @test typeof(q) == Quantity{Int64,typeof(d)}
    @test typeof(q ^ 2) == Quantity{Int64,typeof(d)}
    @test typeof(0.5 * q) == Quantity{Float64,typeof(d)}
    @test typeof(inv(q)) == Quantity{Float64,typeof(d)}

    # Test conversion of unit types 
    @test convert(DEFAULT_DIM_TYPE, u"m") === Dimensions(length=1)
    @test_throws NotDimensionError convert(DEFAULT_DIM_TYPE, u"mm")
    @test convert(DEFAULT_DIM_TYPE, NoDims()) === Dimensions(length=0)
    @test convert(AffineUnits{DEFAULT_DIM_TYPE}, ubase(2u"m")) == AffineUnits(scale=2.0, offset=0.0, dims=dimension(u"m"))
    @test convert(AffineUnits{DEFAULT_DIM_TYPE}, 2u"°C") == AffineUnits(scale=2.0, offset=273.15, dims=dimension(u"K"))
    @test convert(Quantity{Float64, DEFAULT_DIM_TYPE}, 2u"m") === Quantity{Float64, DEFAULT_DIM_TYPE}(2.0, dimension(u"m")) 
    @test_throws NotScalarError convert(Quantity{Float64, DEFAULT_DIM_TYPE}, u"°C") 
    @test promote_type(DEFAULT_DIM_TYPE, NoDims{Int64}) == DEFAULT_DIM_TYPE
    @test promote_type(Quantity{Float32, DEFAULT_DIM_TYPE}, Quantity{Float64, DEFAULT_UNIT_TYPE}) == Quantity{Float64, DEFAULT_DIM_TYPE}

    # Test that adding different dimension subtypes still works
    @test 1*Dimensions{Int64}(length=1) + 1u"m" == 2u"m"

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
    @test typeof(uconvert(u"nm", 5e-9u"m")) <: Quantity{Float64, U} 
    @test typeof(uconvert(u"nm", Quantity(5e-9, u"m"))) <: Quantity{Float64, U}
    @test uconvert(u"nm", Quantity(5e-9, u"m")) ≈ 5u"nm"

    @test dimension(1u"m" |> u"nm")[:length] == 1

   
    # Different types require converting both arguments:
    q = convert(Quantity{Float16}, 1.5u"g")

    # Broadcasting conversions over Arrays
    x = [1.0, 2.0, 3.0] .* Quantity(1, u"kg")
    x2 = x .|> u"g"
    @test typeof(x2) <: Vector{<:Quantity{Float64,<:AffineUnits{<:Any}}}
    @test x2[2] ≈ (2000u"g")
end

@testset "Additional tests of FixedRational" begin
    #Tests were basically copied from DynamicQuantities, since FixedRational isn't its own package
    @test convert(Int64, FixedRational{1000,Int64}(2 // 1)) == 2
    @test convert(Int32, FixedRational{1000,Int64}(3 // 1)) == 3
    @test convert(Bool, FixedRational{6,Int8}(1//1)) === true
    @test convert(Bool, FixedRational{6,Int8}(0//1)) === false
    @test convert(FixedRational{100, Int32}, FixedRational{1000}(10)) === FixedRational{100, Int32}(10) 

    #Test promotion through mathematical operators 
    @test FixedRational{1000}(10)*Float64(0.1) === Float64(1)
    @test FixedRational{1000}(10)*Float32(0.1) === Float32(1)
    @test FixedRational{1000}(10)*BigFloat(0.1) isa BigFloat
    @test FixedRational{1000}(10)*true === FixedRational{1000}(10)

    #Test mathematical operators
    @test FixedRational(0.1)/Float64(100) === Float64(0.001)
    @test FixedRational(0.1) / FixedRational(0.1) === FixedRational(1.0)
    @test FixedRational(0.1) - FixedRational(0.1) === FixedRational(0.0)
    @test FixedRational(0.1) - 0.1 === 0.0
    @test FixedRational(0.1) + FixedRational(0.1) === FixedRational(0.2)
    @test FixedRational(0.1) + 0.1 === 0.2
    @test inv(FixedRational(0.1)) === FixedRational(10.0)

    #Test Rounding 
    @test round(Int32, FixedRational(10.2)) === Int32(10)
    @test round(Missing, FixedRational(10.2)) === missing

    #Test string conversion 
    @test string(FixedRational(10.2)) == "51//5"
    @test string(FixedRational(10)) == "10"

    #Test conversions and rounding
    @test_throws InexactError convert(Int32, FixedRational{6,Int8}(2//3))
    @test_throws InexactError convert(Bool, FixedRational{6,Int8}(2//1))
    @test_throws "Refusing to" promote(FixedRational{10,Int}(2), FixedRational{4,Int}(2))

    f64 = FixedRational{10,Int}(2)
    f8  = FixedRational{10,Int8}(2)
    @test promote(f64, f8) == (2, 2)
    @test typeof(promote(f64, f8)) == typeof((f64, f64))
    @test typeof(promote(FixedRational{10,Int8}(2), FixedRational{10,Int8}(2))) == typeof((f8, f8))

    # Required to hit integer branch (otherwise will go to `literal_pow`)
    f(i::Int) = Dimensions(length=1, mass=-1)^i
    @test f(2) == Dimensions(length=2, mass=-2)

    # Null conversion
    @test typeof(FixedRational{10,Int}(FixedRational{10,Int}(2))) == FixedRational{10,Int}

    # Conversion to Rational without specifying type
    @test convert(Rational, FixedRational{6,UInt8}(2)) === Rational{UInt8}(2)

    # Promotion rules
    @test promote_type(Float64, typeof(f64)) == Float64
    @test promote_type(Bool, typeof(f64)) == typeof(f64)
    @test promote_type(FixedRational{10,Int64},FixedRational{10,BigInt}) == FixedRational{10,BigInt}
    @test promote_type(Rational{Int8}, FixedRational{12345,Int}) == Rational{Int}
    @test promote_type(Int8, FixedRational{12345,Int}) == FixedRational{12345,Int}
    @test promote_type(FixedRational{1000}, Irrational) == Real
    @test promote_type(FixedRational{1000}, BigFloat) == BigFloat

    # Bug where user would create a FixedRational{::Type{Int32}, ::Int64} and get stack overflow,
    # because the stored type was FixedRational{::Type{Int32}, ::Int32}
    x = 10u"m"
    user_quantity = Quantity(10.0, Dimensions{FixedRational{25200,Int32}}(1, 0, 0, 0, 0, 0, 0))
    @test x == user_quantity
end

#Register a new affine unit (and verify re-registering)
register_unit("psig" => AffineUnits(scale=uscale(u"psi"), offset=101.3u"kPa", dims=u"Pa"))

@testset "Registration tests" begin
    #Test re-registering and verify that the unit exists
    register_unit("psig" => AffineUnits(scale=uscale(u"psi"), offset=101.3u"kPa", dims=u"Pa"))
    @test 0*u"psig" == 101.3u"kPa"
    @test q"0psig" == 101.3u"kPa"
    
    #Test registration for a different registry base type
    IntDimType = Dimensions{Int32}
    reg = RegistryTools.PermanentDict{Symbol, AffineUnits{IntDimType}}()
    reg = RegistryTools.registry_defaults!(reg)  
    @test reg[:m]  === AffineUnits(dims=IntDimType(length=1), symbol=:m)
    @test reg[:kg] === AffineUnits(dims=IntDimType(mass=1), symbol=:kg)    
end



@testset "Aqua.jl" begin
    Aqua.test_all(FlexUnits)
end

R = Ref(8.314 * u"J/(mol*K)")
@benchmark let R=R[]
    v_satp = R * (25u"°C") / (101.3u"kPa")
end

#Profiling diagnostic if benchmarks are off 
#=
v = randn(100000).*u""
v2 = deepcopy(v)
@profview for ii in 1:1000
    v2 .= v.*(1.1*u"m/s")
end
=#