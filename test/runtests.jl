#============================================================================================================================
Run these commands at startup to see coverage
julia --startup-file=no --depwarn=yes --threads=auto -e 'using Coverage; clean_folder(\"src\"); clean_folder(\"test\"); clean_folder(\"ext\") '
julia --startup-file=no --depwarn=yes --threads=auto --code-coverage=user --project=. -e 'using Pkg; Pkg.test(coverage=true)'
julia --startup-file=no --depwarn=yes --threads=auto coverage.jl

Run this command for testing invalidations
julia --startup-file=no --depwarn=yes --threads=auto --project=. test/invalidations.jl
============================================================================================================================#

#To see the actual coverage in VSCode, install the Coverage Gutters extension
#https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters

using Revise
using FlexUnits, .UnitRegistry
using TestItems: @testitem
using TestItemRunner
using Test
using BenchmarkTools
using FlexUnits: DEFAULT_RATIONAL, FixedRational, map_dimensions, dimval, FixRat64
using Aqua
using LinearAlgebra
using Statistics
using StaticArrays
import Unitful

const DEFAULT_UNIT_TYPE = typeof(first(values(UnitRegistry.UNITS)))
const DEFAULT_DIM_TYPE  = FlexUnits.dimtype(DEFAULT_UNIT_TYPE)
const AT = AffineTransform{Float64}

@testset "Basic utilities" begin
    #Basic Type operations
    @test UnitRegistry.unittype() === DEFAULT_UNIT_TYPE
    @test UnitRegistry.dimtype() === DEFAULT_DIM_TYPE

    d = Dimensions(length=2, mass=1, time=-2)
    D = typeof(d)
    @test d.length == 2
    @test d.mass == 1
    @test d.time == -2
    @test FlexUnits.dimtype(typeof(d)) == DEFAULT_DIM_TYPE
    @test FlexUnits.dimtype(typeof(u"m/s")) == DEFAULT_DIM_TYPE
    @test uscale(d) == 1
    @test uoffset(d) == 0
    @test FlexUnits.usymbol(d) == FlexUnits.DEFAULT_USYMBOL
    @test eltype(D) == FixRat32
    @test FlexUnits.dimpowtype(D) == FixRat32
    @test FlexUnits.dimpowtype(typeof(ud"m/s")) == FixRat32
    @test FlexUnits.is_dimension(u"m/s")
    @test FlexUnits.is_dimension(dimension(u"m/s"))
    @test FlexUnits.is_dimension(ud"m/s")
    @test !FlexUnits.is_dimension(u"kPa")
    @test FlexUnits.is_scalar(u"m/s")
    @test FlexUnits.is_scalar(dimension(u"m/s"))
    @test FlexUnits.is_scalar(ud"m/s")
    @test !FlexUnits.is_scalar(u"°F")
    @test udynamic(u"m/s") === ud"m/s"
    @test udynamic(dimension(u"m/s")) === dimension(ud"m/s")

    @test FlexUnits.remove_offset(u"°C") == u"K"
    @test FlexUnits.constructorof(typeof(Dimensions())) == Dimensions
    @test FlexUnits.constructorof(typeof(u"m")) == StaticUnits
    @test FlexUnits.constructorof(typeof(1.0*u"m")) == Quantity
    @test FlexUnits.constructorof(typeof((1.0+im)*u"m")) == Quantity
    @test FlexUnits.constructorof(typeof(Quantity("this", ud"m"))) == Quantity
    @test FlexUnits.constructorof(Array{Float64}) == Array

    t_none = NoTransform()
    @test t_none(1) === 1
    @test t_none(AffineTransform()) === AffineTransform()
    @test uscale(t_none) == 1
    @test uoffset(t_none) == 0
    @test FlexUnits.is_identity(t_none)
    @test FlexUnits.is_scalar(t_none)

    @test StaticDims(dimension(ud"m/s")) == StaticDims{dimension(ud"m/s")}()
    @test_throws ArgumentError StaticDims{dimension(ud"m/s")}(dimension(ud"kg/s"))
    
    t_affine = AffineTransform()
    t_affine2 = AffineTransform(scale=2.0, offset=1.0)
    @test t_affine((1,2,3)) == (1.0, 2.0, 3.0)
    @test t_affine2((1,2,3)) == (3.0, 5.0, 7.0)
    @test t_affine(t_affine2) == t_affine2
    @test t_affine2(t_none) == t_affine2
    @test uscale(t_affine2) == 2
    @test uoffset(t_affine2) == 1

    @test Units{Dimensions{FixRat32}}(Dimensions{FixRat64}(length=1), AffineTransform()) == ud"m" 
    @test Units(ud"m/s") == ud"m/s"
    @test Units(u"m/s") == ud"m/s"
    @test FlexUnits.dimtype(Units{Dimensions{FixRat32}}(Dimensions{FixRat64}(length=1), AffineTransform())) <: Dimensions{FixRat32}
    @test FlexUnits.dimtype(Units(Dimensions{FixRat64}(length=1), AffineTransform())) <: Dimensions{FixRat64}
    @test FlexUnits.dimtype(u"m/s") == typeof(dimension(ud"m/s"))
    @test StaticUnits(dimension(ud"m/s"), AffineTransform()) == u"m/s"

    @test ustrip(Quantity{Float32}(1.0, u"m/s")) === Float32(1)
    @test FlexUnits.unittype(Quantity{Float64, typeof(u"m/s")}) == typeof(u"m/s")
    @test FlexUnits.dimtype(Quantity{Float64, typeof(u"m/s")}) == Dimensions{FixRat32}
    @test FlexUnits.dimvaltype(1.0u"kJ") == Dimensions{FixRat32}
    @test FlexUnits.dimtype(1.0u"m") == StaticDims{Dimensions(length=1)}
    @test udynamic(1.0u"kJ") === ubase(1.0ud"kJ")
    @test assert_dimension(dimension(u"m/s")) == dimension(u"m/s")

    #=
    @test MirrorDims() == MirrorDims{Dimensions{FixRat32}}()
    @test eltype([MirrorDims(), Dimensions{FixRat32}()]) == MirrorUnion{Dimensions{FixRat32}}
    @test_throws "MirrorDims should not be a type parameter" Quantity{Float64, MirrorDims{Dimensions{FixRat32}}}(1, u"m/s")
    =#
    
    #Basic utility functions
    vsum = sum(ones(5,3).*[u"rpm" u"kg/hr" u"kPa"], dims=1)

    FlexUnits.pretty_print_units(false)
    @test string(1.0u"kg*m^2/s^2") == "(1.0)(m^2*kg)/s^2"
    @test string(1.0u"1/s^2") == "(1.0)1/s^2"
    @test string(1.0*(ud"km/hr"*ud"m/s")) == "(0.2777777777777778)m^2/s^2"
    @test string(u"m/s"*u"m/s") == "StaticUnits{m^2/s^2, AffineTransform{Float64}}(AffineTransform{Float64}(1.0, 0.0), :_)"
    #@test string(MirrorDims{Dimensions{FixRat32}}()) == "MirrorDims{Dimensions{FixRat32}}()"
    @test string(zero(typeof(1ud"m/s"))) == "(0.0)?/?"
    @test string(vsum) == "Quantity{Float64, Dimensions{FixRat32}}[(0.5235987755982988)1/s (0.001388888888888889)kg/s (5000.0)kg/(m*s^2)]"
    
    FlexUnits.pretty_print_units(true)
    @test string(1.0u"kg*m^2/s^2")  == "1.0 (m² kg)/s²"
    @test string(1.0u"1/s^2")  == "1.0 1/s²"
    @test string(1.0*(ud"km/hr"*ud"m/s")) == "0.2777777777777778 m²/s²"
    @test string(u"m/s"*u"m/s") == "StaticUnits{m²/s², AffineTransform{Float64}}(AffineTransform{Float64}(1.0, 0.0), :_)"
    @test string(zero(typeof(1ud"m/s"))) == "0.0 ?/?"
    @test string(vsum) == "Quantity{Float64, Dimensions{FixRat32}}[0.5235987755982988 1/s 0.001388888888888889 kg/s 5000.0 kg/(m s²)]"

    @test string(FlexUnits.unit_symbols(Dimensions{FixRat32})) == "(:length => :m, :mass => :kg, :time => :s, :current => :A, :temperature => :K, :luminosity => :cd, :amount => :mol)"

    #Vector operations RETRY LATER
    #=
    vq = Quantity([1,2], u"m/s")
    @test [1,2].*u"m/s" == [1u"m/s", 2u"m/s"]
    @test [1,2].*(5u"m/s") == [5.0u"m/s", 10.0u"m/s"]
    @test vq[:] == Quantity([1,2], u"m/s")
    @test vq[1] == 1*u"m/s"
    @test vq[CartesianIndex(1)] == 1*u"m/s"
    @test all([q for q in vq] .== vq)
    @test vq[begin] == 1*u"m/s"
    @test vq[end] == 2*u"m/s"
    =#

    #Size indicators for quantities
    tx = [1, 3.6, 501.3]
    tq = tx.*u"kg/hr"
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

        FlexUnits.pretty_print_units(true)
        if R <: Rational  
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
        @test zero(Units{Dimensions{R}, AT}) === Units{Dimensions{R}, AT}(dims=zero(Dimensions{R}), todims=AT())

        #Static identity transform tests
        @test zero(1u"m/s") + 2.0u"m/s" == 2.0u"m/s"
        @test zero(typeof(1u"m/s")) + 2.0u"m/s" == 2.0u"m/s"
        @test zero(typeof(ubase(1u"m/s"))) + 2.0u"m/s" == ubase(2.0u"m/s")
        @test zero(zero(typeof(ubase(1u"m/s")))) + 2.0u"m/s" == ubase(2.0u"m/s")
        @test zero(typeof(zero(typeof(ubase(1u"m/s"))))) + 2.0u"m/s" == ubase(2.0u"m/s")
        @test max(typemin(1.0u"m/s"), 0.0u"m/s") == 0.0u"m/s"
        @test max(typemin(typeof(1.0u"m/s")), 0.0u"m/s") == 0.0u"m/s"
        @test min(typemax(-1.0u"m/s"), 0.0u"m/s") == 0.0u"m/s"
        @test min(typemax(typeof(1.0u"m/s")), 0.0u"m/s") == 0.0u"m/s"

        @test one(Quantity{T, Units{Dimensions{R}, AT}}) === one(T)
        @test oneunit(typeof(5.0u"m/s")) == 1.0u"m/s"
        @test_throws ArgumentError oneunit(Quantity{T, Units{Dimensions{R}, AT}})

        @test min(0.0u"kg/hr", 1.2u"lb/s") == 0u"kg/s"
        @test max(0.0u"kg/hr", 1.2u"lb/s") == 1.2u"lb/s"
        @test max(0.0u"°C", 0.0u"°F") == ubase(0.0u"°C")
        @test maximum([1.0u"lb/hr", 1.0u"kg/hr", 1.0u"kg/s"]) == 1.0u"kg/s"
        @test minimum([1.0u"lb/hr", 1.0u"kg/hr", 1.0u"kg/s"]) == 1.0u"lb/hr"
        @test sort([1.0u"kg/s", 1.0u"lb/hr", 1.0u"kg/hr"]) == [1.0u"lb/hr", 1.0u"kg/hr", 1.0u"kg/s"]

        #Dynamic identity transform tests
        @test zero(1ud"m/s") + 2.0ud"m/s" == 2.0ud"m/s"
        @test zero(typeof(1ud"m/s")) + 2.0ud"m/s" == 2.0ud"m/s"
        @test zero(typeof(ubase(1ud"m/s"))) + 2.0ud"m/s" == ubase(2.0ud"m/s")
        @test zero(zero(typeof(ubase(1ud"m/s")))) + 2.0ud"m/s" == ubase(2.0ud"m/s")
        @test zero(typeof(zero(typeof(ubase(1ud"m/s"))))) + 2.0ud"m/s" == ubase(2.0ud"m/s")
        @test max(typemin(1.0ud"m/s"), 0.0ud"m/s") == 0.0ud"m/s"
        @test max(typemin(typeof(1.0ud"m/s")), 0.0ud"m/s") == 0.0ud"m/s"
        @test min(typemax(-1.0ud"m/s"), 0.0ud"m/s") == 0.0ud"m/s"
        @test min(typemax(typeof(1.0ud"m/s")), 0.0ud"m/s") == 0.0ud"m/s"

        @test one(Quantity{T, Units{Dimensions{R}, AT}}) === one(T)
        @test_throws ArgumentError oneunit(Quantity{T, Units{Dimensions{R}, AT}})
        @test oneunit(Quantity{T, Dimensions{R}}) === Quantity(one(T), FlexUnits.unknown(Dimensions{R}))

        @test min(0.0ud"kg/hr", 1.2ud"lb/s") == 0ud"kg/s"
        @test max(0.0ud"kg/hr", 1.2ud"lb/s") == 1.2ud"lb/s"
        @test max(0.0ud"°C", 0.0ud"°F") == ubase(0.0ud"°C")
        @test maximum([1.0ud"lb/hr", 1.0ud"kg/hr", 1.0ud"kg/s"]) == 1.0ud"kg/s"
        @test minimum([1.0ud"lb/hr", 1.0ud"kg/hr", 1.0ud"kg/s"]) == 1.0ud"lb/hr"
        @test sort([1.0ud"kg/s", 1.0ud"lb/hr", 1.0ud"kg/hr"]) == [1.0ud"lb/hr", 1.0ud"kg/hr", 1.0ud"kg/s"]

        #Cannot check iszero on non-affine units 
        @test_throws NotScalarError iszero(5u"K" |> u"°C")
        @test_throws NotScalarError iszero(5ud"K" |> ud"°C")
    end

    #Other mathematical operators/functions
    +(5u"m") == 5u"m"
    @test sin(5u"rad") ≈ sin(5)
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
    q0 = [1 2; 3 4]
    qm = q0 .* u"m/s"
    @test sum(qm) ≈ 10u"m/s"
    @test qm[1] == 1u"m/s"

    #Other math on mirror dimensions
    u_unk = FlexUnits.unknown(FlexUnits.dimtype(ud""))
    @test +(zero(typeof(1ud""))) == zero(typeof(1ud""))
    @test zero(typeof(1ud""))^2 == zero(typeof(1ud""))
    @test zero(typeof(1ud"")) + (1ud"m/s") == 1.0ud"m/s"
    @test (1ud"m/s") + zero(typeof(1ud"")) == 1.0ud"m/s"
    @test zero(typeof(1ud""))*(1ud"m/s") == 0.0u_unk
    @test (1ud"m/s")*zero(typeof(1ud"")) == 0.0u_unk
    @test zero(typeof(1ud""))*zero(typeof(1ud"")) == zero(typeof(1ud""))
    @test zero(typeof(1ud""))/(1ud"m/s") == 0.0u_unk
    @test (1ud"m/s")/zero(typeof(1ud"")) == Inf*u_unk
    @test unit(zero(typeof(1ud""))/zero(typeof(1ud""))) == unit(zero(typeof(1ud"")))
    @test zero(Quantity{Float64,Dimensions{FixRat32}}) == zero(typeof(1ud""))

    #Math on unit transforms 
    @test AffineTransform(scale=2, offset=0)*2 == AffineTransform(scale=4, offset=0)
    @test_throws ArgumentError AffineTransform(scale=2, offset=1)*2
    @test AffineTransform(scale=2, offset=0)/2 == AffineTransform(scale=1, offset=0)
    @test_throws ArgumentError AffineTransform(scale=2, offset=1)*2
    @test NoTransform()^60 == NoTransform()
    @test AffineTransform()^60 == AffineTransform()
    @test_throws ArgumentError AffineTransform(scale=2, offset=1)^2
    @test NoTransform()/2 == AffineTransform(scale=0.5, offset=0)

    #Math on units 
    @test ud"m/s"*ud"m/s" == ud"m^2/s^2"
    @test u"m/s"*u"m/s" == u"m^2/s^2"
    @test ud"m/s"/ud"kg/s" == ud"m/kg"
    @test u"m/s"/u"kg/s" == u"m/kg"
    @test (1.0ud"m/s")*u"m/s" == 1.0ud"m^2/s^2"
    @test ud"m/s"*(1.0ud"m/s") == 1.0ud"m^2/s^2"
    @test (1.0ud"m/s")/ud"kg/s" == 1.0ud"m/kg"
    @test ud"m/s"/(1.0ud"kg/s") == 1.0ud"m/kg"

    #Comparisons and math on missing
    @test 1.0u"" ≈ 1.0
    @test 1.0 ≈ 1.0u""
    @test ismissing(missing ≈ 1.0u"kg")
    @test ismissing(1.0u"kg" ≈ missing)
    @test ismissing(1u"m/s" + missing)
    @test ismissing(missing + 2u"kPa")
    @test ismissing(1u"m/s" - missing)
    @test ismissing(missing - 2u"kPa")
    @test ismissing(1u"m/s" * missing)
    @test ismissing(missing * 2u"kPa")
    @test ismissing(1u"m/s" / missing)
    @test ismissing(missing / 2u"kPa")
end

#=
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
=#

@testset "Dynamic unit functionality" begin
    x = Quantity(0.2, Dimensions(length=1, mass=2.5, time=-1))
    u = ustrip(x)

    #Test round-trip parsing
    @test ubase(x) ≈ ustrip(x)*uparse(FlexUnits.ustring(unit(x), pretty=false))

    #Test toggling pretty-print
    pretty_print_units(false)
    xp = ustrip(x)*uparse(string(unit(x)))
    pretty_print_units(true)
    @test ubase(xp) ≈ ubase(x)

    #Ustrip tests
    x = quantity(1.3, ud"km/s^2")
    @test ustrip(x) == 1.3
    @test dstrip(x) == 1300  # SI base units
    @test x == q"1.3km/s^2"
    @test typeof(x) == typeof(q"1.3km/s^2")
    @test abs(x) === ubase(q"1.3km/s^2")

    y = quantity(0.9, ud"sqrt(mΩ)")
    @test typeof(y) == Quantity{Float64, Units{Dimensions{DEFAULT_RATIONAL}, AT}}
    @test typeof(ubase(y)) == Quantity{Float64, Dimensions{DEFAULT_RATIONAL}}
    @test dstrip(y) ≈ 0.02846049894151541
    @test y ≈ q"(0.9*sqrt(mΩ))"

    y = quantity(BigFloat(0.3), ud"mΩ")
    @test typeof(y) == Quantity{BigFloat, Units{Dimensions{DEFAULT_RATIONAL}, AT}}
    @test dstrip(y) ≈ 0.0003

    y32 = convert(Quantity{Float32, Units{Dimensions{DEFAULT_RATIONAL}, AT}}, y)
    @test typeof(y32) == Quantity{Float32, Units{Dimensions{DEFAULT_RATIONAL}, AT}}

    z = quantity(1.0, ud"yr")
    @test dstrip(z) ≈ 60 * 60 * 24 * 365.25
    @test z === quantity(1.0, uparse("yr"))
    @test z === qparse("1yr")
    @test 1/z === ubase(qparse("1/yr"))
    @test ubase(z) === ubase(qparse("1*yr"))
    @test qparse("yr") == 1*ud"yr"

    # Test type stability of extreme range of units
    U = typeof(first(values(UnitRegistry.UNITS)))
    @test typeof(ud"s"^2) == U
    @test typeof(ud"Ω") == U

    @test_throws ArgumentError uparse(":x")
    @test_throws ArgumentError qparse(":x")
    @test_throws "Symbol x not found" uparse("x")
    @test_throws "Unexpected expression" uparse("import ..Units")
    @test_throws "Unexpected expression" uparse("(m, m)")
    @test_throws LoadError eval(:(us"x"))
    @test_throws NotScalarError ud"J/°C"

    #Basic mathematical operations
    xp = 1ud"percent"
    xv = 1ud"m/s"

    @test +(xp, xp, xp) + 1 == 1.03
    @test xp + 1 == 1.01
    @test 1 + xp == 1.01
    @test xv + xv == 2ud"m/s"
    @test (xv + xv) isa Quantity{Float64, <:Dimensions}

    @test 1 - xp == 0.99 
    @test xp - 1 == -0.99
    @test xv - xv == 0ud"m/s"
    @test (xv - xv) isa Quantity{Float64, <:Dimensions}

    @test_throws DimensionError xv + 1
    @test_throws DimensionError 1 + xv

    @test xv*1 == 1ud"m/s"
    @test xv*1 === 1ud"m/s"
    @test xv*1 === ubase(1ud"m/s")
    @test xp*1 == 0.01*ud""
    @test xv*xv == 1ud"m^2/s^2"
    @test *(xv,xv,xv) == 1ud"m^3/s^3"
    @test cbrt(*(xv,xv,xv)) == 1ud"m/s"
    @test ud"m/s" / unit(ubase(1ud"m/s")) == ud""
    @test cbrt(ud"m^3/s^3") == ud"m/s"
    @test abs2(dimension(ud"m/s")) == dimension(ud"m^2/s^2")

    @test xv/1 == 1ud"m/s"
    @test xp/1 == 0.01*ud""
    @test xv/xv == 1ud""

    @test sqrt(4*xv) == 2ud"m^0.5/s^0.5"
    @test (4*xv)^0.5 == 2ud"m^0.5/s^0.5"
    @test (2*xv)^2 == 4ud"m^2/s^2"

    @test (1.0ud"m")^0 == 1.0ud""
    @test (1.0ud"m")^1 == 1.0ud"m"
    @test (1.0ud"m")^2 == 1.0ud"m^2"
    @test (1.0ud"m")^3 == 1.0ud"m^3"
    @test (1.0ud"m")^4 == 1.0ud"m^4"
    @test (1.0ud"m")^(-1) == 1.0ud"1/m"
    @test (1.0ud"m")^(-2) == 1.0ud"m^(-2)"
    @test (1.0ud"m")^(-3) == 1.0ud"1/m^3"

    @test_throws "^ not defined" 2.0^(1ud"m/s")
    @test_throws "^ not defined" (1ud"m/s")^(1ud"m/s")

    @test_throws ArgumentError uparse("s[1]")
    @test_throws ArgumentError uparse("pounds_per_hour")


    #Tests on truly affine units
    °C  = ud"°C"
    °F  = ud"°F"
    K   = ud"K"
    @test dimension(°C).temperature == 1
    @test dimension(°C).length == 0
    @test °C == ud"degC"
    @test °F == ud"degF"
    @test 5°C - 4°C == 1K
    @test 5°C + 4°C == 282.15°C
    @test 5°C + 4°C - 0°C == 9°C
    @test qparse("0°C") == 0°C
    @test 0*uparse("°C") == 0°C

    #Ambiguity tests 
    cplx  = 1.0 + im 
    cplxb = Complex{Bool}(true + im)
    @test xv*cplx  === ubase(Quantity(cplx, ud"m/s"))
    @test xv*cplxb === ubase(Quantity(cplxb, ud"m/s"))
    @test cplx*xv  === ubase(Quantity(cplx, ud"m/s"))
    @test cplxb*xv === ubase(Quantity(cplxb, ud"m/s"))

    @test xv/cplx  === ubase(Quantity(inv(cplx), ud"m/s"))
    @test xv/cplxb === ubase(Quantity(inv(cplxb), ud"m/s"))
    @test cplx/xv  === ubase(Quantity(cplx, inv(ud"m/s")))
    @test cplxb/xv === ubase(Quantity(cplxb, inv(ud"m/s")))

    @test cplx + xp === cplx + 0.01
    @test xp + cplx === cplx + 0.01
    @test cplxb + xp === 0.01 + cplxb
    @test xp + cplxb === 0.01 + cplxb

    @test cplx - xp === cplx - 0.01
    @test xp - cplx === 0.01 - cplx
    @test cplxb - xp === cplxb - 0.01
    @test xp - cplxb === 0.01 - cplxb

    # Constructors
    kelvin = Units(dims=ud"K", todims=AffineTransform())
    @test 1*kelvin == 1K

    rankine = Units(todims=AffineTransform(scale=5/9, offset=0.0), dims=K)
    @test 1*rankine == (5/9)K

    fahrenheit = Units(todims=AffineTransform(scale=5/9, offset=(273.15-32*5/9)), dims=K)
    @test 1*fahrenheit ≈ 1°F

    celsius = Units(todims=AffineTransform(offset=273.15), dims=K)
    @test 1*celsius ≈ 1°C

    # Round-trip sanity checks
    @test -40°C ≈ -40°F
    @test -40.0*celsius ≈ -40.0*fahrenheit

    # Test promotion explicitly for coverage:
    @test eltype([1ud"km/hr", 1.0ud"km/hr"]) <: Quantity{Float64, Dimensions{FixRat32}}
    @test eltype([Units(dimval(u"m/s")), ud"m/s"]) <: Units{Dimensions{FixRat32}, AT}
    @test promote_type(Units{Dimensions{Int16}, AT}, Units{Dimensions{Int32}, AT}) === Units{Dimensions{Int32}, AT}
    @test promote_type(Dimensions{Int16}, Units{Dimensions{Int32}, AT}) === Units{Dimensions{Int32}, AT}
    
    # Test conversions
    @test 1°C |> K isa Quantity{<:Real, <:Units}
    @test 1°C |> unit(ubase(1K)) isa Quantity{<:Real, <:Dimensions}
    @test  °C |> unit(ubase(1K)) isa AffineTransform

    @test 0°C |> ud"K" == 273.15ud"K"
    @test 1ud"K" |> °C isa Quantity{<:Real, <:Units}
    @test 0ud"K" |> °C  == -273.15°C
    @test °C |> °F isa AffineTransform
    @test 0°C |> °F == 32°F
    @test (°C |> °F)(0) ≈ 32

    @test Units(dims=ud"Pa", todims=AffineTransform()) == ud"Pa"
    @test_throws NotDimensionError Units(dims=ud"kPa", todims=AffineTransform())

    # Test display against errors
    celsius = Units(todims=AffineTransform(offset=273.15), dims=ud"K")
    psi = Units(6.89476ud"kPa")
    io = IOBuffer()
    @test isnothing(show(io, (dimension(°F), dimension(ud"K"), psi, celsius, fahrenheit)))

    # Test updating affine units
    @test register_unit("°C" => °C) isa AbstractDict # Updating the same value does nothing for unit
    @test register_unit("K" => Dimensions(temperature=1)) isa AbstractDict # same value yields nothing for dimension
    @test_throws MethodError register_unit(ud"K") # cannot register only a unit

    # Cannot re-register a unit if its value changes nor can we delete a unit
    @test_throws RegistryTools.PermanentDictError register_unit("°C"=>ud"°F")
    @test_throws RegistryTools.PermanentDictError delete!(UnitRegistry.UNITS, :m)

    # Cannot register non-parsable unit
    @test_throws ArgumentError register_unit("m/s"=>ud"m/s")

    # Test map_dimensions
    @test map_dimensions(+, dimension(ud"m/s"), dimension(ud"m/s")) == Dimensions(length=2, time=-2)
    @test map_dimensions(-, dimension(ud"m"), dimension(ud"s"))     == Dimensions(length=1, time=-1)
    @test map_dimensions(Base.Fix1(*,2), dimension(ud"m/s"))       == Dimensions(length=2, time=-2)

    # Parsing tests 
    @test ud"1.0" == ud""
    @test q"5 kg/s" == 5ud"kg/s"
    @test uparse("1.0") == ud""
    @test qparse("10 km/hr") == 10*ud"km/hr"
    @test ud"%" == ud"percent"
    @test q"1%" == 0.01ud""
    @test q"1.0" == 1.0ud""
    @test qparse("1%") == 0.01ud""
    @test qparse("1.0") == 1.0ud""
    @test qparse("5  5kg") == 25ud"kg"

    #Test showerror 
    testio = IOBuffer()

    showerror(testio, ConversionError(ud"m/s", ud"m"))
    @test String(take!(testio)) == "ConversionError: Cannot convert unit 'm' to target unit 'm/s'. Consider multiplying 'm/s' by 's' or similar."

    showerror(testio, DimensionError(dimension(ud"m/s")))
    @test String(take!(testio)) == "DimensionError: m/s is not dimensionless"

    showerror(testio, DimensionError(1ud"m/s"))
    @test String(take!(testio)) == "DimensionError: 1.0 m/s is not dimensionless"

    showerror(testio, DimensionError((dimension(ud"m/s"), dimension(ud"s/m"))))
    @test String(take!(testio)) == "DimensionError: (m/s, s/m) have incompatible dimensions"

    showerror(testio, NotScalarError(ud"°C"))
    @test String(take!(testio)) == "NotScalarError: °C cannot be treated as scalar, operation only valid for scalar units"

    showerror(testio, NotDimensionError(ud"kPa"))
    @test String(take!(testio)) == "NotDimensionError: kPa cannot be treated as dimension, operation only valid for dimension units"
end

@testset "Static unit functionality" begin
    #Test promotion rules
    q1 = 5.0u"km/hr" |> u"km/hr"
    q2 = 20.0*u"m/s"
    q3 = 10.0u"kg/s"

    @test eltype([q1,q1]) <: Quantity{Float64, typeof(u"m/s")} #StaticUnits preserved
    @test eltype([q1,q2]) <: Quantity{Float64, typeof(dimension(u"m/s"))} #StaticUnits and StaticDims promote to StaticDims if dimension is the same
    @test eltype([q1,q2,q3]) <: Quantity{Float64, <:Dimensions} #Different dimensions promote to "Dimensions"
    @test eltype([u"m/s", u"kg/hr"]) <: Units{Dimensions{FixRat32}, AffineTransform{Float64}}
    @test eltype([u"m/s", ud"m/s"]) <: Units{Dimensions{FixRat32}, AffineTransform{Float64}}
    @test eltype([dimension(u"m/s"), dimension(u"kg/hr")]) <: Dimensions{FixRat32}
    @test eltype([dimension(u"m/s"), dimension(ud"m/s")]) <: Dimensions{FixRat32}
    @test eltype([u"m/s", u"m/s"]) <: StaticUnits{dimval(u"m/s"), AffineTransform{Float64}}
    @test eltype([dimension(u"m/s"), dimension(u"m/s")]) <: StaticDims{dimval(u"m/s")}
    

    #Addition/Sybtraction preserves static unit information but Multiplication/Division don't 
    @test (1u"m/s" + 1ud"m/s") isa Quantity{Float64, <:StaticDims}
    @test (1u"m/s" - 1ud"m/s") isa Quantity{Float64, <:StaticDims}
    @test (1u"m/s" * 1ud"m/s") isa Quantity{Float64, <:Dimensions}
    @test (1u"m/s" / 1ud"m/s") isa Quantity{Float64, <:Dimensions}
    @test (1ud"m/s" + 1u"m/s") isa Quantity{Float64, <:StaticDims}
    @test (1ud"m/s" - 1u"m/s") isa Quantity{Float64, <:StaticDims}
    @test (1ud"m/s" * 1u"m/s") isa Quantity{Float64, <:Dimensions}
    @test (1ud"m/s" / 1u"m/s") isa Quantity{Float64, <:Dimensions}


    #Ustrip tests
    x = 1.3u"km/s^2" |> u"km/s^2" 
    @test ustrip(x) == 1.3
    @test dstrip(x) == 1300  # SI base units
    @test x == q"1.3km/s^2"
    @test abs(x) === 1.3u"km/s^2"

    y = 0.9u"sqrt(mΩ)"
    @test FlexUnits.dimval(dimension(y)) == dimension(u"(m*kg^(1/2))/(s^(3/2)*A)")
    @test dstrip(y) ≈ 0.02846049894151541
    @test y ≈ q"(0.9*sqrt(mΩ))"

    y = BigFloat(0.3) * u"mΩ"
    @test typeof(y) == Quantity{BigFloat, StaticDims{FlexUnits.dimval(dimension(y))}}
    @test dstrip(y) ≈ 0.0003

    #Mixed conversions with dynamic units
    y32 = convert(Quantity{Float32, Units{Dimensions{DEFAULT_RATIONAL}, AffineTransform}}, y)
    @test typeof(y32) == Quantity{Float32, Units{Dimensions{DEFAULT_RATIONAL}, AffineTransform}}

    z = 1.0*u"yr"
    @test dstrip(z) ≈ 60 * 60 * 24 * 365.25
    @test z == 1.0*uparse("yr")
    @test z == qparse("1yr")
    @test 1/z == ubase(qparse("1/yr"))
    @test ubase(z) == ubase(qparse("1*yr"))
    @test qparse("yr") == 1*ud"yr"

    #Mathematical operations
    xp = 1u"percent"
    xv = 1u"m/s"
    @test 1u"kg/s" + 1u"kg/s" === 2u"kg/s"
    @test 1u"km/hr" - 1u"km/hr" === 0u"km/hr"
    @test 5.0u"kW"/(5.0u"kW") === 1.0u""
    @test 2.0u"m"*2.0u"m" === 4.0u"m^2"
    @test xv*1 === 1u"m/s"
    @test xv*1 !== 1u"m/s"|>u"m/s"
    @test xv*1 === ubase(1u"m/s")
    @test xp*1 == 0.01*ud""
    @test xv*xv == 1ud"m^2/s^2"
    @test *(xv,xv,xv) == 1u"m^3/s^3"
    @test cbrt(*(xv,xv,xv)) == 1u"m/s"
    @test u"m/s" / unit(1u"m/s") == u""
    @test cbrt(u"m^3/s^3") == u"m/s"
    @test abs2(dimension(u"m/s")) == dimension(u"m^2/s^2")
    @test (1.0u"m/s")' == 1.0u"m/s"

    @test xv/1 == 1u"m/s"
    @test xp/1 == 0.01*u""
    @test xv/xv == 1u""

    @test sqrt(4*xv) == 2u"m^0.5/s^0.5"
    @test (4*xv)^0.5 == 2u"m^0.5/s^0.5"
    @test (2*xv)^2 == 4u"m^2/s^2"

    @test (1.0u"m")^0 === 1.0u""
    @test (1.0u"m")^1 === 1.0u"m"
    @test (1.0u"m")^2 === 1.0u"m^2"
    @test (1.0u"m")^3 === 1.0u"m^3"
    @test (1.0u"m")^4 === 1.0u"m^4"
    @test (1.0u"m")^(-1) === 1.0u"1/m"
    @test (1.0u"m")^(-2) === 1.0u"m^(-2)"
    @test (1.0u"m")^(-3) === 1.0u"1/m^3"
    

end

@testset "Type conversions" begin
    d = Dimensions{Rational{Int16}}(mass=2)
    d32 = convert(Dimensions{Rational{Int32}}, d)
    @test typeof(d) == Dimensions{Rational{Int16}}
    @test typeof(d32) == Dimensions{Rational{Int32}}
    @test convert(Quantity{Float64, typeof(ud"m/s")}, 1.0) isa Quantity{Float64, Units{Dimensions{FixRat32}, AffineTransform{Float64}}}
    @test convert(Quantity{Float64, FlexUnits.dimtype(ud"m/s")}, 1.0) isa Quantity{Float64, Dimensions{FixRat32}}

    # Should not change:
    @test convert(Dimensions{Rational{Int16}}, d) === d

    #Quantity type conversions
    q = Quantity(0.5, d)
    q32_32 = convert(Quantity{Float32,Dimensions{Rational{Int32}}}, q)
    @test typeof(q) == Quantity{Float64,Dimensions{Rational{Int16}}}
    @test typeof(q32_32) == Quantity{Float32,Dimensions{Rational{Int32}}}
    @test ustrip(q) == 0.5
    @test ustrip(q32_32) == 0.5
    @test typeof(ustrip(q)) == Float64
    @test typeof(ustrip(q32_32)) == Float32
    @test dimension(q32_32) == convert(typeof(dimension(q32_32)), dimension(q))
    @test convert(Quantity, q) === q
    @test convert(Quantity{Float64, Units{DEFAULT_DIM_TYPE, AT}}, ubase(1ud"kg")) isa Quantity{Float64, Units{DEFAULT_DIM_TYPE, AT}}
    
    #Static unit type conversions 
    @test convert(Quantity{Int64, typeof(u"m/s")}, 1u"m/s") isa Quantity{Int64, typeof(u"m/s")}
    @test convert(Quantity{Int64, typeof(dimension(u"m/s"))}, 1u"m/s") isa Quantity{Int64, typeof(dimension(u"m/s"))}
    @test_throws ConversionError convert(Quantity{Float64, typeof(u"m/s")}, 1u"kg/hr")
    @test_throws ConversionError convert(Quantity{Int64, typeof(dimension(u"m/s"))}, 1u"kg/hr")

    # Test that regular type promotion applies:
    q = Quantity(2, d)
    @test typeof(q) == Quantity{Int64,typeof(d)}
    @test typeof(q ^ 2) == Quantity{Int64,typeof(d)}
    @test typeof(0.5 * q) == Quantity{Float64,typeof(d)}
    @test typeof(inv(q)) == Quantity{Float64,typeof(d)}

    # Test conversion of unit types 
    @test convert(DEFAULT_DIM_TYPE, ud"m") === Dimensions(length=1)
    @test_throws NotDimensionError convert(DEFAULT_DIM_TYPE, ud"mm")
    @test convert(Units{DEFAULT_DIM_TYPE, AT}, ubase(2ud"m")) == Units(todims=AffineTransform(scale=2.0, offset=0.0), dims=dimension(ud"m"))
    @test_throws ArgumentError convert(Units{DEFAULT_DIM_TYPE, AT}, quantity(2, ud"°C"))
    @test convert(Quantity{Float64, DEFAULT_DIM_TYPE}, 2ud"m") === Quantity{Float64, DEFAULT_DIM_TYPE}(2.0, dimension(ud"m")) 
    @test_throws NotScalarError convert(Quantity{Float64, DEFAULT_DIM_TYPE}, ud"°C") 
    @test promote_type(Quantity{Float32, DEFAULT_DIM_TYPE}, Quantity{Float64, DEFAULT_UNIT_TYPE}) == Quantity{Float64, DEFAULT_DIM_TYPE}
    @test promote_type(Quantity{Float64, typeof(ud"m/s")}, Float64) == Quantity{Float64, typeof(dimension(ud"m/s"))}
    @test promote_type(Quantity{Float32, typeof(ud"m/s")}, Float64) == Quantity{Float64, typeof(dimension(ud"m/s"))}
    @test promote_type(FlexQuant{Matrix{Float32}, typeof(ud"m/s")}, Matrix{Float64}) == FlexQuant{Matrix{Float64}, typeof(dimension(ud"m/s"))}

    # Test that adding different dimension subtypes still works
    @test 1*Dimensions{Int64}(length=1) + 1ud"m" == 2ud"m"

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
    @test uconvert(ud"nm", 5e-9ud"m") ≈ (5e-9ud"m" |> ud"nm") ≈ 5ud"nm"
    @test_throws ConversionError uconvert(ud"nm*J", 5e-9ud"m")

    #Generating transforms 
    @test uconvert(u"kg/hr", u"kg/s") == AffineTransform(3600.0, 0.0)
    @test uconvert(u"kg/s" |> u"kg/hr", 1.0) == 3600.0
    @test uconvert(dimension(u"kg/s"), u"kg/hr") == FlexUnits.todims(u"kg/hr")

    #dconvert and similar 
    @test dconvert(u"kg/s", 5.0u"kg/hr") == ubase(5.0u"kg/hr")
    @test ustrip(u"kg/s", 5.0ud"kg/hr") ≈ 5/3600
    @test ustrip([u"kg/hr", u"m^3/hr"], [1.0u"kg/s",2.0u"L/s"]) ≈ [3600, 7.2]
    @test FlexUnits.ustrip_dimensionless(1.0u"") == 1.0
    @test_throws DimensionError FlexUnits.ustrip_dimensionless(1.0u"m/s")

    # Types:
    @test typeof(uconvert(ud"nm", 5e-9ud"m")) <: Quantity{Float64, U} 
    @test typeof(uconvert(ud"nm", Quantity(5e-9, ud"m"))) <: Quantity{Float64, U}
    @test uconvert(ud"nm", Quantity(5e-9, ud"m")) ≈ 5ud"nm"

    @test dimension(1ud"m" |> ud"nm")[:length] == 1

   
    # Different types require converting both arguments:
    q = Quantity{Float16}(1.5ud"g")
    @test q isa Quantity{Float16}

    # Broadcasting conversions over Arrays
    x = [1.0, 2.0, 3.0] .* Quantity(1, ud"kg")
    x2 = x .|> ud"g"
    @test typeof(x2) <: Vector{<:Quantity{Float64,<:Units{<:Any}}}
    @test x2[2] ≈ (2000ud"g")
end

@testset "Stats and Linear Algebra" begin
    using StaticArrays
    import Random
    using Statistics

    #Nonlinear map
    @kwdef struct PumpInput{T} <: FieldVector{2,T}
        current :: T 
        voltage :: T
    end

    @kwdef struct PumpOutput{T} <: FieldVector{3,T}
        power :: T 
        pressure :: T
        flow :: T 
    end

    function pumpfunc(x::PumpInput)
        p = x.current*x.voltage*0.9   
        return PumpOutput(power = p, pressure = sqrt(p), flow = sqrt(p))
    end
    pumpfunc(x::AbstractVector) = pumpfunc(PumpInput(x))

    Random.seed!(1234)

    #Generate a correlated test set
    U = [ud"kg/s", ud"kW", ud"Hz"]
    X = (rand(30,3)*rand(3,3) .+ 0.001.*randn(30,3))
    Q = ubase.(X .* U')

    #Test summations, stats and some linear algebra
    @test all(sum(Q, dims=1) .≈ sum(X, dims=1).*U')
    @test all(mean(Q, dims=1) .≈ mean(X, dims=1).*U')
    @test all(var(Q, dims=1) .≈ var(X, dims=1).*(U.^2)')
    @test all(cov(Q) .≈ cov(X).*U.*U')
    @test all((cor(Q) .+ 0u"") .≈ cor(X))
    @test sum(Q*inv.(U)) ≈ sum(X)
    @test all(minimum(Q, dims=1, init=typemax(eltype(Q))) .≈ minimum(X, dims=1).*U')
    @test all(maximum(Q, dims=1, init=typemin(eltype(Q))) .≈ maximum(X, dims=1).*U')


    #Quick linear algebra tests 
    u1 = SA[u"lbf*ft", u"kW", u"rpm"]
    u2 = SA[u"kg/s", u"m^3/hr", u"kW"]

    xm = SMatrix{3,3}(randn(3,3))
    qMraw = xm.*u2./u1'
    qM = LinmapQuant(qMraw)

    #Matrix inversion
    x = SVector{3}(randn(3)).*u1
    y = qM*x
    @test all(x .≈ inv(qM)*y)

    #Matrix transpose
    @test all(Matrix(y') .≈ Matrix(x'*qM'))

    #Square matrices
    Σ = cov(randn(20,3)*rand(3,3))
    x = randn(3).*u2

    #Symmetric matrix
    rS = Σ.*inv.(u2).*inv.(u2)'
    qS = LinmapQuant(rS)
    @test all(qS .≈ ubase.(rS))
    @test x'*(rS)*x ≈ x'*qS*x

    #Repeatable matrix
    rR = Σ.*u2.*inv.(u2)'
    qR = LinmapQuant(rR)
    @test all(qR .≈ ubase.(rR))
    @test all((rR^2*x) .≈ (qR*qR*x))

    #Nonlinear mapping
    pumpunits = UnitMap(PumpInput(current=u"A", voltage=u"V"), PumpOutput(power=u"W", pressure=u"Pa", flow=u"m^3/s"))
    upumpfunc = FunctionQuant(pumpfunc, pumpunits)
    qinput = PumpInput(current=500*u"mA", voltage=6u"V")
    @test all(upumpfunc(qinput) .≈ pumpfunc(ustrip.(uinput(pumpunits), qinput)).*uoutput(pumpunits))


    #Matrix and vector operations 
    m  = LinmapQuant(SA[1.0 0.1; 0.2 1.0], UnitMap(u_in = SA[u"kg/s", u"kW"], u_out=SA[u"m^3/s", u"kPa"]))
    mi = inv(m)
    x  = VectorQuant(SA[0.5u"kg/s", 0.5u"kW"])
    y  = m*x

    mraw  = SMatrix{2,2}(m)
    miraw = SMatrix{2,2}(mi)
    xraw  = SVector{2}(x)
    yraw  = SVector{2}(y)
     
    @test x + x ≈ xraw + xraw
    @test m + m ≈ mraw + mraw
    @test x - x ≈ xraw - xraw
    @test m - m ≈ mraw - mraw

    @test m*mi ≈ mraw*miraw
    @test m\m  ≈ miraw*mraw
    @test m\mraw ≈ miraw*mraw
    @test m'\m' ≈ miraw'*mraw'
    @test m'/m' ≈ mraw'*miraw'
    @test m/m  ≈ mraw*miraw
    @test mraw/m ≈ mraw*miraw
    @test m*x ≈ yraw
    @test mi*y ≈ xraw
    @test (x'*m')' ≈ yraw 
    @test (y'*mi')' ≈ xraw
    @test m\y ≈ xraw
    @test (y'/m')' ≈ xraw
    @test transpose(transpose(y)) ≈ yraw

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
    @test FixedRational{1000}(10)*true == FixedRational{1000}(10)

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
    @test convert(Rational, FixedRational{6,Int8}(2)) === Rational{Int8}(2)

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
    x = 10ud"m"
    user_quantity = Quantity(10.0, Dimensions{FixedRational{25200,Int32}}(1, 0, 0, 0, 0, 0, 0))
    @test x == user_quantity
end

#Register a new affine unit (and verify re-registering)
register_unit("psig" => Units(todims=AffineTransform(scale=uscale(ud"psi"), offset=101.3ud"kPa"), dims=ud"Pa"))

@testset "Registration tests" begin
    #Test re-registering and verify that the unit exists
    register_unit("psig" => Units(todims=AffineTransform(scale=uscale(ud"psi"), offset=101.3ud"kPa"), dims=ud"Pa"))
    @test 0*ud"psig" == 101.3ud"kPa"
    @test q"0psig" == 101.3ud"kPa"
    
    #Test registration for a different registry base type
    IntDimType = Dimensions{Int32}
    reg = RegistryTools.PermanentDict{Symbol, Units{IntDimType, AffineTransform{Float64}}}()
    reg = RegistryTools.registry_defaults!(reg)  
    @test reg[:m]  === Units(todims=AffineTransform(), dims=IntDimType(length=1), symbol=:m)
    @test reg[:kg] === Units(todims=AffineTransform(), dims=IntDimType(mass=1), symbol=:kg)    
end

@testset "Integration tests with Unitful" begin
    q1 = 5.0*Unitful.u"km/hr"
    q2 = 10.0*ud"°C"
    q3 = 15*Unitful.u"cd/mol"

    @test uconvert(Unitful.unit(q1), uconvert(ud"m/s", q1)) == q1
    @test Quantity(q1) == 5.0*ud"km/hr"
    @test Quantity{Float64}(q1) == 5.0*ud"km/hr"
    @test convert(Quantity{Float64, UnitRegistry.unittype()}, q1) == 5.0*ud"km/hr"
    @test_throws DimensionError uconvert(ud"kPa", q1)
    @test Unitful.ustrip(uconvert(Unitful.u"K", q2)) == ustrip(uconvert(ud"K", q2))
    @test Unitful.ustrip(Unitful.uconvert(Unitful.u"cd/mol", q3)) == ustrip(uconvert(ud"cd/mol", q3))
end

@testset "Aqua.jl" begin
    Aqua.test_all(FlexUnits)
end

#=
#Benchmark testing
R = Ref(8.314 * u"J/(mol*K)")
@benchmark let R=R[]
    v_satp = R * (25u"°C") / (101.3u"kPa")
end
=#

#Profiling diagnostic if benchmarks are off 
#=
v = randn(100000).*u""
v2 = deepcopy(v)
@profview for ii in 1:1000
    v2 .= v.*(1.1*u"m/s")
end
=#