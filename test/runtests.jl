#============================================================================================================================
Run these commands at startup to see coverage
julia --startup-file=no --depwarn=yes --threads=auto -e 'using Coverage; clean_folder("src"); clean_folder("test"); clean_folder("ext") '
julia --startup-file=no --depwarn=yes --threads=auto --code-coverage=user --project=. -e 'using Pkg; Pkg.test(coverage=true)'
julia --startup-file=no --depwarn=yes --threads=auto coverage.jl

Run this command for testing invalidations
julia --startup-file=no --depwarn=yes --threads=auto --project=. test/invalidations.jl
============================================================================================================================#

#To see the actual coverage in VSCode, install the Coverage Gutters extension
#https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters

using TestItems: @testitem
using TestItemRunner

using Revise

@testitem "Basic Functionality" begin
    include("tests_basic.jl")
end

@testitem "Linear Algebra" begin 
    include("test_linear_algebra.jl")
end

@testitem "Statistics" begin 
    include("test_statistics.jl")
end

@testitem "Integration tests with Unitful" begin
    import Unitful 
    using FlexUnits, .UnitRegistry 

    q1 = 5.0*Unitful.u"km/hr"
    q2 = 10.0*ud"°C"
    q3 = 15*Unitful.u"cd/mol"

    @test uconvert(Unitful.unit(q1), uconvert(ud"m/s", q1)) == q1
    @test Quantity(q1) == 5.0*ud"km/hr"
    @test Quantity{Float64}(q1) == 5.0*ud"km/hr"
    @test convert(Quantity{Float64, UnitRegistry.utype()}, q1) == 5.0*ud"km/hr"
    @test_throws DimensionError uconvert(ud"kPa", q1)
    @test Unitful.ustrip(uconvert(Unitful.u"K", q2)) == ustrip(uconvert(ud"K", q2))
    @test Unitful.ustrip(Unitful.uconvert(Unitful.u"cd/mol", q3)) == ustrip(uconvert(ud"cd/mol", q3))
end

@testitem "Aqua.jl" begin
    using Aqua
    Aqua.test_all(FlexUnits)
end

nothing

