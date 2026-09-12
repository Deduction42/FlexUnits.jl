using Revise
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using StaticArrays 
using FlexUnits, .UnitRegistry
using Test
using TimeRecords

import OrdinaryDiffEqTsit5.SciMLBase.FullSpecialize
import FlexUnits: QuantFieldVector, DimsMod

#============================================================================================
Falling Object Model
============================================================================================#

@kwdef struct FallingObjectState{T} <: QuantFieldVector{2,T}
    v  :: Quantity{T, D"m/s"}
    h  :: Quantity{T, D"m"}
end

@kwdef struct FallingObjectProps{T} <: QuantFieldVector{5,T}
    Cd :: Quantity{T, D""}
    A  :: Quantity{T, D"m^2"}
    ρ  :: Quantity{T, D"kg/m^3"}
    m  :: Quantity{T, D"kg"}
    g  :: Quantity{T, D"m/s^2"}
end 

function acceleration(u0::AbstractVector, p::FallingObjectProps, t)
    u = FallingObjectState(u0)

    #Drag force
    fd = -sign(u.v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    
    #Drag force effect on state (multiply by dt to make units work)
    dv = (fd/p.m - p.g)
    dh = u.v

    return convert(typeof(u0), DimsMod{D"1/s"}(FallingObjectState, (v=dv, h=dh)))
end


@testset "Falling Object Simulation (Explicit Solver)" begin
    u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
    p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")
    abstol = FallingObjectState{Float64}(v=1e-6u"m/s", h=1e-6u"m")
    reltol = SA[1e-6, 1e-6]

    tspan = dstrip.((0.0u"min", 0.25u"min")) #Time span must be in seconds, dstrip takes care of this
    prob = ODEProblem{false, FullSpecialize}(acceleration, u0, tspan, p, abstol=abstol, reltol=reltol)
    sol = solve(prob, Tsit5())

    #=
    #Benchmarking and plotting
    using BenchmarkTools 
    #@btime soln = solve(prob, Tsit5())

    using Plots
    #plt = plot!(plt, sol.t, [dstrip(u.v) for u in sol.u], label="explicit units") #Each element in sol.u is a QuantFieldVector
    =#

    @test sol.retcode == ReturnCode.Success
end


@testset "Falling Object Simulation (Implicit Solver)" begin
    u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
    p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")
    abstol = FallingObjectState{Float64}(v=1e-6u"m/s", h=1e-6u"m")
    reltol = SA[1e-6, 1e-6]

    tspan = dstrip.((0.0u"min", 0.25u"min"))
    prob = ODEProblem{false, FullSpecialize}(acceleration, u0, tspan, p, abstol=abstol, reltol=reltol)
    sol = solve(prob, Rodas5P())

    #=
    #Benchmarking and plotting
    using BenchmarkTools
    @btime soln = solve(prob, Rodas5P())

    using Plots
    plt = plot!(plt, sol.t, [dstrip(u.v) for u in sol.u], label="implicit units") #Each element in sol.u is a QuantFieldVector
    =#

    @test sol.retcode == ReturnCode.Success
end

