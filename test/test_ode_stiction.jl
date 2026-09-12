
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
Valve Model
============================================================================================#
@kwdef struct ValveState{T} <: QuantFieldVector{3,T}
    x :: Quantity{T, D"m"} #Position
    v :: Quantity{T, D"m/s"} #Velocity 
    z :: Quantity{T, D"m"} #Internal state position
end

@kwdef struct ValveInputs{T}
    u  :: RegularTimeSeries{Quantity{T, D"m"}} #Position setpoint
    k  :: Quantity{T, D"N/m"} #Spring constant
    m  :: Quantity{T, D"kg"} #System mass
    Fc :: Quantity{T, D"N"} #Coulomb friction 
    Fs :: Quantity{T, D"N"} #Static ADDITIVE Friction: Total static friction = Fs + Fc
    μD :: Quantity{T, D"Pa*s*m"} #Viscosity * characteristic_length
    σ₀ :: Quantity{T, D"N/m"} #Stiffness
    σ₁ :: Quantity{T, D"Pa*s*m"} #Micro damping (viscosity * characteristic_length)
    vs :: Quantity{T, D"m/s"} #Velocity below which static friction starts to dominate
end

function lugre_diff(xvec::AbstractVector, θ, t)
    x = ValveState(xvec)
    u = interpolate(θ.u, t, order=1)
   
    gv = θ.Fc + θ.Fs*exp(-abs(x.v/θ.vs))
    ż  = x.v - (θ.σ₀/gv)*x.z*abs(x.v)
    Ff = θ.σ₀*x.z + θ.σ₁*ż + θ.μD*x.v
    Fnet = (θ.k*(u-x.x) - Ff) #Force balance

    #Validate units of the results, then return the unitless version
    dx = DimsMod{D"1/s"}(ValveState, (x = x.v, v = Fnet/θ.m, z = ż))
    return SVector(dx)
end

@testset "Valve Stiction Simulation (Implicit Solver)" begin     
    Δt = (0.0, 30.0)
    vt = LinRange(Δt[begin], Δt[end], 10)
    θ  = ValveInputs{Float64}(
        u  = RegularTimeSeries(vt, 0.01.*vt.*u"m"),
        k  = 40u"N/m",
        m  = 1u"kg",
        σ₀ = 1e5u"N/m",
        σ₁ = 3e2u"kg/s",
        μD = 0.4u"kg/s",
        Fc = 1.0u"N",
        Fs = 0.5u"N",
        vs = 0.001u"m/s"
    )

    x0 = SVector(ValveState{Float64}(x=0, v=0, z=0))
    abstol = SVector(ValveState{Float64}(x=1e-6u"m", v=1e-6u"m/s", z=1e-6u"m"))
    reltol = SA[1e-6, 1e-6, 1e-6]

    prob = ODEProblem{false, FullSpecialize}(lugre_diff, x0, Δt, θ, abstol=abstol, reltol=reltol)
    sol = solve(prob, Rodas4P())

    @test sol.retcode == ReturnCode.Success

    #=
    #Benchmarking and plotting
    using BenchmarkTools
    @btime solve(prob, Rodas4P())

    using Plots
    vt = LinRange(Δt[begin], Δt[end], 100)
    ax = plot(vt, dstrip.(interpolate(θ.u, vt, order=1)), label="setpoint")
    plot!(ax, sol.t, [dstrip(x.x) for x in sol.u], label="actual")
    =#

end

