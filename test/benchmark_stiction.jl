using Revise
using FlexUnits, .UnitRegistry
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using StaticArrays 
using TimeRecords

import OrdinaryDiffEqTsit5.SciMLBase.FullSpecialize

#==============================================================================================
This is a demo of the LuGre stiction model
https://www.mathworks.com/matlabcentral/fileexchange/84792-lugre-friction-model-in-matlab-gnu-octave/files/demo2.m
https://github.com/auralius/LuGre/blob/master/README.md
==============================================================================================#
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
    α = 2 # used to approximate Stribeck effect, can range from 0.5 <= α <= 2.0
   
    gv = θ.Fc + θ.Fs*exp(-abs(x.v/θ.vs)^α)
    ż  = x.v - (θ.σ₀/gv)*x.z*abs(x.v)
    Ff = θ.σ₀*x.z + θ.σ₁*ż + θ.μD*x.v
    Fnet = (θ.k*(u-x.x) - Ff) #Force balance

    return dimsmod(D"1/s", ValveState,
        x = x.v,
        v = Fnet/θ.m, 
        z = ż
    )
end


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

x0 = ValveState{Float64}(x=0, v=0, z=0)
abstol = ValveState{Float64}(x=1e-6u"m", v=1e-6u"m/s", z=1e-6u"m")
reltol = SA[1e-6, 1e-6, 1e-6]

prob = ODEProblem{false, FullSpecialize}(lugre_diff, x0, Δt, θ, abstol=abstol, reltol=reltol)
@time sol = solve(prob, Rodas4P())

using Plots
vt = LinRange(Δt[begin], Δt[end], 100)
ax = plot(vt, dstrip.(interpolate(θ.u, vt, order=1)), label="setpoint")
plot!(ax, sol.t, [dstrip(x.x) for x in sol.u], label="actual")




