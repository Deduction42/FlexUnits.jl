using Revise
using FlexUnits, .UnitRegistry
using OrdinaryDiffEq
using StaticArrays
using Plots
using BenchmarkTools

#Make "getindex" return a unitless version of the value
abstract type QuantFieldVector{N,T} <: FieldVector{N,T} end
Base.@propagate_inbounds Base.getindex(a::QuantFieldVector, i::Int) = dstrip(getfield(a, i))

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

function acceleration_raw(u0::AbstractVector, p::FallingObjectProps, t)
    u = FallingObjectState(u0)

    fd = -sign(u.v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    dv = fd/p.m - p.g
    dh = u.v

    #Need to multiply by seconds to match state units (the ODE solver will be unitless)
    return FallingObjectState(v=dv*(1u"s"), h=dh*(1u"s"))
end


# =============================================================================================================
println("\nExplicit Numerical Solution")
# =============================================================================================================
u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"g/L", m=50u"kg", g=9.81u"m/s^2")

tspan = (0.0, 10.0)
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_raw, u0, tspan, p, abstol=[1e-6, 1e-6], reltol=[1e-6, 1e-6])
sol = solve(prob, Tsit5())
@btime solve(prob, Tsit5())
plt = plot(sol.t, [u[1] for u in sol.u], label="explicit")



# =============================================================================================================
println("\nImplicit Numerical Solution")
# =============================================================================================================
u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")

tspan = (0.0, 10.0)
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_raw, u0, tspan, p, abstol=[1e-6, 1e-6], reltol=[1e-6, 1e-6])
sol = solve(prob, Rodas5P())
@btime solve(prob, Rodas5P())
plt = plot!(plt, sol.t, [u[1] for u in sol.u], label="implicit")