using Revise
using FlexUnits, .UnitRegistry
using OrdinaryDiffEq
using StaticArrays
using Plots
using BenchmarkTools

@kwdef struct FallingObjectState{T} <: FieldVector{2,T}
    v  :: T
    h  :: T
end

@kwdef struct FallingObjectProps{T} <: FieldVector{5,T}
    Cd :: T
    A  :: T
    ρ  :: T
    m  :: T
    g  :: T
end 

const STATE_UNITS = FallingObjectState(v=ud"m/s", h=ud"m")

function ustatic(state::FallingObjectState{<:Quantity})
    return (
        v = dconvert(u"m/s", state.v),
        h = dconvert(u"m", state.h)
    )
end

function ustatic(props::FallingObjectProps{<:Quantity})
    return (
        Cd = dconvert(u"", props.Cd),
        A  = dconvert(u"m^2", props.A),
        ρ  = dconvert(u"kg/m^3", props.ρ),
        m  = dconvert(u"kg", props.m),
        g  = dconvert(u"m/s^2", props.g)
    )
end

function acceleration_raw(u, p, t)
    fd = -sign(u.v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    dv = fd/p.m - p.g
    dh = u.v
    return FallingObjectState(v=dv, h=dh)
end

function acceleration_float(u::AbstractVector{<:Real}, p::FallingObjectProps{<:Real}, t) 
    return acceleration_raw(FallingObjectState(ustrip.(u)), FallingObjectProps(ustrip.(p)), t)
end

function acceleration_static(u::AbstractVector{<:Quantity}, p::AbstractVector{<:Quantity}, t)
    du = acceleration_raw(ustatic(FallingObjectState(u)), ustatic(FallingObjectProps(p)), t)
    return FallingObjectState(du)
end

function acceleration_dynamic(u::AbstractVector{<:Quantity}, p::AbstractVector{<:Quantity}, t)
    du = acceleration_raw(FallingObjectState(u), FallingObjectProps(p), t)
    return FallingObjectState(du)
end

# Required additional methods from DiffEqBase's Unitful Extension =============================================
OrdinaryDiffEq.OrdinaryDiffEqCore.DiffEqBase.UNITLESS_ABS2(q::Quantity) = abs2(dstrip(q))
OrdinaryDiffEq.OrdinaryDiffEqCore.DiffEqBase.value(q::Quantity) = dstrip(q)
# =============================================================================================================


# =============================================================================================================
println("\nRaw Numerical Solution")
# =============================================================================================================
u0 = FallingObjectState(v=0.0, h=100)
p  = FallingObjectProps(Cd=1.0, A=0.1, ρ=1.0, m=50, g=9.81)

tspan = (0.0, 10.0)
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_raw, u0, tspan, p, abstol=[1e-6, 1e-6], reltol=[1e-6, 1e-6])
sol = solve(prob, Tsit5())
@btime solve(prob, Tsit5())
plt = plot(sol.t, [u.v for u in sol.u], label="v_unitless")

# =============================================================================================================
println("\nStatic Unit Solution")
# =============================================================================================================
u0 = FallingObjectState(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps(Cd=1.0u"", A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")

tspan = (0.0u"s", 10.0u"s")
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_static, u0, tspan, p, abstol=[1e-6, 1e-6], reltol=[1e-6, 1e-6])

sol = solve(prob, Tsit5())
@btime solve(prob, Tsit5())
plt = plot!(plt, ustrip.(sol.t), [ustrip(u.v) for u in sol.u], label="v_static")

# =============================================================================================================
println("\nDynamic Unit Solution")
# =============================================================================================================
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_dynamic, u0, tspan, p, abstol=[1e-6, 1e-6], reltol=[1e-6, 1e-6])

sol = solve(prob, Tsit5())
@btime solve(prob, Tsit5())
plt = plot!(plt, ustrip.(sol.t), [ustrip(u.v) for u in sol.u], label="v_dynamic")