using Revise
using FlexUnits, .UnitRegistry
using OrdinaryDiffEq
using StaticArrays
using Plots
using BenchmarkTools
using LinearAlgebra

# Required additional methods from DiffEqBase's Unitful Extension =============================================
#OrdinaryDiffEq.OrdinaryDiffEqCore.DiffEqBase.UNITLESS_ABS2(q::Quantity) = abs2(dstrip(q))
#OrdinaryDiffEq.OrdinaryDiffEqCore.DiffEqBase.value(q::Quantity) = dstrip(q)

# =============================================================================================================
include("diffeq_extension.jl")

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

function acceleration_ustatic(u::AbstractVector{<:Quantity}, p::AbstractVector{<:Quantity}, t)
    du = acceleration_raw(ustatic(FallingObjectState(u)), ustatic(FallingObjectProps(p)), t)
    return FallingObjectState(du)
end

#Jacobian needs to be wrapped around a "Real" version to support autodiff
using ForwardDiff
p  = FallingObjectProps(Cd=1.0u"", A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")

function wrapped_jacobian(f, x, p, t)
    u_in  = dimension.(x)
    u_out = dimension.(f(x,p,t))

    function f_unitless(xn)
        return dstrip.(f(xn.*u_in, p, t))
    end
    
    return LinmapQuant(ForwardDiff.jacobian(f_unitless, dstrip.(x)), DimsMap(u_in=u_in, u_out=u_out))
end

wrapped_jacobian(acceleration_ustatic, FallingObjectState(0.0u"m/s", 100.0u"m"), p, 0*u"s")


# =============================================================================================================
println("\nStatic-Unit Solution")
# =============================================================================================================
u0 = FallingObjectState(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps(Cd=1.0u"", A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")
static_jac(u,p,t) = wrapped_jacobian(acceleration_ustatic, u, p, t)
static_tgrad(u,p,t) = [0.0u"m/s^3", 0.0u"m/s^2"]

tspan = (0.0u"s", 10.0u"s")
f_static = ODEFunction{false, OrdinaryDiffEq.SciMLBase.FullSpecialize}(acceleration_ustatic, 
    jac = static_jac, 
    tgrad = static_tgrad,
    mass_matrix = I*1ud""
)
prob = ODEProblem(f_static, u0, tspan, p, abstol=[1e-6, 1e-6], reltol=[1e-6, 1e-6])


sol = solve(prob, Rodas5P())
plt = plot()
plt = plot!(plt, ustrip.(sol.t), [ustrip(u.v) for u in sol.u], label="v_staticu")


#= This triggers an error at 
rosenprock_perform_step.jl line 1240
dT = calc_tderivative(integrator, cache) = Quantity{Float64, Dimensions{FixRat32}}[0.0 m/s², 0.0 m/s]
du = f(uprev, p, t) = Quantity{Float64, Dimensions{FixRat32}}[-9.81 m/s², 0.0 m/s]
dtd = Quantity{Float64, StaticDims{s}}[0.016130215329083333 s, -0.03226043065816667 s, -0.025759833949205672 s, 0.13734855038363272 s, 0.1770147199103226 s, 0.0 s, 0.0 s, 0.0 s]
dtd[1]*dT = Quantity{Float64, Dimensions{FixRat32}}[0.0 m/s, 0.0 m]

du + dtd[1]*dT = [-9.81 m/s², 0.0 m/s] + [0.0 m/s, 0.0 m] => Inconsistent units
=#