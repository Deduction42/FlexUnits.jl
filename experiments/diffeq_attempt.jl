using OrdinaryDiffEq
using FlexUnits, .UnitRegistry
using StaticArrays

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
    fd = -sign(v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    dv = fd - p.g*p.m
    dh = u.v
    return FallingObjectState(v=dv, h=dh)
end

function acceleration_float(u::AbstractVector{<:Real}, p::FallingObjectProps{<:Real}, t) 
    return acceleration_raw(FallingObjectState(ustrip.(u)), FallingObjectProps(ustrip.(p)), t)
end

function acceleration_static(u::AbstractVector{<:Quantity}, p::AbstractVector{<:Quantity}, t)
    return acceleration_raw(ustatic(FallingObjectState(u)), ustatic(FallingObjectProps(ustrip.(p))), t)
end

function acceleration_dynamic(u::AbstractVector{<:Quantity}, p::AbstractVector{<:Quantity}, t)
    return acceleration_raw(FallingObjectState(u)), FallingObjectProps(ustrip.(p), t)
end

# Required additional methods =================================================================================
FlexUnits.Quantity{T,U}(x) where {T,U<:StaticDims} = Quantity{T,U}(convert(T, x), U())
FlexUnits.Quantity{T,U}(x::QuantUnion) where {T,U<:StaticDims} = Quantity{T,U}(convert(T, ustrip(U(),x)), U())
Base.oneunit(::Type{Quantity{T,U}}) where {T,U<:StaticDims} = Quantity{T,U}(one(T), U())
Base.eltype(::Type{Quantity{T}}) where T = T
OrdinaryDiffEq.OrdinaryDiffEqCore.DiffEqBase.UNITLESS_ABS2(q::Quantity) = abs2(dstrip(q))
# =============================================================================================================


u0 = FallingObjectState(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps(Cd=1.0u"", A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")
tspan = (0.0u"s", 100.0u"s")
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_dynamic, u0, tspan, p, abstol=[1e-6u"m/s", 1e-6u"m"], reltol=[1e-6, 1e-6])

# Test that it worked
sol = solve(prob, Tsit5())
using Plots; plot(sol, vars=(1,2))