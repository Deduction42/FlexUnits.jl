using Revise
using FlexUnits, .UnitRegistry
using OrdinaryDiffEq
using StaticArrays
using Plots
using BenchmarkTools


# =============================================================================================================
@info "Raw Numerical Solutions"
# =============================================================================================================
@kwdef struct NumFallingState{T} <: FieldVector{2,T}
    v  :: T
    h  :: T
end

@kwdef struct NumFallingProps{T} <: FieldVector{5,T}
    Cd :: T
    A  :: T
    ρ  :: T
    m  :: T
    g  :: T
end 

function acceleration(u0::AbstractVector, p::NumFallingProps, t)
    u = NumFallingState(u0)

    #Drag force
    fd = -sign(u.v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    
    #Drag force effect on state
    dv = (fd/p.m - p.g)
    dh = u.v

    return NumFallingState(v=dv, h=dh)
end

# =============================================================================================================
@info "Explicit Solution Without Units"
# =============================================================================================================
u0 = NumFallingState{Float64}(v=0.0, h=100)
p  = NumFallingProps{Float64}(Cd=1.0, A=0.1, ρ=1.0, m=50, g=9.81)
abstol = NumFallingState{Float64}(v=1e-6, h=1e-6)
reltol = SA[1e-6, 1e-6]

tspan = (0.0, 15) #Time span must be in seconds
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.FullSpecialize}(acceleration, u0, tspan, p, abstol=abstol, reltol=reltol)
sol = solve(prob, Tsit5())
@btime soln = solve(prob, Tsit5())
plt = plot(sol.t, [u.v for u in sol.u], label="explicit no units") 

# =============================================================================================================
@info "Implicit Solution Without Units"
# =============================================================================================================
u0 = NumFallingState{Float64}(v=0.0, h=100)
p  = NumFallingProps{Float64}(Cd=1.0, A=0.1, ρ=1.0, m=50, g=9.81)

tspan = (0.0, 15)
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.FullSpecialize}(acceleration, u0, tspan, p, abstol=abstol, reltol=reltol)
sol = solve(prob, Rodas5P())
@btime soln = solve(prob, Rodas5P())
plt = plot!(plt, sol.t, [u.v for u in sol.u], label="implicit no units") #Each element in sol.u is a QuantFieldVector




# =============================================================================================================
@info "FlexUnits Solutions"
# =============================================================================================================

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

function acceleration(u0::AbstractVector, p::FallingObjectProps, t)
    u = FallingObjectState(u0)
    dt = D"s"() #Time dimension

    #Drag force
    fd = -sign(u.v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    
    #Drag force effect on state (multiply by dt to make units work)
    dv = (fd/p.m - p.g)*dt
    dh = u.v*dt

    return FallingObjectState(v=dv, h=dh)
end

# =============================================================================================================
@info "Explicit Solution With Units"
# =============================================================================================================
u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")
abstol = FallingObjectState{Float64}(v=1e-6u"m/s", h=1e-6u"m")
reltol = SA[1e-6, 1e-6]

tspan = dstrip.((0.0u"min", 0.25u"min")) #Time span must be in seconds, dstrip takes care of this
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.FullSpecialize}(acceleration, u0, tspan, p, abstol=abstol, reltol=reltol)
sol = solve(prob, Tsit5())
@btime soln = solve(prob, Tsit5())
plt = plot!(plt, sol.t, [dstrip(u.v) for u in sol.u], label="explicit units") #Each element in sol.u is a QuantFieldVector

# =============================================================================================================
@info "Implicit Solution With Units"
# =============================================================================================================
u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")

tspan = dstrip.((0.0u"min", 0.25u"min"))
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.FullSpecialize}(acceleration, u0, tspan, p, abstol=abstol, reltol=reltol)
sol = solve(prob, Rodas5P())
@btime soln = solve(prob, Rodas5P())
plt = plot!(plt, sol.t, [dstrip(u.v) for u in sol.u], label="implicit units") #Each element in sol.u is a QuantFieldVector
