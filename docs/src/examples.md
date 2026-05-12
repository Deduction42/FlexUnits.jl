# Advanced Examples

## Solving Differential Equations
There has been some attempts at passing units through differential equations solvers, and while FlexUnits can do this it comes with certain drawbacks:
1. Passing units through ODE solvers does take a performance hit (often 2x the time)
2. Passing units through PDE solvers isn't rigorously tested by their developers and regressions can happen
3. If done improperly, it's difficult to troubleshoot errors

It also turns out that passing units through an ODE solver is often unnecessary. These solvers already produce correct results as long as the *units are coherent* (that is, they don't require scaling between dimensions). Unit/dimension errors often show up in setting up the ODE equations/functions, so it makes sense to restrict unit validation activity to this space.

In the section below, we introduce a strategy that can validate units inside differential equations, the benefits of this approach include:
1. It solves the system at *essentially zero added cost* 
2. It introduces no added complexity, using the ODE solvers as the developers intend and test for 
3. It validates units in a way that is easy to interpret

This strategy involves defining new objects that behave like a user-defined object with units inside the differential equation system, but behave like numerical vectors inside the differential equation solver. This can be done with the following steps

1. Define a modified `FieldVector` object (from StaticArrays) called `QuantFieldVector` that stores values with static dimensions
2. Override `getindex` so that it that produces raw numbers in base SI units when indexed like a vector (but still produces quantities with field access)
3. Define new objects with clear dimensional representations and subtype it to `QuantFieldVector`

Note that in the near future, FlexUnits may define its own `QuantFieldVector` object, simplifying this process.

#### The problem definition
Let us consider an application where we are using the drag force equation to model a falling object with velocity `v` and position `h`. This equation also requires us to consider the following parameters:
1. `Cd` the drag force coefficient
2. `A` the reference area
3. `ρ` the density of the fluid the object is falling through
4. `m` the mass of the object
5. `g` gravitational acceleration

First, we import necessary packages and define our new `QuantFieldVector` that returns scalars with `getindex`
```julia
using FlexUnits, .UnitRegistry
using OrdinaryDiffEq
using StaticArrays
using Plots
using BenchmarkTools
using LinearAlgebra

#Make "getindex" return a unitless version of the value
abstract type QuantFieldVector{N,T} <: FieldVector{N,T} end
Base.@propagate_inbounds Base.getindex(a::QuantFieldVector, i::Int) = dstrip(getfield(a, i))
```
We can then define the state values with static dimensions attached to them. The `@D_str` macro makes easy work of assigning static dimensions in an intuitive manner.

```julia
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
```

We can write out the equations as we would normally express them; when returning the value, we have to multiply by seconds because this is in differential form. This will validate that the output is in the correct units.
```julia
function acceleration_raw(u0::AbstractVector, p::FallingObjectProps, t)
    u = FallingObjectState(u0)
    dt = D"s"() #Time dimension

    #Drag force
    fd = -sign(u.v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    
    #Drag force effect on state (multiply by dt to make units work)
    dv = (fd/p.m - p.g)*dt
    dh = u.v*dt

    return FallingObjectState(v=dv, h=dh)
end
```

The trick with defining the system this way is that calling `FallingObjectState()` on a unitless vector will apply units to it, and getting indices from it will give you a unitless scalar in SI base units. 
```julia
julia> fos = FallingObjectState([1.0,1.0])
2-element FallingObjectState{Float64} with indices SOneTo(2):
 1
 1

julia> fos.v
1.0 m/s

julia> SVector(fos)
2-element SVector{2, Float64} with indices SOneTo(2):
 1.0
 1.0
```

#### Solving the problem
Defining the problem in this manner will give you internal unit verification of your equations (where you need it most) at *zero runtime cost* regardless as to whether or not you're using an implicit or explicit solver. You can simply solve the system the way you would normally do it

```julia
# =============================================================================================================
println("\nExplicit Solution")
# =============================================================================================================
u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")
abstol = FallingObjectState{Float64}(v=1e-6u"m/s", h=1e-6u"m")
reltol = SA[1e-6, 1e-6]

tspan = dstrip.((0.0u"min", 0.25u"min")) #Time span must be in seconds, dstrip takes care of this
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_raw, u0, tspan, p, abstol=abstol, reltol=reltol)
sol = solve(prob, Tsit5())
plt = plot(sol.t, [dstrip(u.v) for u in sol.u], label="explicit") #Each element in sol.u is a QuantFieldVector

# =============================================================================================================
println("\nImplicit Solution")
# =============================================================================================================
u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")

tspan = dstrip.((0.0u"min", 0.25u"min"))
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_raw, u0, tspan, p, abstol=abstol, reltol=reltol)
sol = solve(prob, Rodas5P())
plt = plot!(plt, sol.t, [dstrip(u.v) for u in sol.u], label="implicit") #Each element in sol.u is a QuantFieldVector
```


## Exact conversions with Rational
This package defaults to using Float64 conversion factors to accomplish conversions. This often results in small but visually annoying round-off errors.
```julia
using FlexUnits, .UnitRegistry
julia> uconvert(u"°C", 32u"°F")
5.684341886080802e-14 °C

julia> uconvert(u"°C", 14u"°F")
-9.999999999999943 °C
```

However, FlexUnits is designed to be registry-agnostic, with simply registry construction so this default Float64 conversion behaviour doesn't have to be the case (which is why it isn't exported by default). A user can simply copy-paste the "UnitRegistry.jl" file and modify one line of code that assigns `const UNITS` to use the transform type `AffineTransform{Rational{Int64}}` instead of `AffineTransform{Float64}`.

```julia
using FlexUnits

module RationalRegistry
    #RegistryTools contains all you need to build a registry in one simple import
    using ..RegistryTools

    const UNIT_LOCK = ReentrantLock()
    const UNITS = PermanentDict{Symbol, Units{Dimensions{FixRat32}, AffineTransform{Rational{Int64}}}}() #Just change the AffineTransform type

    #Fill the UNITS registry with default values
    registry_defaults!(UNITS)

    #Uses a ReentrantLock() on register_unit to prevent race conditions when multithreading
    register_unit(p::Pair{String, <:Any}) = lock(UNIT_LOCK) do
        register_unit!(UNITS, p)
    end

    #Parsing functions that don't require a dictionary argument
    uparse(str::String) = RegistryTools.uparse(str, UNITS)
    qparse(str::String) = RegistryTools.qparse(str, UNITS)

    #String macros are possible now that we are internally referring to UNITS
    macro u_str(str)
        return esc(suparse_expr(str, UNITS))
    end

    macro ud_str(str)
        return esc(uparse_expr(str, UNITS))
    end

    macro q_str(str)
        return esc(qparse_expr(str, UNITS))
    end

    #Functions to facilitate knowing types ahead of time, DO NOT EXPORT IF MULTIPLE REGISTRIES ARE USED
    unittype() = RegistryTools.unittype(UNITS)
    dimtype()  = RegistryTools.dimtype(UNITS)

    #Registry is exported but these functions/macros are not (in case user wants their own version)
    #You can import these by invoking `using .Registry`
    export @u_str, @ud_str, uparse, @q_str, qparse, register_unit
end
```

We can then export all of the macros from our newly created `RationalRegistry` model, and check out the new behaviour.
```julia
using .RationalRegistry

julia> uconvert(u"°C", 32u"°F")
0//1 °C

julia> uconvert(u"°C", 14u"°F")
-10//1 °C
```
This can be used to modify many different behaviours if you don't agree with the design decisions of the default registry. FlexUnits registries are designed to be truly modular and flexible.
