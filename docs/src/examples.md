# Advanced Examples

## Solving Differential Equations
There has been some attempts at passing units through differential equations solvers, and while FlexUnits can do this it comes with certain drawbacks:
1. Passing units through ODE solvers does take a performance hit (often 2x the time)
2. Passing units through ODE solvers isn't rigorously tested by their developers and regressions can happen
3. If done improperly, it's difficult to troubleshoot errors

It also turns out that passing units through an ODE solver is often unnecessary. These solvers already produce correct results as long as the *units are coherent* (that is, they don't require scaling between dimensions). Unit/dimension errors often show up in setting up the ODE equations/functions, so it makes sense to restrict unit validation activity to this space.

In the section below, we introduce a strategy that can validate units inside differential equations, the benefits of this approach include:
1. It solves the system at *essentially zero added cost* 
2. It introduces no added complexity, using the ODE solvers as the developers intended (and test for )
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

First, we import necessary packages and define our object as QuantFieldVectors contain state values with static dimensions attached to them. The `@D_str` macro makes easy work of assigning static dimensions in an intuitive manner.
```julia
using FlexUnits, .UnitRegistry
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using StaticArrays
using Plots
using BenchmarkTools
using LinearAlgebra

import FlexUnits: QuantFieldVector, DimsMod
import OrdinaryDiffEqTsit5.SciMLBase.FullSpecialize

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

We can write out the equations as we would normally express them; when returning the value, we need to modify the dimensions by `1/s` for the derivative. The `ustrip` command ensures the output type is the same as the input type (a requirement for the ode solver)
```julia
function acceleration(u0::AbstractVector, p::FallingObjectProps, t)
    u = FallingObjectState(u0)

    #Drag force
    fd = -sign(u.v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    
    #Drag force effect on state (multiply by dt to make units work)
    dv = (fd/p.m - p.g)
    dh = u.v

    return ustrip(DimsMod{D"1/s"}(FallingObjectState, v=dv, h=dh))
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
prob = ODEProblem{false, FullSpecialize}(acceleration, u0, tspan, p, abstol=abstol, reltol=reltol)
sol = solve(prob, Tsit5())
plt = plot(sol.t, [dstrip(u.v) for u in sol.u], label="explicit") #Each element in sol.u is a QuantFieldVector

# =============================================================================================================
println("\nImplicit Solution")
# =============================================================================================================
u0 = FallingObjectState{Float64}(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps{Float64}(Cd=1.0, A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")

tspan = dstrip.((0.0u"min", 0.25u"min"))
prob = ODEProblem{false, FullSpecialize}(acceleration, u0, tspan, p, abstol=abstol, reltol=reltol)
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

You can change this behaviour by building a new registry with a different `AffineTransform` type that preserves exact rational expressions. FlexUnits is designed around interchangeable unit registries which are simply dictionaries that live inside a module. Due to how many problems can be solved by new unit registries, FlexUnits provides tools to build new registries with only a few lines of code.

```julia
using FlexUnits

module RationalRegistry
    using ..RegistryTools #RegistryTools contains all you need to build a registry in one simple import

    const UNITS = PermanentDict{Symbol, Units{Dimensions{FixRat32}, AffineTransform{Rational{Int64}}}}() #Just change the AffineTransform type
    registry_defaults!(UNITS) #Auto-populate the new registry

    @generate_registry_exports(UNITS) #Use macros to generate the boilerplate code for registry exports
end
```

That's it. We can restart the Julia, add this registry (Module), export all of the macros from our newly created `RationalRegistry`, and check out the new behaviour.
```julia
using .RationalRegistry

julia> uconvert(u"°C", 32u"°F")
0//1 °C

julia> uconvert(u"°C", 14u"°F")
-10//1 °C
```
This can be used to modify many different behaviours if you don't agree with the design decisions of the default registry. FlexUnits registries are designed to be interchangeable.

## Custom dimensions
The default dimensions uses the SI unit system dimensions, but there is some contention around what constitutes a dimension. For example, this unit system has no notion of currency, because it is not a *physical* dimension but merely a fuzzy human notion of value (hence why currencies fluctuate over time, usually downward due to inflation). This does not stop you from being able to include additional dimensions such as currency and angles. In this example below, we will be adding a currency dimension in the form of Euros.

```julia
#===============================================================================================================================
Define your custom dimensions that include Euros, note that the € symbol is a completely valid Julia variable name
===============================================================================================================================#
using FlexUnits 

@kwdef struct MoneyDimensions{P} <: AbstractDimensions{P}
    m   ::P = zero(FixRat32)
    kg  ::P = zero(FixRat32)
    s   ::P = zero(FixRat32)
    A   ::P = zero(FixRat32)
    K   ::P = zero(FixRat32)
    cd  ::P = zero(FixRat32)
    mol ::P = zero(FixRat32)
    €   ::P = zero(FixRat32)
end
MoneyDimensions(args::Real...) = MoneyDimensions{FixRat32}(args...)
```
While this is object is technically usable, most of the FlexUnits API makes use of string macros which look up units from a *type-stable* unit registry. Unfortunately, this new `MoneyDimensions` object is not compatible with the default registry which is built using `Dimensions`. Moreover, the `simplify` function requires looking at a list of preferred units defined in a registry, so unit simplification will not automatically work for this dimension type.
```julia
julia> simplify( 5*Dimensions(kg=1, m=1, s=-2))
5.0 N

julia> simplify( 5*MoneyDimensions(kg=1, m=1, s=-2))
ERROR: Function `preferred_units` not defined for type MoneyDimensions{FixRat32}: This function is usually defined in a unit registry. Perhaps there is no unit registry for MoneyDimensions{FixRat32} or its unit registry was not properly configured
```
However, as seen in the previous example, new unit registries are relatively painless to build. There are a couple of things we should note before starting. Since `MoneyDimensions` is a generalization of `Dimensions`, `registry_defaults!` can be used to populate the registry with all the units normally inside `Dimensions`. Simplification doesn't work with this dimension out of hte box, so we will need to perform the additional step of configuring unit simplification.

```julia
module CurrencyUnits

using FlexUnits.RegistryTools
import ..MoneyDimensions

const UNITS = PermanentDict{Symbol,Units{MoneyDimensions{FixRat32},AffineTransform{Float64}}}()
registry_defaults!(UNITS)
register_unit!(UNITS, "EUR" => UNITS[:€])

const PREFERRED_UNITS = [UNITS[u] for u in [:F, :H, :T, :Ω, :V, :W, :J, :Pa, :N, :C, :L]]

@generate_unit_simplifier(PREFERRED_UNITS)
@generate_registry_exports(UNITS)
end
```
This now gives you the ability to use string macros and look up units. Note that you must "use" this module instead of the default `UnitRegistry`.
```julia
using .CurrencyUnits

fuel_price = 2.04u"€/L"
driver_price = 20u"€/hr"
trip_speed = 5u"km/hr"
trip_distance = 100u"km"
fuel_consumption = 8.1u"L"/100u"km"
estimated_price = fuel_price*trip_distance*fuel_consumption + driver_price*trip_distance/trip_speed
416.524 €
```
Simplification also works as expected
```julia
julia> electricity_price = 0.195u"€/(kW*hr)"
5.416666666666666e-8 (s² €)/(m² kg)

julia> electricity_price = 0.195u"€/(kW*hr)" |> simplify
5.416666666666666e-8 €/J
```
