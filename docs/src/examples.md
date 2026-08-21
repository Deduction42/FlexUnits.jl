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

## Custom dimensions
The default dimensions uses the SI unit system dimensions, but there is some contention around what constitutes a dimension. For example, angles are not a dimension, since they can be described as a ratio. Moreover, this system has no notion of currency, because it is not a *physical* dimension but merely a fuzzy human notion of value (hence why currencies fluctuate over time, usually downward due to inflation). This does not stop you, the user, from being able to include additional dimensions such as currency and angles. In this example below, we will be adding a currency dimension in the form of Euros.

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
While this is technically usable, most of the FlexUnits API makes use of string macros which look up units from a registry. Unfortunately, in order to make units type stable, all units in a registry must have the same (dynamic) dimension type. As shown in the previous RationalRegistry example, FlexUnits provides tooling to make building registries as painless as possible. Since `MoneyDimensions` is a generalization of `Dimensions`, `registry_defaults!` can be used to populate the registry with all the units normally inside `Dimensions`.

```julia
#===============================================================================================================================
Create your unit registry as a module using RegistryTools
===============================================================================================================================#
module CurrencyUnits

using FlexUnits
using FlexUnits.RegistryTools
import ..MoneyDimensions

const UNITS = PermanentDict{Symbol,Units{MoneyDimensions{FixRat32},AffineTransform{Float64}}}()

registry_defaults!(UNITS)
register_unit!(UNITS, "€" => MoneyDimensions{FixRat32}(€=1))
register_unit!(UNITS, "EUR" => UNITS[:€])

uparse(str::String) = RegistryTools.uparse(str, UNITS)
qparse(str::String) = RegistryTools.qparse(str, UNITS)

macro u_str(str); return suparse_expr(str, UNITS); end
macro ud_str(str); return uparse_expr(str, UNITS); end
macro q_str(str); return qparse_expr(str, UNITS); end
macro U_str(str); suexpr = suparse_expr(str, UNITS); return :($typeof($suexpr)); end
macro D_str(str); suexpr = suparse_expr(str, UNITS); return :($dimtype($suexpr)); end

utype() = RegistryTools.regunittype(UNITS)
dtype() = RegistryTools.regdimtype(UNITS)

export @u_str, @ud_str, @q_str, @U_str, @D_str, uparse, qparse, utype, dtype

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
Unfortunately, `simplify` must assume a set of units; if your units are not compatible with `Dimensions` (i.e. fields are not the same), simplify will not work.
```julia
electricity_price = 0.195u"€/(kW*hr)"
5.416666666666666e-8 (s² €)/(m² kg)

electricity_price = 0.195u"€/(kW*hr)" |> simplify
ArgumentError: type does not have a definite number of fields
```
In order to support simplification, you will need to create a constant vector of desired units (ordered by highest complexity first) and extend `FlexUnits.preferred_units(::Type{T})` to map your dimensional type to that vector.

```julia
#===============================================================================================================================
If you want simplify(...) to work with your units, 
    => define a constant vector of preferred units that are compatible with your new type
    => overload FlexUnits.preferred_units(::Type{T}) to refer to that list of preferred units
===============================================================================================================================#
const CURRENCY_UNITS = [ CurrencyUnits.UNITS[k] for k in [:F, :H, :T, :Ω, :V, :W, :J, :Pa, :N, :C, :L, :€] ]
FlexUnits.preferred_units(::Type{<:MoneyDimensions}) = CURRENCY_UNITS

electricity_price = 0.195u"€/(kW*hr)" |> simplify
5.416666666666666e-8 €/J
```

If you create a custom dimension set with a custom registry to go along with it, you will most likely only want to work with the dimension types in that registry. In the rare occasion that you will want to combine operations with the default unit registry (such as creating an extension to FlexUnits), you will need to define (1) promotion rules to prioritize your type (2) a conversion constructor to convert `Dimensions` to your new type
```julia
#===============================================================================================================================
If you want to combine operations with default unit registry
    => Add promotion rules and constructor to build your new type from Dimensions
===============================================================================================================================#
MoneyDimensions{T}(d::Dimensions) where T = MoneyDimensions{T}(m=d.m, kg=d.kg, s=d.s, A=d.A, K=d.K, cd=d.cd, mol=d.mol, €=zero(T))
Base.promote_rule(::Type{MoneyDimensions{T1}}, ::Type{Dimensions{T2}}) where {T1, T2} = MoneyDimensions{promote_type(T1,T2)}
Base.promote_rule(::Type{Dimensions{T1}}, ::Type{MoneyDimensions{T2}}) where {T1, T2} = MoneyDimensions{promote_type(T1,T2)}

(25*u"EUR") / (1*UnitRegistry.u"hr")
0.006944444444444444 €/s
```