# Advanced Examples

## Solving Differential Equations
It is currently possible to solve systems of differential equations using Unitful.jl, but it has some drawbacks:
1. It requires [sectioning arrays with RecursiveArrayTools](https://docs.sciml.ai/DiffEqDocs/stable/features/diffeq_arrays/#Example:-Dynamics-Equations)
2. It only works with explicit ODE solvers
3. Linear algebra operations are less efficient

With FlexUnits, arrays with heterogeneous units are supported and linear algebra operations are efficient. Moreover, with a bit of extra work, it is now possible to use units with explict solvers like `Rodas5P`. Integration with DifferentialEquations.jl is still in its early stages, but so far, FlexUnits support for mixed-unit linear algebra has already proven to make itegration relatively easy.

### The ODE problem (falling object)
This example will model a falling object with the drag force equation.
```julia
using FlexUnits, .UnitRegistry
using OrdinaryDiffEq
using StaticArrays
using Plots
using BenchmarkTools
using LinearAlgebra

#Use named vectors for readability
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

#Convert dynamic units to static units for performance
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

#Main differntial equation (unit-agnostic)
function acceleration_raw(u, p, t)
    fd = -sign(u.v)*0.5*p.ρ*u.v^2*p.Cd*p.A
    dv = fd/p.m - p.g
    dh = u.v
    return FallingObjectState(v=dv, h=dh)
end

#Main differential equation with a static unit wrapper (for performance)
function acceleration_ustatic(u::AbstractVector{<:Quantity}, p::AbstractVector{<:Quantity}, t)
    du = acceleration_raw(ustatic(FallingObjectState(u)), ustatic(FallingObjectProps(p)), t)
    return FallingObjectState(du)
end
```

### Solving with an explict solver (Tsit5)
Solving this differential equation with an explicit solver like `Tsit5` is relatively straightforward.
```julia
u0 = FallingObjectState(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps(Cd=1.0u"", A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")

tspan = (0.0u"s", 10.0u"s")
prob = ODEProblem{false, OrdinaryDiffEq.SciMLBase.NoSpecialize}(acceleration_ustatic, u0, tspan, p, abstol=[1e-6, 1e-6], reltol=[1e-6, 1e-6])
sol = solve(prob, Tsit5())
plt = plot(ustrip.(sol.t), [ustrip(u.v) for u in sol.u], label="Tsit5")
```

### Solving with an implicit solver (Rodas5P)
Unfortunately, implicit ODE solvers are a bit more complicated. The first major challenge is the fact that implicit solvers require differentiation. Automatic differentiation can differntiate *through* quantities, but it can't differentiate *with respect to* quantities. This means that the automatic differentiator needs to differentiate through the function *without units*, and the return a final Jacobian *with units* as shown in the `wrapped_jacboan` function definition below.

``` julia
using ForwardDiff
function wrapped_jacobian(f, x, p, t)
    u_in  = dimension.(x)
    u_out = dimension.(f(x,p,t))

    function f_unitless(xn)
        return dstrip.(f(xn.*u_in, p, t))
    end
    
    return ForwardDiff.jacobian(f_unitless, dstrip.(x)) * DimsMap(u_in=u_in, u_out=u_out)
end

static_jac(u,p,t) = wrapped_jacobian(acceleration_ustatic, u, p, t)
```

In addition, we need to define some other parametrs to make `Rodas5P` work. The first is the `tgrad` parameter which differentiates the function with respect to time. Becuase this set of equations does not explicitly use time, the `tgrad` can simply be set to zero with units of `dx/s`, otherwise, we would neeed to use another wrapped gradient.
```julia
static_tgrad(u,p,t) = [0.0u"m/s^3", 0.0u"m/s^2"]
```

Finally, in order to properly make use of these customized gradient/Jacobian functions, we need to construct an `ODEFunction` object. In order to make it properly function however, we need to overwrite the default mass matrix
```julia
mass_matrix = I*1ud""
```

This is different from the default identy matrix `I` becasue it is unit-agnositc on the off-diagonals, where the default `I` assumes dimensionless values for all elements and causes unit validation failures when elements of `v` have different unit-dimensions. Using unit-agnostic elements ensures that `mass_matrix*v` always returns `v` regardless of its unit-dimensions. One can check this by verifying that off diagonals produce `?/?`.
```julia
julia> mass_matrix[1,2]
0.0 ?/?
```

With this completed, we can now create an appropriate `ODEFunction` object and use the typical steps to solve.
```julia
f_static = ODEFunction{false, OrdinaryDiffEq.SciMLBase.FullSpecialize}(acceleration_ustatic, 
    jac = static_jac, 
    tgrad = static_tgrad,
    mass_matrix = mass_matrix
)

u0 = FallingObjectState(v=0.0u"m/s", h=100u"m")
p  = FallingObjectProps(Cd=1.0u"", A=0.1u"m^2", ρ=1.0u"kg/m^3", m=50u"kg", g=9.81u"m/s^2")

tspan = (0.0u"s", 10.0u"s")
prob = ODEProblem(f_static, u0, tspan, p, abstol=[1e-6, 1e-6], reltol=[1e-6, 1e-6])

sol = solve(prob, Rodas5P())
plt = plot!(plt, ustrip.(sol.t), [ustrip(u.v) for u in sol.u], label="Rodas5P")
```

While moving from `Tsit5` to `Rodas5P` is significantly more effort with quantities than raw numbers, FlexUnits is the only package known to be capable of this functionality. Moreover, these steps, which are already baked into defaults for regular number, could be potentially baked into defaults through an extension on the future.

## Exact conversions with Rational
This package defaults to using Float64 conversion factors to accomplish conversions. This often results in small but visually annoying roundoff errors.
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
    const UNITS = PermanentDict{Symbol, Units{Dimensions{FixRat32}, AffineTransform{Rational{Int64}}}}() #Just change the AffineTrnasform type

    #Fill the UNITS registry with default values
    registry_defaults!(UNITS)

    #Ueses a ReentrantLock() on register_unit to prevent race conditions when multithreading
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

    #Registry is exported but these functions/macros are not (in case user wants their own verison)
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