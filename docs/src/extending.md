# Extending FlexUnits
FlexUnits was built to be modular and extendable, and supports to major approaches for extension:
1. Using local modules to build registries with customized behaviours (FlexUnits is registry-agnostic)
2. Using packages with local modules to customize behaviours and be more opinionated about settings

## Custom registry tooling
A unit registry consists of a local module with the following attributes
1.  A dictionary containing units for each symbol (required)
2.  A list of preferred units for simplification (optional)
3.  Standard functions and macros to be exported (required)

FlexUnits provides the `RegistryTools` module to make it easier to build custom registries and export customized parsing and string macros without all the boilerplate code.

### The RegistryTools module
The FlexUnits.RegistryTools module provides functions and macros to make it easier to build registries. *Because this is the design intent, `RegistryTool`s should only be used inside a module that functions as a unit registry.* An example of a standard unit registry module is shown below.
```julia
module MyUnitRegistry
    using FlexUnits.RegistryTools
    
    #PermanentDict prevents changing units because when a string macro is called, it never updates again
    const UNITS = PermanentDict{Symbol,Units{Dimensions{FixRat32},AffineTransform{Float64}}}()
    registry_defaults!(UNITS) 

    #(Optional) define a custom set of preferred units
    const PREFERRED_UNITS = [UNITS[u] for u in [:F, :H, :T, :Ω, :V, :W, :J, :Pa, :N, :C, :L]]
    @generate_unit_simplifier(PREFERRED_UNITS)

    #Macro to generate the boilerplate code for registry exports
    @generate_registry_exports(UNITS)
end
```
Note that setting `PREFERRED_UNITS` is optional. It is only required if you change the dimension type, but it is still recommended as it provides decoupling from the default registry.

While this module exports a number of functions that are useful, most users only need to know about one special type (`PermanentDict`) and two macros (`@generate_registry_exports`, `@generate_unit_simplifier`). The rest of the macros and functions are primarily provided for the macros to use.

```@docs
RegistryTools.PermanentDict
RegistryTools.@generate_registry_exports
RegistryTools.@generate_unit_simplifier
```

### Unit registry customization examples

#### Exact conversions with Rational
This package defaults to using Float64 conversion factors to accomplish conversions. This often results in small but visually annoying round-off errors.
```julia
using FlexUnits, .UnitRegistry
julia> uconvert(u"°C", 32u"°F")
5.684341886080802e-14 °C

julia> uconvert(u"°C", 14u"°F")
-9.999999999999943 °C
```

You can change this behaviour by building a new registry with a different `AffineTransform` type that preserves exact rational expressions. This can be done with only a few lines of code.

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
using FlexUnits
using .RationalRegistry

julia> uconvert(u"°C", 32u"°F")
0//1 °C

julia> uconvert(u"°C", 14u"°F")
-10//1 °C
```
This can be used to modify many different behaviours if you don't agree with the design decisions of the default registry. FlexUnits is designed to be registry-agnostic.

#### Custom dimensions
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
However, as seen in the previous example, new unit registries are relatively painless to build. There are a couple of things we should note before starting. Since `MoneyDimensions` is a generalization of `Dimensions`, `registry_defaults!` can be used to populate the registry with all the units normally inside `Dimensions`. Simplification doesn't work with this dimension out of the box, so we will need to perform the additional step of configuring unit simplification.

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

## Custom Packages
FlexUnits design is focused around performance and flexibility, but many design decisions made have tradeoffs that might not be the best fit you. FlexUnits disables unit simplification by default, displaying the result in SI base units. Maybe you want unit simplification enabled. FlexUnits does not automatically export string macros because it's registry-agnostic, but maybe you want your package to be opinionated about the registry and export the macros automatically. Due to challenges around changing units, FlexUnits only registers units with unique international definitions. Maybe you want U.S. or U.K. customary units in the registry you want to export by default. FlexUnits doesn't contain constants for units and physical constants, maybe you want to provide some exported constants. These behaviours can be achieved by building your own package that extends FlexUnits. 

The example below is a module that can be used in a package which does the following:
1. Creates a new unit registry with some default registrations and customized preferred units
2. Creates `const` variables and exports them
3. Exports registry-specific parsing functions and macros
4. Automatically enables unit simplification and sets the desired unit simplification basis


```julia
module OpinionatedUnits

import FlexUnits: set_preferred_unit, simplify, display_simplified_units

# Export names from this package
export inch, mm, lb, kg, lbf, °F, °C, psi, ksi, MPa, kPa, bar, atm, R
export simplify, set_preferred_unit

# Export common macros and functions from internal registry (this package is opinionated)
export @u_str, @ud_str, @q_str, @D_str, uparse, qparse, register_unit

# Build internal unit registry with default settings
module InternalRegistry
    using FlexUnits.RegistryTools

    # Define the unit registry as an empty dictionary
    const UNITS = PermanentDict{Symbol,Units{Dimensions{FixRat32},AffineTransform{Float64}}}()

    # Add default units and register new ones to UNITS dictionary
    registry_defaults!(UNITS)
    register_unit!(UNITS, "kip" => 1000 * UNITS[:lbf])
    register_unit!(UNITS, "ksi" => 1000 * UNITS[:psi])
    register_unit!(UNITS, "atm" => 101.325 * UNITS[:kPa])
    register_unit!(UNITS, "mph" => UNITS[:mi] / UNITS[:hr])

    # Define preferred units
    const PREFERRED_UNITS = [UNITS[u] for u in [:F, :H, :T, :Ω, :V, :W, :J, :Pa, :N, :C, :L]]

    # Generate simplifiers and exports for defined units with included macros
    @generate_unit_simplifier(PREFERRED_UNITS)
    @generate_registry_exports(UNITS)
end 

# "use" the internal package to make functions and macros available in the current namespace
using .InternalRegistry

# Define selected units in namespace
const inch = u"inch" # Imperial Length
const mm = u"mm"     # Metric Length
const lb = u"lb"     # Imperial Mass
const kg = u"kg"     # Metric Mass
const lbf = u"lbf"   # Imperial Force
const °F = u"°F"     # Imperial Temperature
const Ra = u"Ra"     # Alternate Imperial Temperature
const °C = u"°C"     # Metric Temperature
const psi = u"psi"   # Imperial Pressure
const ksi = u"ksi"   # Alternate Imperial Pressure
const MPa = u"MPa"   # Metric Pressure
const kPa = u"kPa"   # Alternate Metric Pressure
const bar = u"bar"   # Alternate Metric Pressure
const atm = u"atm"   # Alternate Metric Pressure

# Define important constants in namespace
const R = 8.31446261815324u"J/(K*mol)" #Universal gas constant

# Set preferred units for simplification and turn unit simplification on at startup
# Note the `__init__()` function is required because these statements mutate global variables
function __init__()
    set_preferred_unit(inch)
    set_preferred_unit(lb)
    set_preferred_unit(lbf)
    set_preferred_unit(Ra)  
    set_preferred_unit(ksi)
    set_preferred_unit(lb/inch^3)

    display_simplified_units(true)  # Always convert to simple preferred units
end

end
```