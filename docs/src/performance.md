# Performance

## Why FlexUnits is fast
FexUnits.jl combines techniques from Unitful.jl on "static units" to obtain near zero overhead performance when units can be resolved at parse time, but falls back to "dynamic unit" methods used by DynamicQuantities.jl if units cannot be inferred at parse time. This is done through promotion rules which convert static-dimension quantities to dynamic-dimension quantities if different dimensions are present. This retains the high-performance behaviour of Unitful.jl when units are known at compile time, but often falls back to the performance of DynanicQuantity.jl if they can't be inferred. In the first set of benchmarks, we see that FlexUnits.jl and DynamicQuantities.jl vastly outperform Unitful.jl (by more than 100x) when units cannot be inferred.
```
using FlexUnits
using .UnitRegistry
import DynamicQuantities
import Unitful
using BenchmarkTools

v1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"J/kg", 1.0*Unitful.u"A/V"]
v1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"J/kg", 1.0*DynamicQuantities.u"A/V"]
v1flex = [1.0u"m/s", 1.0u"J/kg", 1.0u"A/V"]

@btime sum(x->x^0.0, $v1uni)
  7.575 μs (86 allocations: 3.92 KiB)
@btime sum(x->x^0.0, $v1dyn)
  41.667 ns (0 allocations: 0 bytes)
@btime sum(x->x^0.0, $v1flex)
  27.209 ns (0 allocations: 0 bytes)
```
In the second example, we see that FlexUnits.jl and Unitful.jl outperform DynanicQuantities.jl when units can be inferred by the compiler.
```
t1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s"]
t1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s"]
t1flex = [1.0u"m/s", 1.0u"m/s", 1.0u"m/s"]

@btime sum(x->x*x, $t1uni)
  2.900 ns (0 allocations: 0 bytes)
@btime sum(x->x*x, $t1dyn)
  7.800 ns (0 allocations: 0 bytes)
@btime sum(x->x*x, $t1flex)
  2.900 ns (0 allocations: 0 bytes)
```
While this performance boost over DynamicQuantities.jl isn't as dramatic as the previous boost over Unitful.jl, it is still significant. In most benchmarks (examples can be found in the test folder of the FlexUnits repo) FlexUnits matches the best performing alternative (DynamicQuantities or Unitful). There are two notable exceptions:
1. FlexUnits.jl is slightly slower than Unitful.jl at `uconvert`, but still much faster than DynamicQuantities.jl
2. FlexUnits.jl can be significantlhy faster than both Unitful.jl and DynamicQuantities.jl for iterative statically-inferred algorithms where unitful values are reassigned (Unitful tends to overspecialize but the FlexUnits design avoids this)

## Why FlexUnits is good at linear algebra
In order for matrix multiplication and other linear algebra operations to be possible, the units need to be consistent. In order to be consistent, the matrix must be constructed to take a vector with units `u_in` and produce a vector of units `u_out`. To do this, the units of an entire N×M matrix can be described by the following factorization
- A scalar dimension factor `u_fact`
- An M-1 vector of input units `u_in`
- An N-1 vector of output units `u_out`
This separation provides the following benefits:
1. Efficient, well-tested, pure-numerical matrix operations can be performed on the raw numerical values
2. Unit dimensions can be solved using efficient O(N) methods on DimsMap objects which only store (M+N+1) values instead of (M×N) values.


## Performance Tips
While FlexUnits is generally fast, there are a few things one may need to watch out for to get the most out of this package.

### 1. Avoid containers with mixed static-unit types, use dynamic units instead
A great deal of the work done in this package was devoted to building a performant type-stable dynamic unit/dimension system that can represent many different unit types and promotion rules that avoid mixed-type containers. While promotion rules help, some Julia functions don't apply conversion (this includes `collect` and `map`, but `vcat(...)` reliably promotes); if such mixed-type containers occur, use `udynamic` to explicitly convert static units to dynamic ones.

### 2. Use dynamic units for high-level code
Dynanmic quantities are always type-stable and are less likely to result in accidental performance-killing dynamic dispatch calls. This is increasingly important if you want to produce small static binaries with Julia because dynamic dispatch can inhibit this. Using dynamic units can also prevent long compile times as it reduces specialization.

### 3. Use LinmapQuant/VectorQuant for linear algebra
These representations allow the use of optimized numerical methods for linear algebra, and employ shortcuts to solve the units of the matrix separately. This is especially important for large matrices. If input values are already numerical matrices, constructors for `LinmapQuant` are more efficient at *attaching units* to said vectors, as only N+M+1 dimension values are stored instead of the full M×N.

### 4. Pay attention to the output types of LinmapQuant/VectorQuant linear algebra
If a shortcut method is implemented, the output should be either a `LinmapQuant` or a `VectorQuant`. Otherwise a slower fallback method has been used. If this happens where not expected, please submit a bug report. It is likely that a method has been overlooked.

### 5. Use `dconvert` to transition from dynamic units to static units in low-level code
Dynamic units are great for achieving effortless type stability, but static units really shine in performance-sensitive low-level code where there's a small number of variables with known dimensions. In such cases, one can simply use `dconvert(u"...", q)` to convert `q` to a quantity with the same dimensions as `u`. Because most of the calcualtion is done using dimensional quantities, no conversion math needs to take place, one simply needs to verify that the units match, thus `dconvert` has very little overhead.

### 6. Use `vcat` or some other promoting method to transition from low-level code to high-level code
As mentioned before, some functions like `collect` and `map` don't use `promote`. However, other methods like explicit vector construction like `v=[x,y,z]` and using `vcat` on splatted tuples properly trigger `promote`. When returning containers with multiple units, make sure they are properly promoted to dynamic units, otherwise overspecialization and dynanmic calls may leak to other parts of your code.

