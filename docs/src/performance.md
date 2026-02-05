# Performance

## Why FlexUnits is fast
FexUnits.jl combines techniques from Unitful.jl on "static units" to obtain near zero overhead performance when units can be resolved at parse time, but falls back to "dynamic unit" methods used by DynamicQuantities.jl if units cannot be inferred at parse time. This is done through promotion rules which convert static-dimension quantities to dynamic-dimension quantities if different dimensions are present. This retains the high-performance behaviour of Unitful.jl when units are known at compile time, but often falls back to the performance of DynanicQuantity.jl if they can't be inferred. In the first set of benchmarks, we see that FlexUnits.jl and DynamicQuantities.jl vastly outperform Unitful.jl (by more than 100x) when units cannot be inferred.
```julia
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
```julia
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

## Why FlexUnits is efficient at linear algebra
In order for matrix multiplication and other linear algebra operations to be possible, the units need to be consistent. In order to be consistent, the matrix must be constructed to take a vector with units `u_in` and produce a vector of units `u_out`. To do this, the units of an entire N×M matrix can be described by the following factorization
- A scalar dimension factor `u_fact`
- An M-1 vector of input units `u_in`
- An N-1 vector of output units `u_out`
This results in a matrix being described by N+M-1 values instead of NxM, this can result in a significant amount of compressions. Addiitonally, these factored units are stored separately from the matrix of numerical values. This factorization and separation provides the two main benefits:
1. Efficient, well-tested, pure-numerical matrix operations can be performed on the raw numerical values
2. Unit dimensions can be solved using efficient O(N) methods on DimsMap objects which only store (M+N+1) values instead of (M×N) values (DimsMap stores two more values than the bare-minimum for validation purposes).

These efficiency gains means that even if units must be dynamic (due to different units in the array), the overhead of resolving these units can be minimal. Let's compare mixed-unit matrix multiplication with differnt packages.

### Matrix Multiplication Benchmarks
The first example consists of multiplying a 200x4 matrix by a 4x4 matrix
```julia
#Use unitless matrices as a benchmark
Nr = 200
X = randn(Nr, 4)
M = rand(4,4)

#Construct unitful matrices
uu = [Unitful.u"kg/s", Unitful.u"kW", Unitful.u"rad/s", Unitful.u"N/m"]
ut = reshape(uu, 1, :)
Xu = X.*ut
Mu = inv.(uu) .* M .* inv.(ut)

#Construct DynamicQuantity matrices
udq = [DynamicQuantities.u"kg/s", DynamicQuantities.u"kW", DynamicQuantities.u"rad/s", DynamicQuantities.u"N/m"]
udqt = reshape(udq, 1, :)
Xdq = X.*udqt
Mdq = inv.(udq) .* M .* inv.(udqt)

#Construct LinmapQuant matrices
ufq = [UnitRegistry.u"kg/s", UnitRegistry.u"kW", UnitRegistry.u"rad/s", UnitRegistry.u"N/m"]
Xfq = LinmapQuant(X, UnitMap(u_out = UnitRegistry.u"", u_in = inv.(ufq)))
Mfq = LinmapQuant(M, UnitMap(u_out = inv.(ufq), u_in=ufq))


julia> @btime X*M #No units
  700.000 ns (3 allocations: 6.35 KiB)

julia> @btime Xu*Mu #Unitful, more than 500x slower
  395.700 μs (5603 allocations: 93.83 KiB)

julia> @btime Xdq*Mdq #DynamicQuantities, about 8x slower
  5.700 μs (3 allocations: 31.34 KiB)

julia> @btime Xfq*Mfq #LinmapQuant, almsot no overhead
  710.000 ns (4 allocations: 6.41 KiB)
```
The main reason why FlexUnits.jl has nearly no overhead is that only the inner product of the units between matrices is considered. Only the first 4-element row of X and the first column of M need to be compared. Unit inference does not touch the other 199 rows of X or the other 3 colums of M.

### Linear Regression Benchmarks
Linear regression benchmarks can only be compared between FlexUnits and the raw numerical methods. Neither Unitful.jl nor DynamicQuantities.jl can handle matrix inversions.
```julia
julia> Mu/Mu #Unitful fails at 'oneunit'
ERROR: MethodError: no method matching oneunit(::Type{Any})
This error has been manually thrown, explicitly, so the method may exist but be intentionally marked as unimplemented.

julia> Mdq/Mdq #DynamicQuantities also fails at 'oneunit'
ERROR: Cannot create a dimensionful 1 from `Type{DynamicQuantities.Quantity}` without knowing the dimensions. Please use `oneunit(::DynamicQuantities.Quantity)` instead.

julia> collect(Mfq)/collect(Mfq) #Matrices of FlexUnit quantities fails further down becuase LU-factorization logic eventually compares quantities of different units
ERROR: DimensionError: (s²/kg², s⁴/(m² kg²)) have incompatible dimensions

julia> Mfq/Mfq #FlexUnits LinmapQuant matrices actually work
4×4 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
          1.0           0.0 m²/s²                 0.0 1/kg           -0.0 1/s
     0.0 s²/m²               1.0   -3.29362e-21 s²/(m² kg)  -2.75549e-20 s/m²
        0.0 kg     0.0 (m² kg)/s²                     1.0           -0.0 kg/s
 6.32944e-16 s  -6.10423e-13 m²/s         2.87226e-16 s/kg               1.0
```

We can use linear algebra to complete the linear regression as follows:
```julia
Nr = 200
XY = randn(Nr, 6) * rand(6, 6)
X = [XY[:, begin:4] ones(Nr)]
Y = XY[:, 5:end]
Xu = LinmapQuant(X, UnitMap(u_out=UnitRegistry.u"", u_in=inv.([UnitRegistry.u"kg/s", UnitRegistry.u"kW", UnitRegistry.u"rad/s", UnitRegistry.u"N/m", UnitRegistry.u""])))
Yu = LinmapQuant(Y, UnitMap(u_out=UnitRegistry.u"", u_in=inv.([UnitRegistry.u"K", UnitRegistry.u"kPa"])))

julia> @btime (X'X)\(X'Y) #No units
  4.880 μs (12 allocations: 992 bytes)

julia> @btime Bu = (Xu'Xu)\(Xu'Yu) #LinmapQuant, about 1.3x slower
  6.400 μs (19 allocations: 1.66 KiB)
```
This time, the overhead from unit inference was noticeable, but still much less than 2x. The reason for this is because `Xu'*Xu` is a long multiplication that compares the 200 columns of `Xu'` vs the 200 rows of `Xu`, with the same thing happening again in `Xu'*Yu`. Thankfully, this 200-row comparison is required only once for the 25 combinations of `Xu'*Xu` and once for the 10 combinations of `Xu'*Yju`. In general, unit inference is about 6x to 8x slower than floating-point operations, so having two unit inferences for every 35 floating-point operations gives a slowdown factor of (35+2*7)/35 = 1.4, almost exactly the slowdown that was seen in the benchmarks.


## Performance Tips
While FlexUnits is generally fast, there are a few things one may need to watch out for to get the most out of this package.

### 1. Avoid containers with mixed static-unit types, use dynamic units instead
A great deal of the work done in this package was devoted to building a performant type-stable dynamic unit/dimension system that can represent many different unit types and promotion rules that avoid mixed-type containers. While promotion rules help, some Julia functions don't apply conversion (this includes `collect` and `map`, but `vcat(...)` reliably promotes); if such mixed-type containers occur, use `udynamic` to explicitly convert static units to dynamic ones.

### 2. Use dynamic units for high-level code
Dynanmic quantities are always type-stable and are less likely to result in accidental performance-killing dynamic dispatch calls. This is increasingly important if you want to produce small static binaries with Julia because dynamic dispatch can inhibit this. Using dynamic units can also prevent long compile times as it reduces specialization.

### 3. Use LinmapQuant/VectorQuant for linear algebra
These representations allow the use of optimized numerical methods for linear algebra, and employ shortcuts to solve the units of the matrix separately. This is especially important for large matrices. If input values are already numerical matrices, constructors for `LinmapQuant` are more efficient at *attaching units* to said vectors, as only N+M+1 dimension values are stored instead of the full M×N.

### 4. Pay attention to the output types of linear algebra operations
If a shortcut method is implemented, the output should be either a `LinmapQuant` or a `VectorQuant`. Otherwise a slower fallback method has been used. If this happens where not expected, please submit a bug report. It is likely that a method has been overlooked.

### 5. Use `dconvert` to transition from dynamic units to static units in low-level code
Dynamic units are great for achieving effortless type stability, but static units really shine in performance-sensitive low-level code where there's a small number of variables with known dimensions. In such cases, one can simply use `dconvert(u"...", q)` to convert `q` to a quantity with the same dimensions as `u`. Because most of the calcualtion is done using dimensional quantities, no conversion math needs to take place, one simply needs to verify that the units match, thus `dconvert` has very little overhead.

### 6. Use `vcat` or some other promoting method to transition from low-level code to high-level code
As mentioned before, some functions like `collect` and `map` don't use `promote`. However, other methods like explicit vector construction like `v=[x,y,z]` and using `vcat` on splatted tuples properly trigger `promote`. When returning containers with multiple units, make sure they are properly promoted to dynamic units, otherwise overspecialization and dynanmic calls may leak to other parts of your code.

