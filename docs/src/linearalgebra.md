# Linear Algebra
FlexUnits provides broader support for mixed-unit operations than other packages to date. This is done in to major ways:
1. Using a sentinel value to denote unknown dimensions and yield a one-time bypass to unit checking
2. Factoring units with LinmapQuant object to separate matrices from units, and provide shortcut methods to infer units on various operations

## Unknown Dimensions

### Much of Julia codebase breaks when dimensions are unknown
Much of the core Julia ecosystem relies the pattern of initializing values with `zero(T::Number)`, initializing matrices with `zeros` or even accessing `zero` elements in a sparse array. This works when the dimensions are statically known, but fails when they are not.
```julia
julia> Diagonal(randn(3).*Unitful.u"m/s")[1,2] #Unitful.jl works when units are statically known
0.0 m s^-1

julia> Diagonal(randn(3).*[Unitful.u"m/s", Unitful.u"kg/hr", Unitful.u"kW"])[1,2] #Unitful.jl fails when units are not statically known
ERROR: ArgumentError: zero(Unitful.Quantity{Float64}) not defined.

julia> Diagonal(randn(3).*DynamicQuantities.u"m/s")[1,2] #DynamicQuantities works because it factors this out as a QuantityArray
0.0 m s⁻¹

julia> Diagonal(randn(3).*[DynamicQuantities.u"m/s", DynamicQuantities.u"kg/hr", DynamicQuantities.u"kW"]) #DynamicQuantities won't let you build a Diagonal matrix with different dimensions
ERROR: DimensionError: -0.18627823028424698 m s⁻¹ and DynamicQuantities.Quantity{Float64, DynamicQuantities.Dimensions{DynamicQuantities.FRInt32}}[-0.18627823028424698 m s⁻¹, -0.00021716734256511895 kg s⁻¹, -262.81003984968476 m² kg s⁻³] have incompatible dimensions

julia> sparse([1,2,3], [1,2,3], [DynamicQuantities.u"m/s", DynamicQuantities.u"kg/hr", DynamicQuantities.u"kW"])[1,2] #Sparse results in same problem as Unitful
ERROR: Cannot create an additive identity from `Type{<:DynamicQuantities.Quantity}`, as the dimensions are unknown. Please use `zero(::DynamicQuantities.Quantity)` instead.
```

### FlexUnits provides an unknown dimension 
FlexUnits works like Unitful when units are statically known, but when they are not, a sentinel value is used to denote unknown dimensions. This is done by setting every dimension element to its `typemax` value (which is a very improbable value of `2147483647//25200` for the default `FixRat32` type). Unknown dimensions are displayed as `?/?`
```julia
julia> Diagonal(randn(3).*u"m/s")[1,2] #FlexUnits correctly returns unitful result when units are statically known
0.0 m/s

julia> Diagonal(randn(3).*[u"m/s", u"kg/hr", u"kW"])[1,2] #FlexUnits returns an object with unknown dimensions 
0.0 ?/?

julia> Diagonal(randn(3).*[u"m/s", u"kg/hr", u"kW"])[1,2] + 5u"m/s" #Unit-validating operations like '+' bypass unit validation once
5.0 m/s
```
This solves the "additive identity" problem mentioned in the DynamicQuantities error. In general, operations that change dimensions (like multiplication, division, and exponentiation) propagate the unknown dimension flag value, while operations that validate units (addition, subtraction, max/min) will return a known dimension.

```julia
julia> zero(Quantity{Float64, Dimensions{FixRat32}}) #zero produces unknown units when fed a type
0.0 ?/?

julia> zero(Quantity{Float64, Units{Dimensions{FixRat32},AffineTransform{Float64}}}) #zero only works on dimensional quantities
ERROR: ArgumentError: This operation only supports dimensional quantities

julia> oneunit(Quantity{Float64, Dimensions{FixRat32}}) * 5u"kW" #Unknown dimensions propagate through multiplication, values are scaled to SI
5000.0 ?/?

julia> max(500*oneunit(Quantity{Float64, Dimensions{FixRat32}}), 2u"J") #Comparison is done on the SI scale and returns known units
500.0 (m² kg)/s²
```

### Mitigating risks of silently wrong results
One potential issue for unkown dimensions is the ability to silently retrun the wrong results. The main way this is mitigated is to only support unknown values on raw dimensions. Since all dimensional units (like SI) can be converted to each other without any scaling factors, the scale numerical results will always be consistent. The only kind of error that can happen is a wrongly bypassed dimensional validation. Since validation can only be bypassed once, initializers (which predate validation) are not going to be an issue.

```
julia> zero(Quantity{Float64, Dimensions{FixRat32}}) + 3u"km/hr" #This makes sense, initializer is ignorant
0.8333333333333334 m/s

julia> zero(Quantity{Float64, Dimensions{FixRat32}}) + 3u"km/hr" + 2u"kg" #Incorrect operations still get flagged
ERROR: DimensionError: (m/s, kg) have incompatible dimensions

julia> julia> Diagonal([1u"m/s", 1u"kg", 1u"mol"])*(ones(3).*u"s") #Matrix multiplication naturally works
3-element Vector{Quantity{Float64, Dimensions{FixRat32}}}:
       1.0 m
  1.0 (kg s)
 1.0 (s mol)
```

The main risk is therefor accidentally deleting known unit information. As seen in the eample above, this doesn't happen as long as the units are full-rank. Issues can still happen with very sparse matrices through, and results clearly show that unknown dimensions are present.
```julia
julia> sparse([1,3], [1,3], [1u"m/s", 1u"mol"])*(ones(3).*u"s")
3-element Vector{Quantity{Float64, Dimensions{FixRat32}}}:
       1.0 m
     0.0 ?/?
 1.0 (s mol)
```
When matrices are this sparse, unkown dimensions in the result are understandable. However if the position in the matrix can be associated with a dimension, this problem can be solved by wrapping the numerical sparse matrix inside a LinmapQuant (explained in the next subsection). 

The other risk from unknown dimensions is confusion between `one` and `oneunit`
```julia
julia> one(Quantity{Float64, Dimensions{FixRat32}}) * 2u"A" #Calling 'one' maintains unit information
2.0 A

julia> oneunit(Quantity{Float64, Dimensions{FixRat32}}) * 2u"A" #Calling 'oneunit' deletes unit information
2.0 ?/?
```
Careful understanding of the nuances between these two functions should be acquired before attempting to use these functions. Even with this potential for confusion, when testing against common Julia functions and patterns, introducing unkown dimensions solves far more problems than it creates.


## The Unit Factorization trick with DimsMap and LinmapQuant

### DimsMap factorization compresses the units as a mapping
In order for matrix multiplication and other linear algebra operations to be possible, the units need to be consistent. In order to be consistent, the matrix must be constructed to take a vector with units `u_in` and produce a vector of units `u_out`. To do this, the units of an entire N×M matrix can be described by the following factorization as a `DimsMap`
- A scalar dimension factor `u_fac`
- An M-1 vector of input units `u_in`
- An N-1 vector of output units `u_out`
The internal object `DimsMap` represents this factorization into base units (although it uses N and M vectors and forces the first element to be dimensionless).

### LinmapQuant separates unit factorization from the numbers
The `LinmapQuant` object separates a numerical matrix from a `DimsMap`. For vectors, the `VectorQuant` object simply has a vector of numbers and a vector of base units. This separation provides the following benefits:
1. Efficient, well-tested, pure-numerical matrix operations can be performed on the raw numerical values
2. Unit dimensions can be solved using efficient O(N) methods on DimsMap objects which only store (M+N+1) values instead of (M×N) values
Thus, in order to perform a matrix operation, one simply performs the operations on the pure numerical values, and then perform the accelerated counterpart operation on the `DimsMap`.

### DimsMap Matrix Multiplication
To solve the `DimsMap` problem for matrix multiplication `m*n`, one simply has to perform these steps:
1. Ensure that `u_in` for `m` is the inverse of the `u_out` for `n` (an O(N) operation)
2. Multiply the `u_fac` values of both matrices (an O(1) operation)
3. Construct a new `DimsMap` object with the `u_out` for `m`, the `u_in` for `n` and the `u_fac` from step 2
The entire `DimsMap` solution is O(N) where N is the number of columns in `m`.

### DimsMap Matrix Inversion
To solve the `DimsMap` problem for a matrix inversion of `m`, one simplyh has to perform these steps:
1. Invert the `u_fac` for `m` (an O(1) operation)
2. Construct a new `DimsMap` object from the `u_fac` result in step 1 and swapping the `u_in` and `u_out` vectors (an O(1) operation)
The entire `DimsMap` solution in this case is O(1)

### DimsMap Matrix Powers
Only certain unit structures can be raised to a power, because in order for an `m*n` multiplication to work the input units of `m` must be dimensionally proprotional to the output units of `n`. Since `DimsMap` factors out these proportions, one can simply check for equality between `u_in` and `u_out` of the two matrices. Thus for a matrix power, one simply has to:
1. Verify that `u_in` is equal to `u_out`
2. Raise the `u_fac` value by the desired power

### More operations coming
Potential shortcut methods for `DimsMap` objects abound and can likely be applied to factorizations like eigendecomposition and cholesky factorization. Broadcasting operations are also likely targets to be short-cutted. To test if a shortcut method currently exists, simply perform it on `LinmapQuant` objects or `VectorQuant` objects. If a `LinmapQuant` or a `VectorQuant` object is returned, a shortcut method has been implemented. Otherwise a generic fallback has been used. This package also has semi-shortcut methods for combined operations between factorized `LinmapQuant`/`VectorQuant` objects and generic vectors and matrices of quantities. These methods tend to be O(N²), which is still better than the likely O(N³) alternative, but the advantage is that the resturn value should be factorized, optimizing long chains of linear algebra operations where possible.