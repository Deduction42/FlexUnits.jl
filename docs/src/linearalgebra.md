# Performance

## Why FlexUnits is good at linear algebra
In order for matrix multiplication and other linear algebra operations to be possible, the units need to be consistent. In order to be consistent, the matrix must be constructed to take a vector with units `u_in` and produce a vector of units `u_out`. To do this, the units of an entire N×M matrix can be described by the following factorization
- A scalar dimension factor `u_fac`
- An M-1 vector of input units `u_in`
- An N-1 vector of output units `u_out`

The internal object `DimsMap` represents this factorization into base units (although it uses N and M vectors and forces the first element to be dimensionless), and the `LinmapQuant` object separates a numerical matrix from a `DimsMap`. For vectors, the `VectorQuant` object simply has a vector of numbers and a vector of base units. This separation provides the following benefits:
1. Efficient, well-tested, pure-numerical matrix operations can be performed on the raw numerical values
2. Unit dimensions can be solved using efficient O(N) methods on DimsMap objects which only store (M+N+1) values instead of (M×N) values
Thus, in order to perform a matrix operation, one simply performs the operations on the pure numerical values, and then perform the accelerated counterpart operation on the `DimsMap`.

### Matrix Multiplication
To solve the `DimsMap` problem for matrix multiplication `m*n`, one simply has to perform these steps:
1. Ensure that `u_in` for `m` is the inverse of the `u_out` for `n` (an O(N) operation)
2. Multiply the `u_fac` values of both matrices (an O(1) operation)
3. Construct a new `DimsMap` object with the `u_out` for `m`, the `u_in` for `n` and the `u_fac` from step 2
The entire `DimsMap` solution is O(N) where N is the number of columns in `m`.

### Matrix Inversion
To solve the `DimsMap` problem for a matrix inversion of `m`, one simplyh has to perform these steps:
1. Invert the `u_fac` for `m` (an O(1) operation)
2. Construct a new `DimsMap` object from the `u_fac` result in step 1 and swapping the `u_in` and `u_out` vectors (an O(1) operation)
The entire `DimsMap` solution in this case is O(1)

### Matrix Powers
Only certain unit structures can be raised to a power, because in order for an `m*n` multiplication to work the input units of `m` must be dimensionally proprotional to the output units of `n`. Since `DimsMap` factors out these proportions, one can simply check for equality between `u_in` and `u_out` of the two matrices. Thus for a matrix power, one simply has to:
1. Verify that `u_in` is equal to `u_out`
2. Raise the `u_fac` value by the desired power

### More operations coming
Potential shortcut methods for `DimsMap` objects abound and can likely be applied to factorizations like eigendecomposition and cholesky factorization. Broadcasting operations are also likely targets to be short-cutted. To test if a shortcut method currently exists, simply perform it on `LinmapQuant` objects or `VectorQuant` objects. If a `LinmapQuant` or a `VectorQuant` object is returned, a shortcut method has been implemented. Otherwise a generic fallback has been used. This package also has semi-shortcut methods for combined operations between factorized `LinmapQuant`/`VectorQuant` objects and generic vectors and matrices of quantities. These methods tend to be O(N²), which is still better than the likely O(N³) alternative, but the advantage is that the resturn value should be factorized, optimizing long chains of linear algebra operations where possible.