
const MatrixOfDims{D} = AbstractMatrix{<:AbstractDimensions}
const VectorOfDims{D} = AbstractVector{<:AbstractDimensions}

#======================================================================================================================
Operators on dimensions objects
======================================================================================================================#

#For matrices, all elements must be checked for equality
Base.:(==)(d1::AbstractDimsMap, d2::AbstractDimsMap) = ufactor(d1) == ufactor(d2) && uinput(d1) == uinput(d2) && uoutput(d1) == uoutput(d2)
Base.:(==)(d1::MatrixOfDims, d2::AbstractDimsMap) = size(d1) == size(d2) && all(i-> d1[i[1]] == d2[i[2]], zip(CartesianIndices(d1), CartesianIndices(d2)))
Base.:(==)(d1::AbstractDimsMap, d2::MatrixOfDims) = size(d1) == size(d2) && all(i-> d1[i[1]] == d2[i[2]], zip(CartesianIndices(d1), CartesianIndices(d2)))

equaldims(u1::AbstractDimsMap, u2::AbstractDimsMap) = (u1==u2) ? u1 : throw(DimensionError(u1, u2))
equaldims(u1::AbstractDimsMap, u2::MatrixOfDims) = (u1==u2) ? u1 : throw(DimensionError(u1, u2))
equaldims(u1::MatrixOfDims, u2::AbstractDimsMap) = (u1==u2) ? u2 : throw(DimensionError(u1, u2))

Base.:+(d1::AbstractDimsMap) = d1
Base.:+(d1::AbstractDimsMap, d2::AbstractDimsMap) = equaldims(d1, d2)
Base.:+(d1::AbstractDimsMap, d2::MatrixOfDims) = equaldims(d1, d2)
Base.:+(d1::MatrixOfDims, d2::AbstractDimsMap) = equaldims(d1, d2)

Base.:-(d1::AbstractDimsMap) = d1
Base.:-(d1::AbstractDimsMap, d2::AbstractDimsMap) = equaldims(d1, d2)
Base.:-(d1::AbstractDimsMap, d2::MatrixOfDims) = equaldims(d1, d2)
Base.:-(d1::MatrixOfDims, d2::AbstractDimsMap) = equaldims(d1, d2)

#Multiplying factored dimensions with dense matrices of dimensions
Base.:*(d1::AbstractDimsMap, d2::AbstractDimsMap) = DimsMap(u_in = uinput(d2), u_out = uoutput(d1), u_fac = ufactor(d1,d2))
Base.:*(d1::MatrixOfDims, d2::AbstractDimsMap) = DimsMap(u_in = uinput(d2), u_out = d1*uoutput(d2), u_fac = ufactor(d2))
Base.:*(d1::AbstractDimsMap, d2::MatrixOfDims) = DimsMap(u_in = inv.(d2'*inv.(uinput(d1))), u_out = uoutput(d1), u_fac = ufactor(d1))

#Multiplying factored dimensions with vectors of dimensions 
Base.:*(d1::AbstractDimsMap, d2::VectorOfDims) = uoutput(d1) .* (dotinv1(uinput(d1), d2)*ufactor(d1))
Base.:*(d1::Adjoint{<:AbstractDimensions, <:VectorOfDims}, d2::AbstractDimsMap) = ((d1*uoutput(d2))*ufactor(d2)./uinput(d1))'

#Multiplying specific factorizations with single dimensions
Base.:*(dm::DimsMap{<:AbstractDimensions}, d::AbstractDimLike) = DimsMap(u_in = dm.u_in, u_out = dm.u_out, u_fac = dm.u_fac*d)
Base.:*(d::AbstractDimLike, dm::DimsMap{<:AbstractDimensions}) = DimsMap(u_in = dm.u_in, u_out = dm.u_out, u_fac = dm.u_fac*d)

#Division of matrices
Base.:/(d1::AbstractDimsMap, d2::AbstractDimsMap) = d1*inv(d2)
Base.:/(d1::MatrixOfDims, d2::AbstractDimsMap) = d1*inv(d2)
Base.:/(d1::AbstractDimsMap, d2::MatrixOfDims) = d1*inv(DimsMap(d2))
Base.:\(d1::AbstractDimsMap, d2::AbstractDimsMap) = inv(d1)*d2
Base.:\(d1::MatrixOfDims, d2::AbstractDimsMap) = inv(DimsMap(d1))*d2
Base.:\(d1::AbstractDimsMap, d2::MatrixOfDims) = inv(d1)*d2

#Division of matrices and vectors
Base.:/(d1::VectorOfDims, d2::AbstractDimsMap) = d1*inv(d2)
Base.:\(d1::AbstractDimsMap, d2::VectorOfDims) = inv(d1)*d2

#Matrix powers 
Base.:^(d::AbstractDimsMap, p::Real) = unsafe_pow(assert_repeatable(d), p)
Base.:^(d::AbstractDimsMap, p::Integer) = unsafe_pow(assert_repeatable(d), p)
unsafe_pow(d::AbstractDimsMap, p::Real) = DimsMap(u_fac=ufactor(d)^p, u_in=uinput(d), u_out=uoutput(d))
unsafe_pow(d::AdjointDmap, p::Real) = adjoint(unsafe_pow(d.parent, p))

#Matrix exponentials and other functions that merely assert idempotence
for op in (:exp, :log)
    @eval Base.$op(d::AbstractDimsMap) = assert_idempotent(d)
end

#======================================================================================================================
Define "q" linear algebra methods that are distinct from LinearAlgebra and don't cause dispatch/ambiguitiy issues
======================================================================================================================#
#Matrices
qinv(q::AbstractMatrix) = inv(LinmapQuant(q))
qadd(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) + dstrip(m2), dimension(m1) + dimension(m2))
qsub(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) - dstrip(m2), dimension(m1) - dimension(m2))
qmul(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) * dstrip(m2), dimension(m1) * dimension(m2))
qmul(m::AbstractMatrix, q::QuantUnion) = LinmapQuant(dstrip(m) * dstrip(q), dimension(m) * dimension(q))
qmul(q::QuantUnion, m::AbstractMatrix) = LinmapQuant(dstrip(q) * dstrip(m), dimension(q) * dimension(m))
qdiv(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) / dstrip(m2), dimension(m1) / dimension(m2))
qldiv(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) \ dstrip(m2), dimension(m1) \ dimension(m2))
qpow(m::AbstractMatrix, p::Real) = LinmapQuant(dstrip(m)^p, DimsMap(dimension(m))^p)
qexp(m::AbstractMatrix) = LinmapQuant(exp(dstrip(m)), exp(DimsMap(dimension(m))))
qlog(m::AbstractMatrix) = LinmapQuant(log(dstrip(m)), log(DimsMap(dimension(m))))
qadjoint(m::AbstractMatrix) = adjoint(LinmapQuant(m))
qtranspose(m::AbstractMatrix) = transpose(LinmapQuant(m))
qisapprox(m1::AbstractMatrix, m2::AbstractMatrix) = dstrip(m1) ≈ dstrip(m2) && dimension(m1) == dimension(m2)

#Vectors
qadd(v1::AbstractVector, v2::AbstractVector) = VectorQuant(dstrip(v1) + dstrip(v2), dimension(v1) + dimension(v2))
qsub(v1::AbstractVector, v2::AbstractVector) = VectorQuant(dstrip(v1) - dstrip(v2), dimension(v1) - dimension(v2))
qmul(m::AbstractMatrix, v::AbstractVector) = VectorQuant(dstrip(m) * dstrip(v), dimension(m) * dimension(v))
qmul(vt::Adjoint{<:Any, <:AbstractVector}, m::AbstractMatrix) = qmul(m', qadjoint(vt))'
qldiv(m::AbstractMatrix, v::AbstractVector) = VectorQuant(dstrip(m) \ dstrip(v), dimension(m) \ dimension(v))
qdiv(vt::Adjoint{<:Any, <:AbstractVector}, m::AbstractMatrix) = qldiv(m', qadjoint(vt))'
qadjoint(v::AbstractVector) = adjoint(VectorQuant(v))
qadjoint(vt::Adjoint{<:Any, <:AbstractVector}) = VectorQuant(adjoint(vt))
qtranspose(v::AbstractVector) = transpose(VectorQuant(v))
qtranspose(vt::Transpose{<:Any, <:AbstractVector}) = transpose(VectorQuant(transpose(vt)))
qisapprox(v1::AbstractVector, v2::AbstractVector) = dstrip(v1) ≈ dstrip(v2) && dimension(v1) == dimension(v2)

#Factorizations
qdiv(mq::AbstractMatrix, fq::FactorQuant)  = LinmapQuant(dstrip(mq)/dstrip(fq), dimension(mq)/dimension(fq))
qldiv(fq::FactorQuant, mq::AbstractMatrix) = LinmapQuant(dstrip(fq)\dstrip(mq), dimension(fq)\dimension(mq))
qldiv(fq::FactorQuant, v::AbstractVector) = VectorQuant(dstrip(fq)\dstrip(v), dimension(fq)\dimension(v))
qdiv(vt::Adjoint{<:Any, <:AbstractVector}, fq::FactorQuant) = qldiv(fq', qadjoint(vt))'

#Overload the accelerated LinmapQuant/VectorQuant methods
Base.:+(m1::LinmapQuant, m2::LinmapQuant) = qadd(m1, m2)
Base.:+(v1::VectorQuant, v2::VectorQuant) = qadd(v1, v2)
Base.:-(m1::LinmapQuant, m2::LinmapQuant) = qsub(m1, m2)
Base.:-(v1::VectorQuant, v2::VectorQuant) = qsub(v1, v2)
Base.:*(m1::LinmapQuant, m2::LinmapQuant) = qmul(m1, m2)
Base.:*(m::LinmapQuant, v::VectorQuant)   = qmul(m, v)
Base.:*(vt::Adjoint{<:Any, <:VectorQuant}, m::LinmapQuant) = qmul(vt, m)
Base.:*(m::LinmapQuant, q::QuantUnion) = qmul(m, q)
Base.:*(q::QuantUnion, m::LinmapQuant) = qmul(q, m)
Base.:/(m1::LinmapQuant, m2::LinmapQuant) = qdiv(m1, m2)
Base.:/(vt::Adjoint{<:Any, <:VectorQuant}, m::LinmapQuant) = qdiv(vt, m)
Base.:\(m1::LinmapQuant, m2::LinmapQuant) = qldiv(m1, m2)
Base.:\(m::LinmapQuant, v::VectorQuant)   = qldiv(m, v)
Base.:^(m::LinmapQuant, p::Real) = qpow(m, p)
Base.:^(m::LinmapQuant, p::Integer) = qpow(m, p)
Base.:exp(m::LinmapQuant) = qexp(m)
Base.:log(m::LinmapQuant) = qlog(m)
Base.:(≈)(m1::LinmapQuant, m2::LinmapQuant) = qisapprox(m1, m2)
Base.:(≈)(v1::VectorQuant, v2::VectorQuant) = qisapprox(v1, v2)

#======================================================================================================================
Specify Base methods combingin AbstractMatrix/AbstractVector subtypes with LinmapQuant and VectorQuant
======================================================================================================================#

#List of matrices we want to overload when using bivariate operations
const COMB_MATRIX_TYPES = [Matrix, DenseMatrix, AbstractSparseMatrixCSC, Diagonal, Hermitian, Symmetric, SymTridiagonal, Tridiagonal, 
                            UpperHessenberg, SMatrix, MMatrix, SizedMatrix, FieldMatrix]

#List of vectors we want to overload when using bivariate operations
const COMB_VECTOR_TYPES = [Vector, DenseVector, AbstractCompressedVector, SVector, SizedVector, FieldVector]                       

#List out quantity matrix types we want to explicitly overload for univariate operations
const QUANT_MATRIX_TYPES = [:(Matrix{<:Quantity}), :(Diagonal{<:Quantity}), :(Hermitian{<:Quantity}), :(Symmetric{<:Quantity}),
                            :(SymTridiagonal{<:Quantity}), :(Tridiagonal{<:Quantity}), :(UpperHessenberg{<:Quantity}), :(SMatrix{<:Any,<:Any,<:Quantity}), 
                            :(MMatrix{<:Any,<:Any,<:Quantity}), :(SizedMatrix{<:Any,<:Any,<:Quantity}), :(FieldMatrix{<:Any,<:Any,<:Quantity})]

#Apply the mixed methods with various kinds of matrices
for MU in COMB_MATRIX_TYPES
    for M in [MU, Adjoint{<:Any, <:MU}, Transpose{<:Any, <:MU}] #Disambiguate Transpose and Adjoint
        @eval Base.:+(m1::$M, m2::LinmapQuant) = qadd(m1, m2)
        @eval Base.:+(m1::LinmapQuant, m2::$M) = qadd(m1, m2)
        @eval Base.:-(m1::$M, m2::LinmapQuant) = qsub(m1, m2)
        @eval Base.:-(m1::LinmapQuant, m2::$M) = qsub(m1, m2)

        @eval Base.:*(m1::$M, m2::LinmapQuant) = qmul(m1, m2)
        @eval Base.:*(m1::LinmapQuant, m2::$M) = qmul(m1, m2)
        @eval Base.:*(m::$M, v::VectorQuant) = qmul(m, v)
        @eval Base.:*(vt::Adjoint{<:Any, <:VectorQuant}, m::$M) = qmul(vt, m)

        @eval Base.:/(m1::$M, m2::LinmapQuant) = qdiv(m1, m2)
        @eval Base.:/(m1::LinmapQuant, m2::$M) = qdiv(m1, m2)
        @eval Base.:/(vt::Adjoint{<:Any, <:VectorQuant}, m::$M) = qdiv(vt, m)

        @eval Base.:\(m1::$M, m2::LinmapQuant) = qldiv(m1, m2)
        @eval Base.:\(m1::LinmapQuant, m2::$M) = qldiv(m1, m2)
        @eval Base.:\(m::$M, v::VectorQuant) = qldiv(m, v)

        @eval Base.:≈(m1::$M, m2::LinmapQuant) = qisapprox(m1, m2)
        @eval Base.:≈(m1::LinmapQuant, m2::$M) = qisapprox(m1, m2)
    end
end 

#Apply mixed methods with various kinds of vectors
for V in COMB_VECTOR_TYPES
    @eval Base.:+(v1::$V, v2::VectorQuant) = qadd(v1, v2)
    @eval Base.:+(v1::VectorQuant, v2::$V) = qadd(v1, v2)
    @eval Base.:-(v1::$V, v2::VectorQuant) = qsub(v1, v2)
    @eval Base.:-(v1::VectorQuant, v2::$V) = qsub(v1, v2)

    @eval Base.:*(v::Adjoint{<:Any, <:$V}, m::LinmapQuant) = qmul(v, m)
    @eval Base.:*(m::LinmapQuant, v::$V) = qmul(m, v)

    @eval Base.:/(v::Adjoint{<:Any, <:$V}, m::LinmapQuant) = qdiv(v, m)
    @eval Base.:\(m::LinmapQuant, v::$V) = qldiv(m, v)

    @eval Base.:≈(v1::$V, v2::VectorQuant) = qisapprox(v1, v2)
    @eval Base.:≈(v1::VectorQuant, v2::$V) = qisapprox(v1, v2)
end

#Apply the quantity-specific methods on single-argument matrix functions
for M in QUANT_MATRIX_TYPES
    @eval Base.:^(m::$M, p::Real) = qpow(m, p)
    @eval Base.:^(m::$M, p::Integer) = qpow(m, p)
    @eval Base.:exp(m::$M) = qexp(m)
    @eval Base.:log(m::$M) = qlog(m)
end

#Add FactorQuant methods 
Base.:/(mq::AbstractArray, fq::FactorQuant) = qdiv(mq, fq)
Base.:\(fq::FactorQuant, mq::AbstractArray) = qldiv(fq, mq)

#Special case ambiguities
Base.:\(m::Diagonal{T, SVector{N,T}}, v::VectorQuant) where {N,T} = qldiv(m, v)

#======================================================================================================================
Utility functions
======================================================================================================================#
"""
    ufactor(d1::AbstractDimsMap, d2::AbstractDimsMap)

Multiplies the first row of d1 with the first column of d2 in order to solve the multiplication scale
"""
ufactor(d1::AbstractDimsMap, d2::AbstractDimsMap) = ufactor(d1)*ufactor(d2)*dotinv1(uinput(d1), uoutput(d2))
ufactor(d1::AdjointDmap, d2::AbstractDimsMap) = ufactor(d1)*ufactor(d2)*dot(uoutput(d1.parent), uoutput(d2))
ufactor(d1::AbstractDimsMap, d2::AdjointDmap) = ufactor(d1)*ufactor(d2)*dotinv(uinput(d1), uinput(d2.parent))
ufactor(d1::AdjointDmap, d2::AdjointDmap)     = ufactor(d1)*ufactor(d2)*dotinv2(uoutput(d1.parent), uinput(d2.parent))


LinearAlgebra.dot(d1::AbstractDimensions, d2::AbstractDimensions) = d1*d2

#dot products with the second inverse
dotinv(d1::AbstractDimensions, d2::AbstractDimensions) = inv(d1)*inv(d2)
dotinv(d1::AbstractVector, d2::AbstractVector) = _sumfunc(dotinv, d1, d2)

dotinv1(d1::AbstractDimensions, d2::AbstractDimensions) = inv(d1)*d2
dotinv1(d1::AbstractVector, d2::AbstractVector) = _sumfunc(dotinv1, d1, d2)

dotinv2(d1::AbstractDimensions, d2::AbstractDimensions) = d1*inv(d2)
dotinv2(d1::AbstractVector, d2::AbstractVector) = _sumfunc(dotinv2, d1, d2)

#Dot products after applying a two-argument function
_sumfunc(f, v1::AbstractVector, v2::AbstractVector) = sum(f(x1, x2) for (x1, x2) in strictzip(v1, v2))



