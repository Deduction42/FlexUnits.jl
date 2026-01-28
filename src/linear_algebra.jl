#Preamble (delete when finished)
#=
include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")
include("linalg_types.jl")
=#

const MatrixOfDims{D} = AbstractMatrix{<:AbstractDimensions}
const VectorOfDims{D} = AbstractVector{<:AbstractDimensions}


#======================================================================================================================
Operators on dimensions objects
======================================================================================================================#

#AbstractDimsMap only needs to check the first row and column for equality
function Base.:(==)(d1::AbstractDimsMap, d2::AbstractDimsMap)
    equal_element(ind) = d1[begin-1+ind[1], begin-1+ind[2]] == d2[begin-1+ind[1], begin-1+ind[2]]

    size(d1) == size(d2) || return false
    all(equal_element, 1:size(d1)[1] .=> 1) || return false
    return all(equal_element, 1 .=> 1:size(d1)[2])
end

#For matrices, all elements must be checked for equality
Base.:(==)(d1::MatrixOfDims, d2::AbstractDimsMap) = size(d1) == size(d2) && all(i-> d1[i[1]]==d2[i[2]], zip(CartesianIndices(d1), CartesianIndices(d2)))
Base.:(==)(d1::AbstractDimsMap, d2::MatrixOfDims) = size(d1) == size(d2) && all(i-> d1[i[1]]==d2[i[2]], zip(CartesianIndices(d1), CartesianIndices(d2)))

function equaldims(u1::AbstractDimsMap, u2::AbstractDimsMap)
    if u1 == u2
        return u1
    else
        throw(DimensionError(u1, u2))
    end
end

Base.:+(d1::AbstractDimsMap) = d1
Base.:+(d1::AbstractDimsMap, d2::AbstractDimsMap) = equaldims(d1, d2)
Base.:-(d1::AbstractDimsMap) = d1
Base.:-(d1::AbstractDimsMap, d2::AbstractDimsMap) = equaldims(d1, d2)

#Multiplying factored dimensions with dense matrices of dimensions
Base.:*(d1::AbstractDimsMap, d2::AbstractDimsMap) = canonical!(DimsMap(u_in = uinput(d2), u_out = uoutput(d1).*dotinv1(uinput(d1), uoutput(d2))))
Base.:*(d1::MatrixOfDims, d2::AbstractDimsMap) = canonical!(DimsMap(u_in = uinput(d2), u_out = d1*uoutput(d2)))
Base.:*(d1::AbstractDimsMap, d2::MatrixOfDims) = canonical!(DimsMap(u_in = inv.(d2'*inv.(uinput(d1))), u_out = uoutput(d1)))

#Multiplying factored dimensions with vectors of dimensions 
Base.:*(d1::AbstractDimsMap, d2::VectorOfDims) = uoutput(d1) .* dotinv1(uinput(d1), d2)
Base.:*(d1::Adjoint{<:AbstractDimensions, <:VectorOfDims}, d2::AbstractDimsMap) = (d2'*d1')'

#Multiplying specific factorizations with single dimensions
Base.:*(dm::DimsMap{<:AbstractDimensions}, d::AbstractDimensions)    = canonical!(DimsMap(u_in = dm.u_in, u_out = dm.u_out.*d))
Base.:*(d::AbstractDimensions, dm::DimsMap{<:AbstractDimensions})    = canonical!(DimsMap(u_in = dm.u_in, u_out = dm.u_out.*d))
Base.:*(dm::RepDimsMap{<:AbstractDimensions}, d::AbstractDimensions) = canonical!(RepDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(d::AbstractDimensions, dm::RepDimsMap{<:AbstractDimensions}) = canonical!(RepDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(dm::SymDimsMap{<:AbstractDimensions}, d::AbstractDimensions) = canonical!(SymDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(d::AbstractDimensions, dm::SymDimsMap{<:AbstractDimensions}) = canonical!(SymDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))

#Division of matrices
Base.:/(d1::MatrixOfDims, d2::AbstractDimsMap) = d1*inv(d2)
Base.:/(d1::AbstractDimsMap, d2::MatrixOfDims) = d1*inv(DimsMap(d2))
Base.:\(d1::MatrixOfDims, d2::AbstractDimsMap) = inv(DimsMap(d1))*d2
Base.:\(d1::AbstractDimsMap, d2::MatrixOfDims) = inv(d1)*d2

#Division of matrices and vectors
Base.:/(d1::VectorOfDims, d2::AbstractDimsMap) = d1*inv(d2)
Base.:\(d1::AbstractDimsMap, d2::VectorOfDims) = inv(d1)*d2

#Matrix powers 
Base.:^(d::Union{AbstractDimsMap,MatrixOfDims}, p::Real) = RepDimsMap(d)^p 
Base.:^(d::Union{AbstractDimsMap,MatrixOfDims}, p::Integer) = RepDimsMap(d)^p 
Base.:^(d::RepDimsMap, p::Real) = RepDimsMap(u_scale = d.u_scale^p, u_in=d.u_in)
Base.:^(d::RepDimsMap, p::Integer) = RepDimsMap(u_scale = d.u_scale^p, u_in=d.u_in)

#Matrix exponentials and other functions that merely assert idempotence
assert_idempotent(d::RepDimsMap) = isone(d.u_scale) ? d : throw(ArgumentError("Cannot exponentiate dimension mapping unless it is idempotent"))
assert_idempotent(d::Union{AbstractDimsMap,MatrixOfDims}) = assert_idempotent(RepDimsMap(d))

for op in (:exp, :log)
    @eval Base.$op(d::Union{AbstractDimsMap,MatrixOfDims}) = assert_idempotent(d)
end

#======================================================================================================================
Define "q" linear algebra methods that are distinct from LinearAlgebra and don't cause dispatch issues
======================================================================================================================#

#Matrices
qinv(q::AbstractMatrix) = inv(LinmapQuant(DimsMap, q))
qadd(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) + dstrip(m2), dimension(m1) + dimension(m2))
qsub(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) - dstrip(m2), dimension(m1) - dimension(m2))
qmul(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) * dstrip(m2), dimension(m1) * dimension(m2))
qdiv(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) / dstrip(m2), dimension(m1) / dimension(m2))
qldiv(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) \ dstrip(m2), dimension(m1) \ dimension(m2))
qpow(m::AbstractMatrix, p::Real) = LinmapQuant(dstrip(m)^p, RepDimsMap(dimension(m))^p)
qexp(m::AbstractMatrix) = LinmapQuant(exp(dstrip(m)), exp(dimension(m)))
qlog(m::AbstractMatrix) = LinmapQuant(log(dstrip(m)), log(dimension(m)))
qadjoint(m::AbstractMatrix) = adjoint(LinmapQuant(m))
qtranspose(m::AbstractMatrix) = transpose(LinmapQuant(m))

#Vectors
qadd(v1::AbstractVector, v2::AbstractVector) = VectorQuant(dstrip(v1) + dstrip(v2), dimension(m1) + dimension(m2))
qsub(v1::AbstractVector, v2::AbstractVector) = VectorQuant(dstrip(v1) - dstrip(v2), dimension(v1) - dimension(v2))
qmul(m::AbstractMatrix, v::AbstractVector) = VectorQuant(dstrip(m) * dstrip(v), dimension(m) * dimension(v))
qmul(vt::Adjoint{<:Any, <:AbstractVector}, m::AbstractMatrix) = qmul(m', qadjoint(vt))'
qldiv(m::AbstractMatrix, v::AbstractVector) = VectorQuant(dstrip(m) \ dstrip(v), dimension(m) \ dimension(v))
qdiv(vt::Adjoint{<:Any, <:AbstractVector}, m::AbstractMatrix) = qldiv(m', qadjoint(vt))'
qadjoint(v::AbstractVector) = adjoint(VectorQuant(v))
qadjoint(vt::Adjoint{<:Any, <:AbstractVector}) = VectorQuant(adjoint(vt))
qtranspose(v::AbstractVector) = transpose(VectorQuant(v))
qtranspose(vt::Transpose{<:Any, <:AbstractVector}) = transpose(VectorQuant(transpose(vt)))

#Factorizations
qdiv(mq::AbstractMatrix, fq::FactorQuant) = LinmapQuant(dstrip(mq)/dstrip(fq), dimension(mq)/dimension(fq))
qldiv(fq::FactorQuant, mq::AbstractMatrix) = LinmapQuant(dstrip(fq)\dstrip(mq), dimension(fq)\dimension(mq))

#Overload the accelerated LinmapQuant/VectorQuant methods
Base.:+(m1::LinmapQuant, m2::LinmapQuant) = qadd(m1, m2)
Base.:+(v1::VectorQuant, v2::VectorQuant) = qadd(v1, v2)
Base.:-(m1::LinmapQuant, m2::LinmapQuant) = qsub(m1, m2)
Base.:-(v1::VectorQuant, v2::VectorQuant) = qsub(v1, v2)
Base.:*(m1::LinmapQuant, m2::LinmapQuant) = qmul(m1, m2)
Base.:*(m::LinmapQuant, v::VectorQuant)   = qmul(m, v)
Base.:*(vt::Adjoint{<:Any, <:VectorQuant}, m::LinmapQuant) = qmul(vt, m)
Base.:/(m1::LinmapQuant, m2::LinmapQuant) = qdiv(m1, m2)
Base.:/(vt::Adjoint{<:Any, <:VectorQuant}, m::LinmapQuant) = qdiv(vt, m)
Base.:\(m1::LinmapQuant, m2::LinmapQuant) = qldiv(m1, m2)
Base.:\(m::LinmapQuant, v::VectorQuant)   = qldiv(m, v)
Base.:^(m::LinmapQuant, p::Real) = qpow(m, p)
Base.:^(m::LinmapQuant, p::Integer) = qpow(m, p)
Base.:exp(m::LinmapQuant) = qexp(m)
Base.:log(m::LinmapQuant) = qlog(m)

#List of matrices we want to overload when using bivariate operations
const COMB_MATRIX_TYPES = [Matrix, DenseMatrix, AbstractSparseMatrixCSC, Diagonal, Hermitian, Symmetric, SymTridiagonal, Tridiagonal, 
                            UpperHessenberg, SMatrix, MMatrix, SizedMatrix, FieldMatrix]

const COMB_VECTOR_TYPES = [Vector, DenseVector, AbstractCompressedVector, SVector, SizedVector, FieldVector]                       

#List out quantity matrix types we want to explicitly overload for univariate operations
const QUANT_MATRIX_TYPES = [:(Matrix{<:Quantity}), :(Diagonal{<:Quantity}), :(Hermitian{<:Quantity}), :(Symmetric{<:Quantity}),
                            :(SymTridiagonal{<:Quantity}), :(Tridiagonal{<:Quantity}), :(UpperHessenberg{<:Quantity}), :(SMatrix{<:Any,<:Any,<:Quantity}), 
                            :(MMatrix{<:Any,<:Any,<:Quantity}), :(SizedMatrix{<:Any,<:Any,<:Quantity}), :(FieldMatrix{<:Any,<:Any,<:Quantity})]

#Apply the mixed methods with various kinds of matrices
for M0 in COMB_MATRIX_TYPES
    for M in [M0, Adjoint{<:Any, <:M0}, Transpose{<:Any, <:M0}] #Disambiguate Transpose and Adjoint
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
end
Base.:\(m::Diagonal{T, SVector{N,T}}, v::FlexUnits.VectorQuant) where {N,T} = qldiv(m, v) #Random ambiguity

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


#======================================================================================================================
Utility functions
======================================================================================================================#
#Dot products of dimensions and dotinv (i.e. dot(x, inv.(y)))
LinearAlgebra.dot(d1::AbstractDimensions, d2::AbstractDimensions) = d1*d2
dotinv2(d1::AbstractDimensions, d2::AbstractDimensions) = d1*inv(d2)
function dotinv2(d1::AbstractVector, d2::AbstractVector)
    length(d1) == length(d2) || throw(DimensionMismatch("Inputs had different lengths $((length(d1), length(d2)))"))
    return sum(dotinv2(v1,v2) for (v1,v2) in zip(d1, d2))
end

#Dot products of dimensions and invdot (i.e. dot(inv.(x), y))
dotinv1(d1::AbstractDimensions, d2::AbstractDimensions) = inv(d1)*d2
function dotinv1(d1::AbstractVector, d2::AbstractVector)
    length(d1) == length(d2) || throw(DimensionMismatch("Inputs had different lengths $((length(d1), length(d2)))"))
    return sum(dotinv1(v1,v2) for (v1,v2) in zip(d1, d2))
end

