#Preamble (delete when finished)

include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")
include("linalg_types.jl")


const ArrayOfDims{D} = Union{QuantArrayDims{D}, SMatrix{N,M,D} where {N,M}} where D<:AbstractDimensions

#======================================================================================================================
Operators on dimensions objects
======================================================================================================================#

#AbstractDimsMap only needs to check the first row and column for equality
function Base.:(==)(d1::AbstractDimsMap, d2::AbstractDimsMap)
    equal_element(ind) = d1[begin-1+ind[1], begin-1+ind[2]] == d2[begin-1+ind[1], begin-1+ind[2]]

    size(d1) == size(d2) || return false
    all(equal_element, 1:size(d1,1) .=> 1) || return false
    return all(equal_element, 1 .=> 1:size(d1,2))
end

#For matrices, all elements must be checked for equality
Base.:(==)(d1::ArrayOfDims, d2::AbstractDimsMap) = size(d1) == size(d2) && all(==, zip(d1,d2))
Base.:(==)(d1::AbstractDimsMap, d2::ArrayOfDims) = size(d1) == size(d2) && all(==, zip(d1,d2))

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

#Multiplying generic factorizations with dense matrices of dimensions
Base.:*(d1::AbstractDimsMap, d2::AbstractDimsMap) = canonical!(DimsMap(u_in = uinput(d2), u_out = uoutput(d1).*dotinv(d2.uoutput, d1.uinput)))
Base.:*(d1::ArrayOfDims, d2::AbstractDimsMap) = canonical!(DimsMap(u_in = uinput(d2), u_out = d1*uoutput(d2)))
Base.:*(d1::AbstractDimsMap, d2::ArrayOfDims) = canonical!(DimsMap(u_in = inv.(d2'*inv.(uinput(d1))), u_out = uoutput(d1)))

#Multiplying specific factorizations with single dimensions
Base.:*(dm::DimsMap{<:AbstractDimensions}, d::AbstractDimensions)    = canonical!(DimsMap(u_in = dm.u_in, u_out = dm.u_out.*d))
Base.:*(d::AbstractDimensions, dm::DimsMap{<:AbstractDimensions})    = canonical!(DimsMap(u_in = dm.u_in, u_out = dm.u_out.*d))
Base.:*(dm::RepDimsMap{<:AbstractDimensions}, d::AbstractDimensions) = canonical!(RepDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(d::AbstractDimensions, dm::RepDimsMap{<:AbstractDimensions}) = canonical!(RepDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(dm::SymDimsMap{<:AbstractDimensions}, d::AbstractDimensions) = canonical!(SymDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(d::AbstractDimensions, dm::SymDimsMap{<:AbstractDimensions}) = canonical!(SymDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))

#Division
Base.inv(d::ArrayOfDims) = inv(DimsMap(d))
Base.:/(d1::Union{AbstractDimsMap,ArrayOfDims}, d2::Union{AbstractDimsMap,ArrayOfDims}) = d1*inv(d2)
Base.:\(d1::Union{AbstractDimsMap,ArrayOfDims}, d2::Union{AbstractDimsMap,ArrayOfDims}) = inv(d1)*d2

#Matrix powers 
Base.:^(d::Union{AbstractDimsMap,ArrayOfDims}, p::Real) = RepDimsMap(d)^p 
Base.:^(d::RepDimsMap, p::Real) = RepDimsMap(u_scale = d.u_scale^p, u_in=d.u_in)

#Matrix exponentials and other functions that merely assert idempotence
assert_idempotent(d::RepDimsMap) = isone(d.u_scale) ? d : throw(ArgumentError("Cannot exponentiate dimension mapping unless it is idempotent"))
assert_idempotent(d::Union{AbstractDimsMap,ArrayOfDims}) = assert_idempotent(RepDimsMap(d))

for op in (:exp, :log)
    @eval Base.$op(d::Union{AbstractDimsMap,ArrayOfDims}) = assert_idempotent(d)
end

#======================================================================================================================
Define "q" linear algebra methods that are distinct from LinearAlgebra and don't cause dispatch issues
======================================================================================================================#
qinv(q::AbstractMatrix) = inv(LinmapQuant(DimsMap, q))
qadd(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) + dstrip(m2), dimension(m1) + dimension(m2))
qsub(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) - dstrip(m2), dimension(m1) - dimension(m2))
qmul(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) * dstrip(m2), dimension(m1) * dimension(m2))
qdiv(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) / dstrip(m2), dimension(m1) / dimension(m2))
qldiv(m1::AbstractMatrix, m2::AbstractMatrix) = LinmapQuant(dstrip(m1) \ dstrip(m2), dimension(m1) \ dimension(m2))
qpow(m::AbstractMatrix, p::Real) = LinmapQuant(dstrip(m)^p, RepDimsMap(dimension(m))^p)
qexp(m::AbstractMatrix) = LinmapQuant(exp(dstrip(m)), exp(dimension(m)))
qlog(m::AbstractMatrix) = LinmapQuant(log(dstrip(m)), log(dimension(m)))
qadjoint(m::AbstractMatrix) = LinmapQuant(adjoint(dstrip(m)), adjoint(dimesnion(m)))
qtranspose(m::AbstractMatrix) = LinmapQuant(transpose(dstrip(m)), adjoint(dimension(m)))

#Overload the base methods for pure LinmapQuant methods
Base.:+(m1::LinmapQuant, m2::LinmapQuant) = qadd(m1, m2)
Base.:-(m1::LinmapQuant, m2::LinmapQuant) = qsub(m1, m2)
Base.:*(m1::LinmapQuant, m2::LinmapQuant) = qmul(m1, m2)
Base.:/(m1::LinmapQuant, m2::LinmapQuant) = qdiv(m1, m2)
Base.:\(m1::LinmapQuant, m2::LinmapQuant) = qldiv(m1, m2)
Base.:^(m::LinmapQuant, p::Real) = qpow(m, p)
Base.:exp(m::LinmapQuant) = qexp(m)
Base.:log(m::LinmapQuant) = qlog(m)

#List of matrices we want to overload when using bivariate operations
const COMB_MATRIX_TYPES = [:Matrix, :Diagonal, :Hermitian, :Symmetric, :SymTridiagonal, :Tridiagonal, 
                            :UpperHessenberg, :SMatrix, :MMatrix, :SizedMatrix, :FieldMatrix]

#List out quantity matrix types we want to explicitly overload for univariate operations
const QUANT_MATRIX_TYPES = map(Symbol, String["Matrix{<:Quantity}", "Diagonal{<:Quantity}", "Hermitian{<:Quantity}", "Symmetric{<:Quantity}",
                            "SymTridiagonal{<:Quantity}", "Tridiagonal{<:Quantity}", "UpperHessenberg{<:Quantity}", "SMatrix{<:Any,<:Any,<:Quantity}", 
                            "MMatrix{<:Any,<:Any,<:Quantity}", "SizedMatrix{<:Any,<:Any,<:Quantity}", "FieldMatrix{<:Any,<:Any,<:Quantity}"])

#Apply the mixed methods with various kinds of matrices
for M in COMB_MATRIX_TYPES
    @eval Base.:+(m1::$M, m2::LinmapQuant) = quadd(m1, m2)
    @eval Base.:+(m1::LinmapQuant, m2::$M) = quadd(m1, m2)

    @eval Base.:-(m1::$M, m2::LinmapQuant) = qsub(m1, m2)
    @eval Base.:-(m1::LinmapQuant, m2::$M) = qsub(m1, m2)

    @eval Base.:*(m1::$M, m2::LinmapQuant) = qmul(m1, m2)
    @eval Base.:*(m1::LinmapQuant, m2::$M) = qmul(m1, m2)

    @eval Base.:/(m1::$M, m2::LinmapQuant) = qdiv(m1, m2)
    @eval Base.:/(m1::LinmapQuant, m2::$M) = qdiv(m1, m2)

    @eval Base.:\(m1::$M, m2::LinmapQuant) = qldiv(m1, m2)
    @eval Base.:\(m1::LinmapQuant, m2::$M) = qldiv(m1, m2)
end 

#Apply the quantity-specific methods on single-argument matrix functions
for M in QUANT_MATRIX_TYPES
    @eval Base.:^(m::$M, p::Real) = qpow(m, p)
    @eval Base.:exp(m::$M) = qexp(m)
    @eval Base.:log(m::$M) = qlog(m)
    @eval Base.:transpose(m::$M) = qtranspose(m)
    @eval Base.:adjoint(m::$M) = qadjoint(m)
end



#======================================================================================================================
Utility functions
======================================================================================================================#
#Dot products of dimensions and dotinv (i.e. dot(x, inv.(y)))
LinearAlgebra.dot(d1::AbstractDimensions, d2::AbstractDimensions) = d1*d2
dotinv(d1::AbstractDimensions, d2::AbstractDimensions) = d1/d2
function dotinv(d1::AbstractVector, d2::AbstractVector)
    length(d1) == length(d2) || throw(DimensionMismatch("Inputs had different lengths $((length(d1), length(d2)))"))
    return sum(dotinv, zip(d1, d2))
end


#======================================================================================================================
Shortcut multiplication strategies
*(U::DimTransform, M::AbstractMatrix) = U.u_out * (U.u_in'*M)
*(M::AbstractMatrix, U::DimTransform) = (M*U.u_out) * U.u_in'
*(U1::DimTransform, U2::DimTransform) = U1.u_out.*dot(U1.u_in, U2.u_out).*U2.u_in'
      = DimTransform(u_out=U1.u_out.*dot(U1.u_in, U2.u_out), u_in=U2.u_in) 
*(U1::ScaleDimTransform, U2::ScaleDimTransform) = DimTransform(u_scale=U1.u_scale*U2.u_scale, u_out=U1.u_out) iff U1.u_out == U2.u_out
======================================================================================================================#
using Test
import .UnitRegistry.@u_str
import .UnitRegistry.@ud_str

using StaticArrays
import Random
using Statistics

#Nonlinear map
@kwdef struct PumpInput{T} <: FieldVector{2,T}
    current :: T 
    voltage :: T
end

@kwdef struct PumpOutput{T} <: FieldVector{3,T}
    power :: T 
    pressure :: T
    flow :: T 
end

function pumpfunc(x::PumpInput)
    p = x.current*x.voltage*0.9   
    return PumpOutput(power = p, pressure = sqrt(p), flow = sqrt(p))
end
pumpfunc(x::AbstractVector) = pumpfunc(PumpInput(x))


@testset "Linear Mapping Basics" begin
    Random.seed!(1234)

    #Quck tests 
    u1 = SA[u"lbf*ft", u"kW", u"rpm"]
    u2 = SA[u"kg/s", u"m^3/hr", u"kW"]

    xm = SMatrix{3,3}(randn(3,3))
    qMraw = xm.*u2./u1'
    qM = LinmapQuant(DimsMap, qMraw)

    #Matrix inversion
    x = SVector{3}(randn(3)).*u1
    y = qM*x
    @test all(x .≈ inv(qM)*y)

    #Matrix transpose
    @test all(Matrix(y') .≈ Matrix(x'*qM'))

    #Square matrices
    Σ = cov(randn(20,3)*rand(3,3))
    x = randn(3).*u2

    #Symmetric matrix
    rS = Σ.*inv.(u2).*inv.(u2)'
    qS = LinmapQuant(SymDimsMap, rS)
    @test all(qS .≈ ubase.(rS))
    @test x'*(rS)*x ≈ x'*qS*x

    #Repeatable matrix
    rR = Σ.*u2.*inv.(u2)'
    qR = LinmapQuant(RepDimsMap, rR)
    @test all(qR .≈ ubase.(rR))
    @test all((rR*rR)*x .≈ (qR*qR)*x)

    #Nonlinear mapping
    pumpunits = UnitMap(PumpInput(current=u"A", voltage=u"V"), PumpOutput(power=u"W", pressure=u"Pa", flow=u"m^3/s"))
    upumpfunc = FunctionQuant(pumpfunc, pumpunits)
    qinput = PumpInput(current=500*u"mA", voltage=6u"V")
    @test all(upumpfunc(qinput) .≈ pumpfunc(ustrip.(uinput(pumpunits), qinput)).*uoutput(pumpunits))

end
