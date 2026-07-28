"""
    abstract type QuantFieldArray{N,T,D} <: FieldArray{N,T,D} end

Inheriting from this object will make it easier to define your own rank-D tensor types
with known units (only supports fundamental units).
```julia
@kwdef struct Stiffness{T<:Real} <: QuantFieldArray{Tuple{2,2,2}, Float64, 3}
    xxx::Quantity{T, D"N/m"}
    yxx::Quantity{T, D"N/m"}
    xyx::Quantity{T, D"N/m"}
    yyx::Quantity{T, D"N/m"}
    xxy::Quantity{T, D"N/m"}
    yxy::Quantity{T, D"N/m"}
    xyy::Quantity{T, D"N/m"}
    yyy::Quantity{T, D"N/m"}
end
```
Calling `getproperty(qa::QuantFieldArray, fn)` will produce a Quantity,
```julia
stiffness = Stiffness(0.01u"lbf/inch" .* ones(2,2,2))
stiffness.xxx
1.7512677165354331 kg/s²
```
while calling `getindex(qa::QuantFieldArray, ind)` will produce a pure numerical scalar in 
the coherent base units (fundamental SI units by default).
```julia
stiffness[1]
1.7512677165354331
```
This means that linear algebra operations will only "see" a numerical vector but engineering formulas 
that index by field will have (coherent) units attached to them ensuring unit correctness for such formulas. 
This results in near zero-overhead unit verification cost in applications like ODE solving
"""
abstract type QuantFieldArray{N, T, D} <: FieldArray{N, T, D} end

abstract type QuantFieldMatrix{N1, N2, T} <: QuantFieldArray{Tuple{N1, N2}, T, 2} end 

abstract type QuantFieldVector{N, T} <: QuantFieldArray{Tuple{N}, T, 1} end 

Base.getindex(qa::QuantFieldArray{N,T}, ind::Int) where {N,T} = convert(T, dstrip(getfield(qa, ind)))
fieldunit(::Type{T}, ind::Union{Symbol, Int}) where T<:QuantFieldArray = unittype(fieldtype(T, ind))()
fieldunits(::Type{QA}) where QA<: QuantFieldArray = map(Base.Fix1(fieldunit, QA), fieldnames(QA))


struct DimsMod{D<:StaticDims, A<:QuantFieldArray}
    parent :: A 
end

DimsMod{D}(arr::QA) where {D,QA} = DimsMod{D,A}(arr)
DimsMod{D}(::Type{QA}; kwargs...) where {D,QA} = DimsMod{D}(dimsmod(D, QA; kwargs...))
ustrip(dm::DimsMod) = getfield(dm, :parent)


@generated function dimsmod(::Type{D}, ::Type{QA}; kwargs...) where {D<:StaticDims, QA<:QuantFieldArray}
    fns = fieldnames(QA)
    fus = fieldunits(QA).*D()
    expr = :($QA())
    expr_qa = expr.args

    for (fn, fu) in zip(fns, fus)
        expr_kw = :($(Expr(:kw, fn, :(ustrip($fu, kwargs[$(Expr(:quote, fn))])))))
        push!(expr_qa, expr_kw)
    end

    return expr
end

Base.getproperty(dm::DimsMod{D}, fn::Symbol) where D = getproperty(ustrip(dm), fn)*D()
Base.getindex(dm::DimsMod, ind::Int) = getindex(ustrip(dm), ind)

