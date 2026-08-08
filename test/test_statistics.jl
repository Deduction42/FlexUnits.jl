using FlexUnits, .UnitRegistry
using Statistics 
using LinearAlgebra
using Random

@testset "Statistics" begin 
    Random.seed!(1234)

    #Generate a correlated test set
    UX  = [ud"kg/s", ud"kW", ud"Hz"]
    UY  = [ud"kPa", ud"K"]
    X   = (rand(30,3)*rand(3,3) .+ 0.001.*randn(30,3))
    Y   = X*rand(3,2) .+ 0.001.*randn(30,2)
    QX  = X * UnitMap(u_in=inv.(UX), u_out=u"")
    QY  = Y * UnitMap(u_in=inv.(UY), u_out=u"")

    #Test summations, stats and some linear algebra
    @test sum(x->1, QX, dims=1) isa Matrix{<:Integer}
    @test sum(x->1, QX, dims=2) isa Matrix{<:Integer}
    @test sum(x->1, QX, dims=:) isa Integer
    @test FlexUnits.ureduce(identity, dimension(QX), dims=1) == dimension(sum(QX, dims=1))
    @test FlexUnits.ureduce(identity, dimension(QX'), dims=2) == dimension(sum(QX', dims=2))
    @test QX .+ QX .+ QX isa LinmapQuant

    @test sum(QX, dims=1) isa LinmapQuant
    @test all(sum(QX, dims=1) .≈ sum(X, dims=1).*UX')
    @test all(sum(QX', dims=2) .≈ sum(X', dims=2).*UX)

    @test mean(QX, dims=1) isa LinmapQuant
    @test all(mean(QX, dims=1) .≈ mean(X, dims=1).*UX')
    @test all(mean(QX', dims=2) .≈ mean(X', dims=2).*UX)
    @test mean(abs2, QX, dims=1) isa LinmapQuant
    @test all(mean(abs2, QX, dims=1) .≈ mean(abs2, X, dims=1).*(UX.^2)')

    @test median(QX, dims=1) isa LinmapQuant
    @test all(median(QX, dims=1) .≈ median(X, dims=1).*UX')
    @test all(median(QX', dims=2) .≈ median(X', dims=2).*UX)

    @test var(QX, dims=1) isa LinmapQuant
    @test all(var(QX, dims=1) .≈ var(X, dims=1).*(UX.^2)')
    @test all(varm(QX, mean(QX, dims=1), dims=1) .≈ varm(X, mean(X, dims=1), dims=1).*(UX.^2)')
    @test all(var(QX, mean=mean(QX, dims=1), dims=1) .≈ var(X, mean=mean(X, dims=1), dims=1).*(UX.^2)')
    

    @test std(QX, dims=1) isa LinmapQuant
    @test all(std(QX, dims=1) .≈ std(X, dims=1).*UX')
    @test all(stdm(QX, mean(QX, dims=1), dims=1) .≈ stdm(X, mean(X, dims=1), dims=1).*UX')
    @test all(std(QX, mean=mean(QX, dims=1), dims=1) .≈ std(X, mean=mean(X, dims=1), dims=1).*UX')

    @test cov(QX) isa LinmapQuant
    @test all(cov(QX) .≈ cov(X).*UX.*UX')
    @test cov(QX, QY) isa LinmapQuant
    @test all(cov(QX, QY) .≈ cov(X,Y).*(UX.*UY'))
    @test all(cov(QX', dims=2) .≈ cov(QX, dims=1))
    @test_throws ArgumentError cov(QX, dims=3)

    @test cor(QX) isa Matrix{<:Real}
    @test all(cor(QX) .≈ cor(X))
    @test cor(QX, QY) isa Matrix{<:Real}
    @test all(cor(QX, QY) .≈ cor(QX, QY))

    QZ = X*UnitMap(u_in=u"", u_out=u"") 
    @test sum(QZ) ≈ sum(X)
    @test sum(identity, QZ) ≈ sum(X)
    @test FlexUnits.ureduce(dimension(QZ)) == D""()
    @test FlexUnits.ureduce(x->x*D"m"(), dimension(QZ)) == D"m"()
    @test all(minimum(QX, dims=1, init=typemax(eltype(QX))) .≈ minimum(X, dims=1).*UX')
    @test all(maximum(QX, dims=1, init=typemin(eltype(QX))) .≈ maximum(X, dims=1).*UX')

    #Test shortcut methods 
    MQ = X * UnitMap(u_out=u"", u_in=inv.(UX))
    @test sum(MQ, dims=1) isa LinmapQuant
    @test all(sum(MQ, dims=1) .≈ sum(X, dims=1).*UX')
    @test minimum(MQ, dims=1) isa LinmapQuant
    @test all(minimum(MQ, dims=1) .≈ minimum(X, dims=1).*UX')
    @test maximum(MQ, dims=1) isa LinmapQuant
    @test all(maximum(MQ, dims=1) .≈ maximum(X, dims=1).*UX')
    @test all(sum(MQ', dims=2) .≈ sum(X', dims=2).*UX)
end
