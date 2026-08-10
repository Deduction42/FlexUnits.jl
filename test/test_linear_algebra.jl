using FlexUnits, .UnitRegistry
using LinearAlgebra
using Statistics
using StaticArrays
using Random

@testset "Linear Algebra" begin

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

    Random.seed!(1234)

    #Generate a correlated test set
    U = [ud"kg/s", ud"kW", ud"Hz"]
    X = (rand(30,3)*rand(3,3) .+ 0.001.*randn(30,3))
    Q = ubase.(X .* U')
    S = cov(Q)

    #Test Eigenvalue decompsoition
    R = dconvert.(u"", cor(Q)) #PCA must use correlation matrix
    eig = eigen(R) #symmetric eigendecomposition
    @test eig isa FactorQuant
    @test all(R .≈ (eig.vectors * Diagonal(eig.values) * eig.vectors'))

    a = rand(3,3)
    A = (a'*a) * UnitMap(u_in=U, u_out=U.*u"s")
    eig = eigen(A) #repeatable eigendecomposition
    @test all(A^2 .≈ eig.vectors * Diagonal(eig.values.^2) * inv(eig.vectors))

    #Test cholesky decomposition
    ch  = cholesky(S)
    @test ch isa FactorQuant
    @test all(S .≈ (ch.L * ch.U))

    #DimsMap constructor with units 
    dm = DimsMap(u_fac=u"", u_in=[u"m/s", u"K", u"Pa"], u_out=[u"N"])
    @test dm == DimsMap(u_fac=dimension(u""), u_in=dimension.([u"m/s", u"K", u"Pa"]), u_out=dimension.([ud"N"]))
    @test axes(dm) == (Base.OneTo(1), Base.OneTo(3))
    @test Base.IndexStyle(typeof(dm)) == Base.IndexStyle(typeof(dm'))
    @test axes(dm') == (Base.OneTo(3), Base.OneTo(1))
    @test axes(dm,1) == Base.OneTo(1)
    @test axes(dm,2) == Base.OneTo(3)
    @test axes(dm',1) == Base.OneTo(3)
    @test axes(dm',2) == Base.OneTo(1)
    @test +(dm) == dm 
    @test -dm == dm

    #Quick linear algebra tests 
    u1 = SA[u"lbf*ft", u"kW", u"rpm"]
    u2 = SA[u"kg/s", u"m^3/hr", u"kW"]

    xm = SMatrix{3,3}(randn(3,3))
    qMraw = xm.*(u2./u1')
    qM = LinmapQuant(qMraw)
    x = SVector{3}(randn(3)).*u1
    y = qM*x
    d = dimension(qM)

    @test all(d .== collect(d))
    @test all(qM .== collect(qM))

    #Test various constructors 
    @test all(qM .≈ LinmapQuant(xm, UnitMap(u_in=u1, u_out=u2)))
    qC = xm.*u1'
    @test all(qC .≈ LinmapQuant(xm, UnitMap(u_in=inv.(u1), u_out=u"")))
    @test all(qC .≈ LinmapQuant(xm, UnitMap(u_in=inv.(u1), u_out=ud"")))
    @test all(qC .≈ LinmapQuant(xm, UnitMap(u_in=inv.(u1), u_out=fill(u"", length(u1)))))
    @test all(qC .≈ LinmapQuant(Matrix(xm), UnitMap(u_in=inv.(u1), u_out=fill(u"", length(u1)))))
    @test !FlexUnits.ArrayInterface.can_setindex(dimension(Matrix(qMraw)))
    @test !FlexUnits.ArrayInterface.can_setindex(dstrip(Matrix(qMraw)))
    @test ubase(qM) == qM

    #Test alternate constructors
    @test qM ≈ LinmapQuant(xm, UnitMap(u_in = u1, u_out = u2))
    @test qM ≈ LinmapQuant(xm, UnitMap(u_in = Vector(u1), u_out = Vector(u2)))
    @test qM ≈ LinmapQuant(dstrip(Matrix(qMraw)), dimension(Matrix(qMraw)))
    @test qM == LinmapQuant(qM)
    dm = dimension(qM)
    @test qM ≈ LinmapQuant(dstrip(qMraw), UnitMap(u_in = dm.u_in, u_out=dm.u_out.*dm.u_fac))

    qx = VectorQuant(dstrip(Vector(x)), dimension(Vector(x)))
    @test VectorQuant(ustrip.(x), u1) == VectorQuant(ustrip.(x).*u1)
    @test all(qx .≈ x)
    @test ubase(qx) == qx
    @test IndexStyle(typeof(VectorQuant(x))) == IndexStyle(typeof(x))

    #LU factorization
    luQ = lu(qM)
    luR = lu(dstrip.(qM))
    lup = MVector{3}(luQ.p)
    luM = lu(LinmapQuant(Matrix(qMraw)))
    @test luQ isa FactorQuant
    @test inv(luQ.factor) ≈ inv(luR)
    @test inv(luQ) ≈ inv(qM)
    @test LinearAlgebra.inv!(lu(LinmapQuant(Matrix(qMraw)))) ≈ inv(qM)    
    @test all( (luQ.L * luQ.U)[invperm(lup),:] .≈ qM )
    @test all( luQ.P * qM .≈ luQ.L*luQ.U )
    @test qM/luQ ≈ qM/qM
    @test luQ\qM ≈ qM\qM
    @test all(x .≈ luQ\y)
    @test all(x' .≈ y'/lu(LinmapQuant(Matrix(qMraw)')))
    @test ustrip(luM') == ustrip(luM)'
    @test dimension(luM') == dimension(luM)'
    @test ustrip(transpose(luM)) == transpose(ustrip(luM))
    @test dimension(transpose(luM)) == transpose(dimension(luM))
    #@test_throws "Unsupported method" dstrip(lu(qM, Val(true))) == lu(dstrip(qM), Val(true))
    #@test_throws "Unsupported method" dstrip(lu(qM, Val(false))) == lu(dstrip(qM), Val(false))
    
    #Inverses and transposes
    @test all(x .≈ inv(qM)*y)
    @test all(x .≈ FlexUnits.qinv(qMraw)*y)
    @test all(x' .≈ y'*inv(LinmapQuant(collect(qM'))))
    @test all(x' .≈ identity.(y)'*inv(LinmapQuant(collect(qM'))))
    @test all(Matrix(y') .≈ Matrix(x'*qM'))
    @test all(Matrix(transpose(y)) .≈ Matrix(transpose(x)*transpose(qM)))
    @test all(FlexUnits.qtranspose(x) .≈ x')
    @test all(transpose(FlexUnits.qtranspose(x)) .≈ x)
    @test all(FlexUnits.qadjoint(x) .≈ x')
    @test all((FlexUnits.qadjoint(x)') .≈ x)
    @test all(FlexUnits.qtranspose(qMraw) .≈ qM')
    @test all(FlexUnits.qtranspose(transpose(x)) .≈ x)
    @test all(transpose(FlexUnits.qtranspose(qMraw)) .≈ qM)
    @test all(FlexUnits.qadjoint(qMraw) .≈ qM')
    @test all((FlexUnits.qadjoint(qMraw)') .≈ qM)
    @test inv(qM')*LinmapQuant(collect(qM')) ≈ inv(qM')*qM'
    @test LinmapQuant(collect(qM'))*inv(qM') ≈ qM'*inv(qM')
    @test FlexUnits.qinv(qMraw) ≈ inv(qM)

    #Square matrices
    Σ = cov(randn(20,3)*rand(3,3))
    x = randn(3).*u2

    #Symmetric matrix
    rS = Σ.*inv.(u2).*inv.(u2)'
    qS = LinmapQuant(rS)
    @test all(qS .≈ ubase.(rS))
    @test x'*(rS)*x ≈ x'*qS*x
    @test FlexUnits.assert_symmetric(dimension(qS)) == dimension(qS)
    @test FlexUnits.assert_symmetric(dimension(qS')) == dimension(qS')

    #Repeatable matrix
    rR = Σ.*u2.*inv.(u2)'
    qR = LinmapQuant(rR)
    @test all(qR .≈ ubase.(rR))
    @test all((rR^2*x) .≈ (qR*qR*x))
    @test FlexUnits.assert_repeatable(dimension(qR)) == dimension(qR)
    @test FlexUnits.assert_repeatable(dimension(qR')) == dimension(qR')
    @test FlexUnits.assert_idempotent(dimension(qR)) == dimension(qR)
    @test FlexUnits.assert_idempotent(dimension(qR')) == dimension(qR')
    @test all(exp(qR) .≈ exp(dstrip(qR)) .* dimension(qR)) 
    @test all(log(qR) .≈ log(dstrip(qR)) .* dimension(qR))

    #Indexing
    @test UniformScaling(2ud"m/s")[2,2] == 2ud"m/s"
    @test FlexUnits.isunknown(dimension(UniformScaling(2ud"m/s")[1,2]))
    @test !FlexUnits.isunknown(dimension(UniformScaling(2u"m/s")[1,2]))

    @test qR[:,1] ≈ rR[:,1]
    @test qR[:,1] isa VectorQuant
    @test qR[2,:] ≈ rR[2,:]
    @test qR[2,:] isa VectorQuant
    @test qR[1:2, 1:2] ≈ rR[1:2, 1:2]
    @test qR[1:2, 1:2] isa LinmapQuant
    @test all(qR'[1:2, 1:3] .≈ rR'[1:2, 1:3])
    @test qR'[1,2] ≈ rR'[1,2]
    vR = qR[:,1]
    @test vR[1:2] ≈ rR[1:2, 1]
    @test vR[1:2] isa VectorQuant
    @test vR[1] ≈ rR[1,1]
    @test ustrip(vR) == dstrip(vR)
    @test unit(vR) == dimension(vR)
    @test ustrip(qR) == dstrip(qR)
    @test unit(qR) == dimension(qR)

    #Concatenation 
    A = SA[1.0 0.0; 1.0 1.0]*UnitMap(u_in=SA[u"m/s", u"m"], u_out=SA[u"m/s^2", u"m/s"])
    B = SA[0.2; 0]*UnitMap(u_in=SA[u"W"], u_out=SA[u"m/s^2", u"m/s"])
    C = SA[2.5 0.0; 0.0 1.0]*UnitMap(u_in=SA[u"m/s", u"m"], u_out=SA[u"N", u"m"])
    D = SA[0.0; 0.0]*UnitMap(u_in=SA[u"W"], u_out=SA[u"N", u"m"])

    MAB = [Matrix(A) Matrix(B)]
    AB = [A B]
    @test AB isa LinmapQuant
    @test AB ≈ MAB

    AB = [A Matrix(B)]
    @test AB isa LinmapQuant
    @test AB ≈ MAB

    AB = [A B[:]]
    @test AB isa LinmapQuant
    @test AB ≈ MAB

    BA = [B[:] A]
    @test BA isa LinmapQuant 
    @test BA ≈ [Matrix(B) Matrix(A)]


    MAC = [Matrix(A); Matrix(C)]
    AC = [A;C]
    @test AC isa LinmapQuant 
    @test AC ≈ MAC

    AC = [A;Matrix(C)]
    @test AC isa LinmapQuant
    @test AC ≈ MAC
    
    AC = [Matrix(A);C]
    @test AC isa LinmapQuant
    @test AC ≈ MAC

    MBB = [Matrix(B); Matrix(B)] 
    BB = [B; B[:]]
    @test BB isa LinmapQuant 
    @test BB ≈ MBB 

    BB = [B[:]; B]
    @test BB isa LinmapQuant 
    @test BB ≈ MBB 

    ABCD = [[A;C] [B;D]]
    @test ABCD isa LinmapQuant 
    @test [Matrix(A) Matrix(B); Matrix(C) Matrix(D)] ≈ ABCD

    AAB = [A A B]
    @test AAB isa LinmapQuant
    @test AAB ≈ [Matrix(A) Matrix(A) Matrix(B)]

    AAC = [A; A; C]
    @test AAC isa LinmapQuant
    @test AAC ≈ [Matrix(A); Matrix(A); Matrix(C)]

    @test Base.hcat(A) === A
    @test Base.vcat(A) === A
    @test FlexUnits.dm_hcat(dimension(A)) === dimension(A)
    @test FlexUnits.dm_vcat(dimension(A)) === dimension(A)
    @test Base.hcat([1 0; 1 0]*UnitMap(u_in=inv(u"kg"), u_out=u""), [1,0]*u"kPa") isa LinmapQuant

    #Concatenationg FlexQuant matrices/vectors 
    A = randn(3,3)
    B = randn(3)

    ABh = [A*u"m" B*u"kg"] 
    @test ABh isa LinmapQuant
    @test ABh ≈ [A.*u"m" B.*u"kg"]

    ABh_dims = [DimsMap(A*u"m") DimsMap(B*u"kg")]
    @test ABh_dims isa DimsMap
    @test ABh_dims == dimension(ABh)

    ABh_dims = FlexUnits.dm_hcat(dimension(A.*u"m"), dimension(B.*u"kg"))
    @test ABh_dims isa DimsMap
    @test ABh_dims == dimension(ABh)

    AAv = [A*u"m"; A*u"kg"]
    @test AAv isa LinmapQuant 
    @test AAv ≈ [A.*u"m"; A.*u"kg"]

    AAv_dims = [DimsMap(A*u"m"); DimsMap(A*u"kg")]
    @test AAv_dims isa DimsMap 
    @test AAv_dims == dimension(AAv)

    AAv_dims = FlexUnits.dm_vcat(dimension(A.*u"m"), dimension(A.*u"kg"))
    @test AAv_dims isa DimsMap 
    @test AAv_dims == dimension(AAv)

    BBh = [B*u"m" B*u"kg"] 
    @test BBh isa LinmapQuant 
    @test BBh ≈ [B.*u"m" B.*u"kg"]


    #Linear Regression Example
    XYRaw = randn(500,4)*rand(4,4)

    X = [XYRaw[:,1]*u"kg" XYRaw[:,2]*u"Pa" ones(500)*u""]
    @test X isa LinmapQuant

    Y = [XYRaw[:,3]*u"K" XYRaw[:,4]*u"W"]
    @test Y isa LinmapQuant

    B = (X'*X)\(X'*Y)
    @test B isa LinmapQuant 

    Yh = X*B
    @test dimension(Yh) == dimension(Y)

    #Matrix attributes
    @test adjoint(adjoint(qR)) == qR
    @test transpose(transpose(qR)) == qR
    @test adjoint(adjoint(vR)) == vR 
    @test transpose(transpose(qR)) == qR
    dR = dimension(qR)
    @test FlexUnits.dimtype(dR) == typeof(dimension(ud""))
    @test eltype(dR) == typeof(dimension(ud""))
    @test IndexStyle(typeof(dR)) == IndexCartesian()
    @test dR[4] == dR[CartesianIndices(dR)[4]]
    @test length(dR) == 9
    @test collect(dR) == dR[:,:]
    @test IndexStyle(typeof(dimension(rR))) == IndexStyle(rR)
    @test IndexStyle(typeof(dstrip(rR))) == IndexStyle(rR)
    @test DimsMap(qR) == dimension(qR)

    #Nonlinear mapping
    pumpunits = UnitMap(PumpInput(current=u"A", voltage=u"V"), PumpOutput(power=u"W", pressure=u"Pa", flow=u"m^3/s"))
    upumpfunc = FunctionQuant(pumpfunc, pumpunits)
    qinput = PumpInput(current=500*u"mA", voltage=6u"V")
    @test all(upumpfunc(qinput) .≈ pumpfunc(ustrip.(uinput(pumpunits), qinput)).*uoutput(pumpunits))

    #Matrix and vector operations 
    m  = LinmapQuant(SA[1.0 0.1; 0.2 1.0], UnitMap(u_in = SA[u"kg/s", u"kW"], u_out=SA[u"m^3/s", u"kPa"]))
    mi = inv(m)
    x  = VectorQuant(SA[0.5u"kg/s", 0.5u"kW"])
    y  = m*x

    mraw  = SMatrix{2,2}(m)
    miraw = SMatrix{2,2}(mi)
    xraw  = SVector{2}(x)
    yraw  = SVector{2}(y)
     
    @test x + x ≈ xraw + xraw
    @test x + xraw ≈ xraw + x
    @test m + m ≈ mraw + mraw
    @test m + mraw ≈ mraw + m

    @test x - x ≈ xraw - xraw
    @test x - xraw ≈ xraw - x
    @test m - m ≈ mraw - mraw
    @test mraw - mraw ≈ m - m
    @test m - mraw ≈ mraw - m

    @test m*mi ≈ mraw*miraw
    @test m*miraw ≈ mraw*mi

    @test m\m  ≈ miraw*mraw
    @test m\mraw ≈ miraw*m
    @test m'\m' ≈ miraw'*mraw'
    @test m'\mraw' ≈ miraw'*m'
    @test mraw\m ≈ miraw*m

    @test all(m*(2u"m/s") .≈ mraw.*(2u"m/s"))
    @test all((2u"m/s")*m .≈ (2u"m/s").*mraw)
    @test m*(2u"m/s") isa LinmapQuant
    @test (2u"m/s")*m isa LinmapQuant

    @test m/m ≈ mraw*miraw
    @test mraw/m ≈ m*miraw
    @test m/mraw ≈ m*miraw
    @test m'/m' ≈ mraw'*miraw'
    @test mraw'/m' ≈ m'*miraw'

    @test m*x ≈ yraw
    @test mraw*x ≈ y
    @test all(x'*mraw' .≈ y')
    @test mi*y ≈ xraw
    @test (x'*m')' ≈ yraw 
    @test (y'*mi')' ≈ xraw
    @test m\y ≈ xraw
    @test (y'/m')' ≈ xraw
    @test (y'/mraw')' ≈ xraw
    @test (yraw'/LinmapQuant(mraw'))' ≈ xraw
    @test m\yraw ≈ xraw
    @test xraw ≈ mraw\y 
    @test transpose(transpose(y)) ≈ yraw
    @test all(Diagonal(xraw)\x .≈ xraw./xraw)

    mrep = LinmapQuant(SA[1.0 0.1; 0.2 1.0], UnitMap(u_in = SA[u"kg/s", u"kW"], u_out=SA[u"kg/s", u"kW"]))
    mraw = SMatrix{2,2}(mrep)

    @test mrep^2 ≈ mraw*mraw
    @test mrep^3 ≈ mrep*mrep*mrep
    @test mrep^2.0 ≈ mraw*mraw
    @test (mrep')^2 ≈ mraw'*mraw'
    @test dimension(exp(mrep)) == dimension(mrep)
    @test mraw^2 ≈ mraw*mraw
    @test mraw^2.0 ≈ mraw*mraw

    @test dstrip(exp(mrep)) ≈ exp(dstrip(mrep))
    @test dstrip(exp(mraw)) ≈ exp(dstrip(mrep))
    @test dstrip(log(mraw)) ≈ log(dstrip(mrep))

    @test FlexUnits.assert_allequal(identity, [dimension(ud"m"), dimension(ud"m")]) isa Dimensions
    @test FlexUnits.assert_allequal(identity, [dimension(u"m"), dimension(u"m")]) isa StaticDims
    @test FlexUnits.assert_allequal([dimension(ud"m"), dimension(ud"m")]) isa Dimensions
    @test FlexUnits.assert_allequal([dimension(u"m"), dimension(u"m")]) isa StaticDims


end
