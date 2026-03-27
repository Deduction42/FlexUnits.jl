module StatisticsExt 
    import Statistics 
    import FlexUnits
    import FlexUnits: LinmapQuant, Quantity, dstrip, dimension, ureduce
    import FlexUnits: DimsMap, uinput, uoutput, ufactor

    #===========================================================================================================================================
    Main API extension
    ===========================================================================================================================================#
    Statistics.mean(f, q::LinmapQuant; dims=:) = LinmapQuant(
        Statistics.mean(f, dstrip(q), dims=dims), 
        ureduce(f, dimension(q), dims=dims)
    )
    Statistics.mean(q::LinmapQuant; dims=:) = LinmapQuant(
        Statistics.mean(dstrip(q), dims=dims), 
        ureduce(dimension(q), dims=dims)
    )
    Statistics.median(q::LinmapQuant; dims=:) = LinmapQuant(
        Statistics.median(dstrip(q), dims=dims),
        ureduce(dimension(q), dims=dims)
    )

    Statistics.std(q::LinmapQuant; corrected=true, mean=nothing, dims=:) = _with_mean_kwarg(Statistics.std, q; corrected=corrected, mean=mean, dims=dims)
    Statistics.stdm(q::LinmapQuant, mean; corrected=true, dims=:) = _with_mean_arg(Statistics.stdm, q, mean, corrected=corrected, dims=dims)

    Statistics.var(q::LinmapQuant; corrected=true, mean=nothing, dims=:) = _usquared(_with_mean_kwarg(Statistics.var, q; corrected=corrected, mean=mean, dims=dims))
    Statistics.varm(q::LinmapQuant, mean; corrected=true, dims=:) = _usquared(_with_mean_arg(Statistics.varm, q, mean; corrected=corrected, dims=dims))

    function Statistics.cov(q::LinmapQuant; corrected=true, dims=1)
        d = _dimcov(q, dims=dims)
        x = Statistics.cov(dstrip(q), corrected=corrected, dims=dims)
        return LinmapQuant(x, d)
    end

    function Statistics.cov(q1::LinmapQuant, q2::LinmapQuant; corrected=true, dims=1)
        d = _dimcov(q1, q2, dims=dims)
        x = Statistics.cov(dstrip(q1), dstrip(q2), corrected=corrected, dims=dims)
        return LinmapQuant(x, d)
    end 

    #Correlation matrices are always dimensionless, so simply dstrip them
    Statistics.cor(q::LinmapQuant; corrected=true, dims=1) = Statistics.cor(dstrip(q), corrected=corrected, dims=dims)
    Statistics.cor(q1::LinmapQuant, q2::LinmapQuant; corrected=true, dims=1) = Statistics.cor(dstrip(q1), dstrip(q2), corrected=corrected, dims=dims)


    #===========================================================================================================================================
    Helper functions
    ===========================================================================================================================================#
    #Conditionally applies a vector stat with or without mean
    function _with_mean_kwarg(f, q::LinmapQuant; mean=nothing, corrected=true, dims=:)
        d = ureduce(dimension(q), dims=dims)
        x = if isnothing(mean)
            f(dstrip(q), corrected=corrected, dims=dims)
        else
            f(dstrip(q), mean=_dstrip(d, mean), corrected=corrected, dims=dims)
        end
        return LinmapQuant(x, d)
    end

    function _with_mean_arg(f, q::LinmapQuant, mean; corrected=true, dims=:)
        d = ureduce(dimension(q), dims=dims)
        x = f(dstrip(q), _dstrip(d, mean), corrected=corrected, dims=dims)
        return LinmapQuant(x, d)
    end 

    #Covariancce of a DimsMap
    function _dimcov(d1::DimsMap, d2::DimsMap; dims=1)
        if dims == 1
            return d1'*d2 
        elseif dims == 2 
            return d1*d2'
        else
            throw(ArgumentError("Dimension argument for LinmapQuant can only be 1 or 2"))
        end
    end
    _dimcov(d::DimsMap; dims=1) = _dimcov(d, d, dims=dims)
    _dimcov(q1::LinmapQuant, q2::LinmapQuant; dims=1) = _dimcov(dimension(q1), dimension(q2), dims=dims)    
    _dimcov(q::LinmapQuant; dims=1) = _dimcov(dimension(q), dims=dims)

    #Produces a DimsMap that contains squared elements
    _usquared(d::DimsMap) = DimsMap(u_fac=abs2(ufactor(d)), u_in=abs2.(uinput(d)), u_out=abs2.(uoutput(d)))
    _usquared(q::LinmapQuant) = LinmapQuant(dstrip(q), _usquared(dimension(q)))
    
    #Checks dimension of vector to see if it matches the ureduce results
    function _dstrip(d::DimsMap, v::AbstractMatrix{<:Quantity})
        size(v) == size(d) || throw(DimensionMismatch("Quantity matrix had size $(size(v)) but dimensions had size $(size(d)), sizes must be the same"))
        return dstrip.(d, v)
    end

end