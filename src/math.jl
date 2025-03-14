"""
    map_dimensions(f::F, args::AbstractDimensions...)

Similar to the `map` function, but specifically iterates over dimensions
"""
function map_dimensions(f::F, args::AbstractDimensions...) where {F<:Function}
    D = promote_type(typeof(args).parameters...)
    return  D(
        ( f((getproperty(arg, dim) for arg in args)...) for dim in dimension_names(D) )...
    )
end


