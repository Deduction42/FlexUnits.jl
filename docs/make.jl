#julia --project make.jl

using Documenter, FlexUnits

DocMeta.setdocmeta!(FlexUnits, :DocTestSetup, :(using FlexUnits, .UnitRegistry))
makedocs(
    sitename="FlexUnits.jl",
    pages = [
        "Home" => "index.md",
        "Performance" => "performance.md",
        "Linear Algebra" => "linearalgebra.md",
        "Advanced Examples" => "examples.md",
        "Types" => "types.md"
    ]
)
