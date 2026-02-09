#julia --project make.jl

using Documenter, FlexUnits
import FlexUnits.QuantUnion

DocMeta.setdocmeta!(FlexUnits, :DocTestSetup, :(using FlexUnits, .UnitRegistry))
makedocs(
    sitename="FlexUnits.jl",
    pages = [
        "Home" => "index.md",
        "Performance" => "performance.md",
        "Unit Manipulation" => "manipulation.md",
        "Linear Algebra" => "linearalgebra.md",
        "Advanced Examples" => "examples.md",
        "Types" => "types.md"
    ]
)
