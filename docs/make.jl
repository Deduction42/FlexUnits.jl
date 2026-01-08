#julia --project make.jl

using Documenter, FlexUnits

DocMeta.setdocmeta!(FlexUnits, :DocTestSetup, :(using FlexUnits, .UnitRegistry))
makedocs(
    sitename="FlexUnits.jl",
    pages = [
        "Home" => "index.md",
        "Performance" => "performance.md",
        "Types" => "types.md"
    ]
)
