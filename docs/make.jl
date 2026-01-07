#julia --project make.jl

using Documenter, FlexUnits

DocMeta.setdocmeta!(FlexUnits, :DocTestSetup, :(using FlexUnits, .UnitRegistry))
makedocs(
    sitename="FlexUnits Documentation",
    pages = [
            "Home" => "index.md",
            "Types" => "types.md"
    ]
)
