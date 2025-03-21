using DimensionalUnits
using TestItems: @testitem
using TestItemRunner

#Run these commands at startup to see coverage
#julia --startup-file=no --depwarn=yes --threads=auto -e 'using Coverage; clean_folder("src"); clean_folder("test")'
#julia --startup-file=no --depwarn=yes --threads=auto --code-coverage=user --project=. -e 'using Pkg; Pkg.test(coverage=true)'
#julia --startup-file=no --depwarn=yes --threads=auto coverage.jl

#To see the actual coverage in VSCode, install the Coverage Gutters extension
#https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters

@testitem "Assorted unittests" begin
    include("unittests.jl")
end
