using BEM
using Documenter

DocMeta.setdocmeta!(BEM, :DocTestSetup, :(using BEM); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [BEM],
    authors = "Lucas",
    repo = "https://github.com/l-s-campos/BEM.jl/blob/{commit}{path}#{line}",
    sitename = "BEM.jl",
    format = Documenter.HTML(; canonical = "https://l-s-campos.github.io/BEM.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/l-s-campos/BEM.jl")
