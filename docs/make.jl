using Documenter, ImageUtils

makedocs(
    modules = [ImageUtils],
    format = :html,
    checkdocs = :exports,
    sitename = "ImageUtils.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/tknopp/ImageUtils.jl.git",
)
