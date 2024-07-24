using PhyloTraits
using Documenter

DocMeta.setdocmeta!(PhyloTraits, :DocTestSetup, :(using PhyloTraits); recursive=true)

makedocs(;
    modules=[PhyloTraits],
    authors="Cecile Ane <cecileane@users.noreply.github.com>, Paul Bastide <pbastide@users.noreply.github.com>, and contributors",
    sitename="PhyloTraits.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaPhylo.github.io/PhyloTraits.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaPhylo/PhyloTraits.jl",
    devbranch="main",
)
