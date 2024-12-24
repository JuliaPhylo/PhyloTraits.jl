using PhyloTraits
using Documenter

# Interlink with PhyloNetworks
using DocumenterInterLinks
links = InterLinks(
    "PhyloNetworks" => "https://juliaphylo.github.io/PhyloNetworks.jl/stable/objects.inv"
)

# NOTE: default loading of PhyloNetworks in all docstring examples
DocMeta.setdocmeta!(
    PhyloTraits,
    :DocTestSetup,
    :(using PhyloNetworks; using PhyloTraits);
    recursive=true)

makedocs(;
    modules=[PhyloTraits],
    authors="Cecile Ane <cecileane@users.noreply.github.com>, Paul Bastide <pbastide@users.noreply.github.com>, and contributors",
    sitename="PhyloTraits.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true", # easier local build
        size_threshold = 600 * 2^10,
        size_threshold_warn = 500 * 2^10, # 600 KiB
        canonical="https://JuliaPhylo.github.io/PhyloTraits.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Continuous trait analysis" => "man/phyloregression.md",
            "Continuous trait simulation" => "man/simulate_continuous.md",
            "Discrete trait analysis" => "man/fitDiscrete.md",
            "DNA evolutionary models" => "man/fitdiscreteDNA.md",
            "Discrete trait simulation" => "man/simulate_discrete.md",
        ],
        "Library" => [
            "public" => "lib/public.md",
            "internals" => "lib/internal.md",
        ]
    ],
    plugins=[links],
)

deploydocs(;
    repo="github.com/JuliaPhylo/PhyloTraits.jl",
    devbranch="main",
)
