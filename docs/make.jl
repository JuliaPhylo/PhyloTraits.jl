using PhyloTraits
using Documenter

# NOTE: this installs the dev versions of PhyloNetworks and PhyloPlots for compatibility.
# To be edited when update to non dev versions.
using Pkg
Pkg.add(PackageSpec(name="PhyloNetworks", rev="dev"))
Pkg.add(PackageSpec(name="PhyloPlots", rev="dev11"))

# Interlink with PhyloNetworks
using DocumenterInterLinks
links = InterLinks(
    "PhyloNetworks" => "https://juliaphylo.github.io/PhyloNetworks.jl/stable/objects.inv"
)

# NOTE: default loading of PhyloNetworks in all docstring examples
DocMeta.setdocmeta!(PhyloTraits, :DocTestSetup, :(using PhyloNetworks, PhyloTraits); recursive=true)


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
            "Continuous Trait Evolution" => "man/trait_tree.md",
            "Discrete Trait Evolution" => "man/fitDiscrete.md",
        ]
    ],
    plugins=[links],
)

deploydocs(;
    repo="github.com/JuliaPhylo/PhyloTraits.jl",
    devbranch="main",
)
