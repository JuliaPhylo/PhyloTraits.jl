# PhyloTraits: trait evolution along phylogenies

|**Documentation**| **Build Status** | **Code Coverage**   | **Style Guide** |
|:---------------:|:----------------:|:-------------------:|:----------------|
|[![stable][docs-stable-img]][docs-stable-url] [![dev][docs-dev-img]][docs-dev-url] | [![build][build-img]][build-url] [![PkgEval][pgkeval-img]][pgkeval-url] [![aqua][aqua-img]][aqua-url] | [![coverage][codecov-img]][codecov-url] | [![Code Style: Blue][style-img]][style-url] [![collaborative][colprac-img]][colprac-url]

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaPhylo.github.io/PhyloTraits.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaPhylo.github.io/PhyloTraits.jl/dev/

[build-img]: https://github.com/JuliaPhylo/PhyloTraits.jl/actions/workflows/CI.yml/badge.svg?branch=main
[build-url]: https://github.com/JuliaPhylo/PhyloTraits.jl/actions/workflows/CI.yml?query=branch%3Amain
[pgkeval-img]: https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/P/PhyloTraits.svg
[pgkeval-url]: https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/P/PhyloTraits.html
[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[codecov-img]: https://codecov.io/gh/JuliaPhylo/PhyloTraits.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaPhylo/PhyloTraits.jl

[style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[style-url]: https://github.com/invenia/BlueStyle
<!-- ColPrac: Contributor's Guide on Collaborative Practices for Community Packages -->
[colprac-img]: https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet
[colprac-url]: https://github.com/SciML/ColPrac

[Julia](http://julialang.org/) package for the analysis of trait evolution along
a phylogeny, including phylogenetic networks. It depends on utilities from
[PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl),
from the [JuliaPhylo](https://github.com/JuliaPhylo) project.

## Quick start

To install & use the package, install [Julia](https://julialang.org/downloads/)
then type, within a Julia session:
```julia-repl
julia> using PhyloTraits
```
The first time, julia will prompt you for package installation.

## Tutorials

- Follow PhyloTraits's [latest documentation][docs-dev-url].
  It includes a manual and examples of empirical analyses from the literature.
- This [tutorial](https://cecileane.github.io/networkPCM-workshop/)
  (from a 2023 workshop) is on comparative methods
  and includes a tutorial for network calibration.
  Caveat: it does *not* use PhyloTraits. It uses older versions of the
  functions from PhyloNetworks v0.16 (see Note below).
  Some function names and options will need to be adjusted.

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
Below is a short list.

For continuous traits, analyses based on the Brownian motion process,
with or without transgressive evolution after reticulations:

- Bastide, Solís-Lemus, Kriebel, Sparks & Ané (2018).
  Phylogenetic Comparative Methods for Phylogenetic Networks with Reticulations.
  Systematic Biology, 67(5):800–820.
  [doi:10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).
  SI on [dryad](http://dx.doi.org/10.5061/dryad.nt2g6)
  including a tutorial for trait evolution
  and a tutorial for network calibration.

Continuous traits, accounting for within-species variation:

- Teo, Rose, Bastide & Ané (2023).
  Accounting for intraspecific variation in continuous trait evolution
  on a reticulate phylogeny.
  Bulletin of the Society of Systematic Biologists, 2(3):1-29.
  doi: [10.18061/bssb.v2i3.8977](https://doi.org/10.18061/bssb.v2i3.8977)

For a discrete trait (influence of gene flow on the trait,
ancestral state reconstruction, rates):

- Karimi, Grover, Gallagher, Wendel, Ané & Baum (2020). Reticulate evolution
  helps explain apparent homoplasy in floral biology and pollination in baobabs
  (*Adansonia*; Bombacoideae; Malvaceae).
  Systematic Biology,
  69(3):462-478. doi: [10.1093/sysbio/syz073](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syz073/5613901?guestAccessKey=a32e7dd3-27fd-4a13-b171-7ff5d6da0e01).

> [!NOTE]
> Much of this package was formerly part of PhyloNetworks v0.16.4 (and prior).
> PhyloNetworks v0.17, 1.0 (and later) are stripped of functions for trait evolution.
