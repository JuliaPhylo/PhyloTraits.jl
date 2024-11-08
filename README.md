# PhyloTraits: trait evolution along phylogenies

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaPhylo.github.io/PhyloTraits.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaPhylo.github.io/PhyloTraits.jl/dev/)
[![Build Status](https://github.com/JuliaPhylo/PhyloTraits.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaPhylo/PhyloTraits.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaPhylo/PhyloTraits.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaPhylo/PhyloTraits.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/P/PhyloTraits.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/P/PhyloTraits.html)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[Julia](http://julialang.org/) package for the analysis of trait evolution along
a phylogeny, including phylogenetic networks. It depends on utilities from
[PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl).

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
Below is a short list.

For continuous traits, analyses based on the Brownian motion process,
with or without transgressive evolution after reticulations:

- Bastide, Solís-Lemus, Kriebel, Sparks, Ané (2018).
  Phylogenetic Comparative Methods for Phylogenetic Networks with Reticulations.
  Systematic Biology, 67(5):800–820.
  [doi:10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).
  SI on [dryad](http://dx.doi.org/10.5061/dryad.nt2g6)
  including a tutorial for trait evolution
  and a tutorial for network calibration.

Continuous traits, accounting for within-species variation:

- Benjamin Teo, Jeffrey P. Rose, Paul Bastide & Cécile Ané (2022).
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
> PhyloNetworks v0.17 will be stripped of functions for trait evolution.
