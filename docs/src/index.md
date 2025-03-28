```@meta
CurrentModule = PhyloTraits
```

# PhyloTraits

[PhyloTraits](https://github.com/JuliaPhylo/PhyloTraits.jl)
is a [Julia](http://julialang.org) package with core utilities for
phylogenetic networks.
See the [PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl)
package, which PhyloTraits depends on, for background on phylogenetic networks
and for how to get help.

---

## References

See them in [bibtex format](https://github.com/juliaphylo/PhyloTraits.jl/blob/master/CITATION.bib).

Quick links for methods about trait evolution on networks:
- Teo, Rose, Bastide & Ané (2023).
  Accounting for intraspecific variation in continuous trait evolution
  on a reticulate phylogeny.
  Bulletin of the Society of Systematic Biologists, 2(3):1-29.
  doi: [10.18061/bssb.v2i3.8977](https://doi.org/10.18061/bssb.v2i3.8977)
- Karimi, Grover, Gallagher, Wendel, Ané & Baum (2020). Reticulate evolution
  helps explain apparent homoplasy in floral biology and pollination in baobabs
  (*Adansonia*; Bombacoideae; Malvaceae).
  Systematic Biology, 69(3):462-478.
  [doi:10.1093/sysbio/syz073](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syz073/5613901?guestAccessKey=a32e7dd3-27fd-4a13-b171-7ff5d6da0e01).
- Bastide, Solís-Lemus, Kriebel, Sparks, Ané (2018).
  Phylogenetic Comparative Methods for Phylogenetic Networks with Reticulations.
  Systematic Biology, 67(5):800–820.
  [doi:10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).

## Manual

The manual pages contain detailed tutorials on how to use the functions of the package.

```@contents
Pages = [
    "man/phyloregression.md",
    "man/simulate_continuous.md",
    "man/fitDiscrete.md",
    "man/fitdiscreteDNA.md",
    "man/simulate_discrete.md",
]
Depth = 3
```

## Examples

This section contains example of empirical analyses from the literature.

```@contents
Pages = [
    "man/example_fishes.md",
]
Depth = 3
```

## Library

For help on individual functions, see this library.

```@contents
Pages = [
    "lib/public.md",
    "lib/internal.md",
]
Depth = 3
```
