# notes to maintain documentation

- built with [Documenter](https://documenter.juliadocs.org/stable/).
- deployed [here](https://juliaphylo.github.io/PhyloTraits.jl/)
  using GitHub and files committed to the `gh-pages` branch:
  go to `dev/`, or manually type the url to preview the docs from a pull request:
  [https://juliaphylo.github.io/PhyloTraits.jl/previews/PR#/] with `#` replaced
  by the pull request number.

## quick start: create a local version of the website

generally: to see deprecation warnings, add the option `--depwarn=yes`.
go into the `docs/` folder, load the project environment, then run `make.jl`:

```shell
cd docs
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=dirname(pwd())))'
julia --project --color=yes make.jl
```

or interactively in `docs/`, if needing to have finer control of loaded dependencies:

```shell
pkg> activate .
pkg> status # just to check
pkg> status --manifest
pkg> instantiate # after deleting Manifest.toml
# do what you need to do to update dependencies / dev packages
julia> include("make.jl")
```

or, after project & manifest setup, from the main repo instead of `docs/`:
```shell
julia --project=docs/ --color=yes docs/make.jl
```

it will:
- test the `jldoctest` blocks of examples in the docstrings
- create or update a `build/` directory with html files.

To see the result, just open `docs/build/index.html` in a browser and follow the links.

## how it works: overview

- The `docs` job (named "Documentation") of `.github/workflows/CI.yml`
  * installs R, installs Julia
  * runs Julia commands that asks to start the doc project, which will trigger
    the installation of dependencies listed in `docs/Project.toml`,
  * then run `doctest(PhyloTraits)`, which will call `./docs/make.jl`.
- the julia script `docs/make.jl` includes these steps:
  1. run `makedocs()` to create the documentation website.
     also runs all `jldoctest` blocks in the source files, to check that
     the output in the blocks matches the actual output.
     This steps also translate the "Documenter md" documentation files
     into html files.
  2. run `deploydocs(...)` also from Documenter:
     to push the files on github, gh-pages branch.

for now, docstrings are automatically used to build an entry for
- each internal thing that has a docstring (e.g. not exported in `src/PhyloTraits.jl`)
- each public thing.

## The "Documenter md" format

The documentation pages are all written in this format. It is a regular md, but
with extra blocks of codes (as `@example` and `@setup`) that contain Julia
commands. These lines will be executed during the `makedoc()` process. See the
`Documenter` [doc](https://juliadocs.github.io/Documenter.jl/stable/man/syntax/)
for more details on the syntax. For instance, @example blocks with the same "name"
are run in the same session. Otherwise, an @example blocks with no name
is run in its own anonymous Module.

## documentation figures & plots

Some of these blocks may contain plots, which are going to be drawn during the
process, requiring the use of `PhyloPlots` along with `RCall`. Hence,
before the doc is built, `.github/workflows/ci.yml` installs `R` on the server,
sets up the julia environment with dependencies like `PhyloPlots` before
starting the build in itself.

### directory to save figures

We chose to group all the output plots in the directory `assets/figures`.
Hence, the typical setup in a documentation page containing plots is:

    ```@setup name
    using PhyloPlots, RCall
    mkpath("../assets/figures")
    figname(x) = joinpath("..", "assets", "figures", x) # hide
    ```

The `mkpath` command is there to ensure that the target directory does indeed
exist. In theory, it only needs to be called once (by the first documentation
page being built). However, as this order might be subject to change over time,
it could be safer to include it on every such page.

### figure format

After trial and error and discussions, we chose to use only the *SVG* format
for plots. This format should ensure that when a plot is drawn again,
identical, in a new build, then Git will recognize that it has not changed, and
hence not save a new version of it. This should ensure that the repository does
not grow unreasonably at each build of the doc, i.e. after each push to
master. The typical commands to save and display a plot should hence be:

    ```@example name
    R"svg"(figname("my_useful_name.svg"), width=4, height=4) # hide
    plot(net);
    R"dev.off()" # hide
    nothing # hide
    ```
    ![my_useful_name](../assets/figures/my_useful_name.svg)

**Warning**: this is not like an interactive session. If the same file name
is re-used by some other documentation page for some other plot, only the
final version of the plot will be committed by git, with possible unintended
consequences. Make sure to use different file names for plots that are supposed
to look different (across the whole site).

## references

big difference:

    [blabla](@ref)
    [`blabla`](@ref)

The first version will look for a *section* header "blabla", to link to that section.
The secon version will look for a *function* named "blabla",
to link to the documentation for that function.

### external references

References to functions from external packages are handled with
[`DocumenterInterLinks.jl`](https://juliadocs.org/DocumenterInterLinks.jl/stable/).
To set up external references, `make.jl` needs to be modified to include
the package for referencing. For example, to make references to functions in
`PhyloNetworks`, `make.jl` has these additional lines of code:

```julia
# Interlink with PhyloNetworks
using DocumenterInterLinks
links = InterLinks(
    "PhyloNetworks" => "https://juliaphylo.github.io/PhyloNetworks.jl/stable/objects.inv"
)
```

The `links` object is then used in the `plugins` argument of `makedocs()` ran
at the end of `make.jl`.
Now `@ref` can be used to reference things in `PhyloNetworks`!
We can also use `@extref` to specify that a reference is an external one.
