@testset "covariances under gaussian with coalescent" begin

v0=0.1

# 1-taxon network with 1 hybrid (2 parallel edges)
ell = 2.0; # in coalescent units
net = readnewick("((t:0.0)#H1:$ell::0.6,#H1:$ell)r;")
cM, eV = PhyloTraits.gaussiancoalescent_covariancematrix(net,v0/1) # λ = v0/σ2eq
#= theoretically:
v = ell + v0 = 2.1 # variance at the tip t
q = 1-exp(-ell) = 0.8646647167633873 # prob(coalescence) along each hybrid edge
r = 1-q/ell = 0.5676676416183064 # proportion of shared coalescence time
c = (0.6^2 +0.4^2)*(q*v0 + ell*r) = 0.6353369125547349 # covariance at tip t
=#
c = 0.6353369125547349
v1 = 2.1 - c # expected within-species variance = total variance - cov(2 ind.)
@test cM[:all] ≈ [0 0 0; 0 c c; 0 c c]
@test eV[:all] ≈ [v0, v1, v1]

# 5 taxa, level-2, subnet on A,B,C is a tree, d1 below 1 hybrid, d2 below 2
net = readnewick("(((C:2,#H1:.1):0.3,(((d1:1,#H2:.1):.8)#H1:.7::.6,(d2:.5)#H2:1::.7):.4):.3,(B:1,A:.5):2);");
# plot(net, showedgelength=true, showgamma=true);
cM, eV = PhyloTraits.gaussiancoalescent_covariancematrix(net,v0/1)
@test tiplabels(eV) == tiplabels(cM) == ["C","d1","d2","B","A"]
@test eV[:tips] + diag(cM[:tips]) ≈ [2.7,3.02,2.396,3.1,2.6]
v11=1.7668462203929007; v12=0.11761402816197718; v13=0.08199968747807541
v22=2.038420771136349; v23=0.443417290678788; v33=1.247326232249089
v44=2.144808361531078; v45=1.2218017549129516; v55=1.6738764987615091
@test cM[:tips] ≈ [v11 v12 v13 0 0; v12 v22 v23 0 0; v13 v23 v33 0 0;
    0 0 0 v44 v45; 0 0 0 v45 v55]
end

@testset "fit σ2 under the gaussian with coalescent" begin

m = PhyloTraits.GaussianCoalescent(0.1,0.05) # v0=0.1, σ2=0.05
@test PhyloTraits.lambda(m) ≈ 2
PhyloTraits.lambda!(m, 4)
@test m.v0 ≈ 0.2
PhyloTraits.setNe!(m, 10)
s = IOBuffer(); show(s, m)
@test String(take!(s)) == """Gaussian with coalescent:
v0: 0.2
σ2: 0.005
Ne: 10
σ2 equilibrium (within pop): 0.05
λ: 4
"""

nwk = "(A:2.5,((B:1,#H1:0.5::0.1):1,(C:1,(D:0.5)#H1:0.5::0.9):1):0.5);"
net = readnewick(nwk)
Y = [8.60,10.56,11.9,9.96,11.24,]
X = ones(5, 1)
df = DataFrame(trait = Y, tipnames = ["B","C","A","D","A",])
f0 = phylolm(@formula(trait ~ 1),df,net; reml=true, model="gaussiancoalescent",
        paramlist=Dict(:lambda => (start=0.1, fixed=true)))
@test f0.evomodel.Ne == 1
@test lambda_estim(f0) == 0.1
@test sigma2_phylo(f0) ≈ 0.5470357584927545
@test f0.evomodel.v0 ≈ 0.05470357584927545
@test loglikelihood(f0) ≈ -6.704652757076957
@test dof(f0) == 2 # intercept, s2, but not lambda (fixed)
@test aic(f0) ≈ 17.409305514153914

f0 = phylolm(@formula(trait ~ 1),df,net; reml=false, model="gaussiancoalescent",
        paramlist=Dict(:lambda => (start=0.1,)))
@test f0.evomodel.Ne == 1
@test lambda_estim(f0) ≈ 13.525028685269309
@test sigma2_phylo(f0) ≈ 0.08974694768102379
@test f0.evomodel.v0 ≈ 1.2138300418012107
@test loglikelihood(f0) ≈ -6.907945500519682
@test dof(f0) == 3 # intercept, s2, lambda
@test aic(f0) ≈ 19.815891001039365

# bounds
f0 = phylolm(@formula(trait ~ 1),df,net; reml=false, model="gaussiancoalescent",
        paramlist=Dict(:lambda => (start=0.1, upper=11.3)))
@test lambda_estim(f0) ≈ 11.3
f0 = phylolm(@formula(trait ~ 1),df,net; reml=false, model="gaussiancoalescent",
        paramlist=Dict(:lambda => (start=18.4, lower=15.3)))
@test lambda_estim(f0) ≈ 15.3

f0 = phylolm(@formula(trait ~ 1),df,net; reml=false, model="gaussiancoalescent",
        paramlist=Dict(:lambda => (start=0, fixed=true)))
@test f0.evomodel.Ne == 1
@test lambda_estim(f0) == 0.0
@test sigma2_phylo(f0) ≈ 0.4539365904733537
@test f0.evomodel.v0 == 0.0
@test loglikelihood(f0) ≈ -7.116223253873285
@test dof(f0) == 2

for e in net.edge e.length *= 2; end
f0 = phylolm(@formula(trait ~ 1),df,net; reml=true, model="gaussiancoalescent",
        paramlist=Dict(:lambda => (start=0.1, fixed=true), :Ne => (start=2,)))
tmp = """Model: Gaussian with coalescent

Parameter Estimates, using REML:
phylogenetic variance rate: 0.547036
v0: 0.0547036
σ2: 0.273518
Ne: 2
σ2 equilibrium (within pop): 0.547036
λ: 0.1
"""
s = repr("text/plain", f0)
@test occursin(tmp, s)
@test dof(f0) == 2

@test_throws "estimate either λ or v0" phylolm(@formula(trait~1),df,net;
  model="gaussiancoalescent", paramlist=Dict(
    :lambda => (start=1,fixed=true), :v0 => (start=0.1,fixed=true)))
@test_throws "optimization of λ under a fixed v0" phylolm(@formula(trait~1),df,net;
  model="gaussiancoalescent", paramlist=Dict(
    :lambda => (start=1,), :v0 => (start=0.1,fixed=true)))


## On a tree : comparison with phylolm in R
nwk = "(A:2.5,(B:2,(C:1.5,D:1.5):0.5):0.5);"
net = readnewick(nwk)
Y = [8.60,10.56,11.9,9.96,]
X = ones(4, 1)
df = DataFrame(trait = Y, tipnames = ["B","C","A","D",])

f0 = phylolm(@formula(trait ~ 1),df,net; reml=false, model="gaussiancoalescent",
        paramlist=Dict(:lambda => (start=0.1, fixed=true)))
@test f0.evomodel.Ne == 1
@test lambda_estim(f0) ≈ 0.1
@test sigma2_phylo(f0) ≈ 0.5452923 atol=1e-6
@test f0.evomodel.v0 ≈ 0.05452923 atol=1e-6
@test coef(f0) ≈ [10.29722] atol=1e-6
@test loglikelihood(f0) ≈ -6.357265 atol=1e-6

## R code to reproduce the results
# remotes::install_github("pbastide/phylolm", ref = "ils")
# library(phylolm)
# tree <- read.tree(text = "(A:2.5,(B:2,(C:1.5,D:1.5):0.5):0.5);")
# trait <- c(8.60,10.56,11.9,9.96)
# names(trait) <- c("B","C","A","D")
# 
# fit <- phylolm(trait ~ 1, phy = tree, model = "ILS",
#                lower.bound = list(lambda_ILS = 0.1),
#                upper.bound = list(lambda_ILS = 0.1),
#                starting.value = list(lambda_ILS = 0.1))
# fit$sigma2
# fit$coefficients
# fit$logLik
end
