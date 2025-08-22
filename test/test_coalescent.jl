@testset "covariances under gaussian with coalescent" begin

v0=0.1

# 1-taxon network with 1 hybrid (2 parallel edges)
ell = 2.0; # in coalescent units
net = readnewick("((t:0.0)#H1:$ell::0.6,#H1:$ell)r;")
cM, eV = PhyloTraits.gaussiancoalescent_covariancematrix(net,v0,1)
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
cM, eV = PhyloTraits.gaussiancoalescent_covariancematrix(net,v0,1)
@test tiplabels(eV) == tiplabels(cM) == ["C","d1","d2","B","A"]
@test eV[:tips] + diag(cM[:tips]) ≈ [2.7,3.02,2.396,3.1,2.6]
v11=1.7668462203929007; v12=0.11761402816197718; v13=0.08199968747807541
v22=2.038420771136349; v23=0.443417290678788; v33=1.247326232249089
v44=2.144808361531078; v45=1.2218017549129516; v55=1.6738764987615091
@test cM[:tips] ≈ [v11 v12 v13 0 0; v12 v22 v23 0 0; v13 v23 v33 0 0;
    0 0 0 v44 v45; 0 0 0 v45 v55]
end

@testset "fit σ2 under the gaussian with coalescent" begin

nwk = "(A:2.5,((B:1,#H1:0.5::0.1):1,(C:1,(D:0.5)#H1:0.5::0.9):1):0.5);"
net = readnewick(nwk)
Y = [8.60,10.56,11.3,9.96,11.24,]
X = ones(5, 1)
df = DataFrame(trait = Y, tipnames = ["B","C","A","D","A",])
f0 = phylolm(@formula(trait ~ 1), df, net; model="gaussiancoalescent",
        reml=false, startingValue=1.0, fixedValue=0.1)

end
