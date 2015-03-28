library(ape)
library(phangorn)
#reading the data matrix and the trees 
matriz <- read.phyDat("yeti.phy", format="phylip", type = "DNA")
arbolparsimony <- read.tree('parsimony.tre')
arbolnepal <- read.tree('nepal.tre')
#Optimizing based on genetic distance.
dm <- dist.logDet(matriz)
arbolparsimony <- rtree(13,tip.label = arbolparsimony$tip.label)
arbolparsimony <- NJ(dm)
arbolnepal <- rtree(13,tip.label = arbolnepal$tip.label)
arbolnepal <- unroot(upgma(dm))
#Optimizing edge weights
fit1 <- pml(arbolparsimony, matriz) 
fit2 <- pml(arbolnepal, matriz )
fit1 <- optim.pml(fit1)
fit2 <- optim.pml(fit2)
#Running the SH test.
SH.test(fit1, fit2, B = 100)
#optimizing edge weights from a model, in this case HKY85 
fitHKY85pars <- update(fit1, k=4)
fitHKY85pars = optim.pml(fitHKY85pars, TRUE,TRUE, TRUE, TRUE, TRUE, control = pml.control(trace = 0))
fitHKY85nepa <- update(fit2, k=4)
fitHKY85nepa = optim.pml(fitHKY85nepa, TRUE,TRUE, TRUE, TRUE, TRUE, control = pml.control(trace = 0))

SH.test(fitHKY85pars, fitHKY85nepa, B = 100)

