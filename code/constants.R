# **********************
#
# Global variables used throughout project code
#
# **********************


# install and/or load necessary packages
local.packages <- installed.packages()[,"Package"]
if(!("INLA" %in% local.packages))
  install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/stable")
local.packages <- installed.packages()[,"Package"]
req.packages <- c("xtable", "data.table", "matrixStats", "RColorBrewer", "fields", "animation")
needed.packages <- req.packages[!(req.packages %in% local.packages)]
if(length(needed.packages) > 0) 
  install.packages(needed.packages)

library(xtable)
library(data.table)
library(matrixStats)
library(INLA)
library(RColorBrewer)
library(fields)
library(animation)

#Creates functional basis "mesh"
mesh_boundary <- inla.mesh.segment(matrix(c(0,0, 47,0, 47,50, 0,50), 4, 2, byrow=T))
mesh <- inla.mesh.create.helper(points=as.matrix(rbind(expand.grid(2:32, seq(2, 48, 2)))),
                                boundary=mesh_boundary,
                                offset=c(1,25),
                                cutoff=2, max.edge=c(10,40))
mesh.proj <- inla.mesh.projector(mesh, dims=c(200,200))

# number of basis functions
n.basis <- 10

# formulas for hierarchical models (covariates only)
formul.take <- y ~ dribble + ndef + ball.lastsec
formul.TO <- y ~ dribble + ndef + ball.lastsec
formul.make <- y ~ dribble + ndef
formul.pass <- y ~ dribble + ndef + ball.lastsec + doff1 + doff2 + doff3 + ddef

# COARSENED STATES
# note, entityPlace covariate mapping function needs to correspond with these states
states <- expand.grid(c("behind", "per", "rest", "key", "cor3", "cen3", "other"), c(TRUE, FALSE))
state_nms <- paste(states[,1], states[,2], sep="-")
cont.limit <- log(1 + 5) #5 feet qualifies as "contested"

# macrotransitoin models fit
inla.names <- c("TAKE", "MAKE", "PASS1", "PASS2", "PASS3", "PASS4", "TO")
inlas.covariates <- list(labels(terms(formul.take)), 
                         labels(terms(formul.make)), 
                         labels(terms(formul.pass)), 
                         labels(terms(formul.pass)), 
                         labels(terms(formul.pass)), 
                         labels(terms(formul.pass)), 
                         labels(terms(formul.TO)))
temp <- grep("PASS", inla.names)
for(i in 1:length(temp)) {
  temp2 <- inlas.covariates[[temp[i]]]
  inlas.covariates[[temp[i]]][length(temp2)] <- sprintf("ddef%i", i)
}

# function bases for ech macro entroy model
load(sprintf("%s/spatial_bases.Rdata", data.dir))
bases <- list(take.basis, take.basis, pass1.basis, pass2.basis, pass3.basis, pass4.basis, TO.basis)
