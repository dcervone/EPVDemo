library(fields)

#requires R-INLA... visit website for details
library(INLA)

#Creates functional basis "mesh"
#DO NOT TOUCH THESE LINES!!!!!
#
# 4/25 MESH
boundary <- matrix(c(0,0, 47,0, 47,50, 0,50), 4, 2, byrow=T)
mesh_boundary <- inla.mesh.segment(boundary)
pts <- as.matrix(rbind(expand.grid(2:32, seq(2, 48, 2))))
mesh <- inla.mesh.create.helper(points=pts,boundary=mesh_boundary,offset=c(1,25),
                                cutoff=2,max.edge=c(10,40))
materncov <- inla.spde.create(mesh, model="matern")
mesh.proj <- inla.mesh.projector(mesh, dims=c(200,200))
