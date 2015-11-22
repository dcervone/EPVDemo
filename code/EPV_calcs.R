# **********************
#
# Functions to calculate EPV and related quantities given parameter
# estimates
#
# **********************

tdatCovars <- function(tdat) {
  # Makes covariate and spatial effects multiplier matrix
  # Inputs:
  #   tdat: 'tdat' formatted data frame
  # Outputs:
  #   dat.fix: covariate matrix W_j(t)
  #   dat.spat: matrix of spatial basis loadings
  #   dat.state: matrix of coarsened state information
  
  # build multiplier matrices from data and hyperparams
  dat.fix <- vector("list", length(inla.names))
  dat.spat <- vector("list", length(inla.names))
  for(j in 1:length(inla.names)) {
    dat.fix[[j]] <- as.matrix(tdat[, inlas.covariates[[j]]])
    if (j %in% grep("PASS", inla.names)) {
      temp <- grep("[0-9]", strsplit(inla.names[j], NULL)[[1]], value=TRUE)
      spat <- cBind(inla.mesh.projector(mesh, loc=as.matrix(tdat[,c("x", "y")]))$proj$A,
                    inla.mesh.projector(mesh, loc=as.matrix(tdat[,c(sprintf("off%s_x", temp), sprintf("off%s_y", temp))]))$proj$A)
    } else {
      spat <- inla.mesh.projector(mesh, loc=as.matrix(tdat[,c("x", "y")]))$proj$A      
    }
    dat.spat[[j]] <- spat %*% t(bases[[j]])
  }
  
  # coarsened state information for transition options
  eP_pass1 <- as.vector(apply(tdat[,c("off1_x", "off1_y")], 1, function(r) entityPlace(r[1], r[2])))
  eP_pass2 <- as.vector(apply(tdat[,c("off2_x", "off2_y")], 1, function(r) entityPlace(r[1], r[2]))) 
  eP_pass3 <- as.vector(apply(tdat[,c("off3_x", "off3_y")], 1, function(r) entityPlace(r[1], r[2]))) 
  eP_pass4 <- as.vector(apply(tdat[,c("off4_x", "off4_y")], 1, function(r) entityPlace(r[1], r[2])))
  
  this_cont <- tdat$ndef < cont.limit
  cont_pass1 <- getNdef(tdat, "off1_x", "off1_y") < cont.limit
  cont_pass2 <- getNdef(tdat, "off2_x", "off2_y") < cont.limit
  cont_pass3 <- getNdef(tdat, "off3_x", "off3_y") < cont.limit
  cont_pass4 <- getNdef(tdat, "off4_x", "off4_y") < cont.limit
  
  dat.state <- data.frame(eP_pass1, eP_pass2, eP_pass3, eP_pass4, 
                          this_cont, 
                          cont_pass1, cont_pass2, cont_pass3, cont_pass4)
  return(list(dat.fix=dat.fix, dat.spat=dat.spat,
              dat.state=dat.state))
}


fvToEPV <- function(tdat, tdat.covars, fv, ev.out) {
  # transforms hazards into instataneous probabilities and EVs
  # Inputs:
  #   tdat: 'tdat' formatted data frame
  #   tdat.covars: multiplier matrices from tdatCovars()
  #   fv: hazards data frame from fvMatToDF()
  #   ev.out: coarsened state EVs from evLineups()
  # Outputs:
  #   list with data frame of instataneous transition probabilities
  #   and their (instantaneous) expected values
  
  p.do <- epv.do <- rep(NA, nrow(tdat))
  ppass1 <- epass1 <- rep(NA, nrow(tdat))
  ppass2 <- epass2 <- rep(NA, nrow(tdat))
  ppass3 <- epass3 <- rep(NA, nrow(tdat))
  ppass4 <- epass4 <- rep(NA, nrow(tdat))
  pmake <- emake <- rep(NA, nrow(tdat))
  pmiss <- emiss <- rep(NA, nrow(tdat))
  pTO <- eTO <- rep(NA, nrow(tdat))
  enow <- rep(NA, nrow(tdat))
  
  ppass1 <- fv$fvs_pass1
  ppass2 <- fv$fvs_pass2
  ppass3 <- fv$fvs_pass3
  ppass4 <- fv$fvs_pass4
  pmake <- fv$fvs_take * fv$fvs_make
  pmiss <- fv$fvs_take * (1 - fv$fvs_make)
  pTO <- fv$fvs_TO
  
  # coarsened state EV of each possible transition
  evs <- ev.out$evs
  for(i in 1:nrow(ev.out$teammates.all)){
    these <- ev.out$lineup.inds[[i]]
    adder <- ev.out$adder.mat[i, ]
    statepass1 <- match(paste(tdat.covars$dat.state$eP_pass1[these], tdat.covars$dat.state$cont_pass1[these], sep="-"), state_nms)
    epass1[these] <- evs[[i]][adder[1]+statepass1]
    
    statepass2 <- match(paste(tdat.covars$dat.state$eP_pass2[these], tdat.covars$dat.state$cont_pass2[these], sep="-"), state_nms)
    epass2[these] <- evs[[i]][adder[2]+statepass2]
    
    statepass3 <- match(paste(tdat.covars$dat.state$eP_pass3[these], tdat.covars$dat.state$cont_pass3[these], sep="-"), state_nms)
    epass3[these] <- evs[[i]][adder[3]+statepass3]
    
    statepass4 <- match(paste(tdat.covars$dat.state$eP_pass4[these], tdat.covars$dat.state$cont_pass4[these], sep="-"), state_nms)
    epass4[these] <- evs[[i]][adder[4]+statepass4]
    
    state.now <- match(paste(tdat$eP[these], tdat$ndef[these] < cont.limit, sep="-"), state_nms)
    enow[these] <- evs[[i]][adder[5]+state.now]  
  }
  
  emake <- 2 * (tdat$threept == 0) + 3 * (tdat$threept == 1)
  emiss <- rep(.15, nrow(tdat))
  eTO <- rep(0, nrow(tdat))
  return(list(probs=data.frame(ppass1, ppass2, ppass3, ppass4, 
                               pmake, pmiss, pTO),
              vals=data.frame(epass1, epass2, epass3, epass4,
                              emake, emiss, eTO, enow)))
}

fitVals <- function(hyper, tdat, tdat.covars) {
  # Calculates hazards given parameter estimates and multiplier matrices
  # Inputs:
  #   hyper: list of hyperparameters from getHyperParams()
  #   tdat: 'tdat' formatted data frame
  #   tdat.covars: multiplier matrices from tdatCovars()
  # Outputs:
  #   fvs: matrix of hazards
  
  fvs <- matrix(NA, nrow=nrow(tdat), ncol=length(inla.names))
  off0.idx <- match(tdat$entity, hyper$player.ids)
  for (i in 1:length(inla.names)) {
    for(j in na.omit(unique(off0.idx))) {
      inds <- which(off0.idx == j)
      dat.fix <- tdat.covars$dat.fix[[i]][inds, ]
      dat.spat <- tdat.covars$dat.spat[[i]][inds, ]
      y <- hyper$macro.means[[j]][[i]]
      fvs[inds, i] <- as.numeric(as.vector(y[1] + dat.fix %*% y[1 + 1:ncol(dat.fix)] + dat.spat %*% tail(y, ncol(dat.spat))))
    }
  }
  return(fvs)
}

fvMatToDF <- function(fv) {
  # Formats hazard matrix to data frame

  fv <- as.data.frame(exp(fv))
  names(fv) <- c("fvs_take", "fvs_make", "fvs_pass1", 
                 "fvs_pass2", "fvs_pass3", "fvs_pass4", 
                 "fvs_TO")
  fv$fvs_make <- fv$fvs_make / (1 + fv$fvs_make)
  over1 <- which(fv > .99, arr.ind=T)
  if(nrow(over1) > 0) {
    fv[over1] <- .99
  }
  return(fv)
}

multiresDraw <- function(tdat, hyper, micro.def.mod, ev.out, nmic=50, save.positions=FALSE) {
  # Multiresolution transition simulation for a single EPV draw
  # Inputs:
  #   tdat: 'tdat' formatted data frame
  #   hyper: list of hyperparameters from getHyperParams()
  #   micro.def.mod: defensive micro model output of microDefModel()
  #   ev.out: coarsened state EVs from evLineups()
  #   nmic: number of microtransitions to simulate ahead (25 = 1 second)
  #   save.positions: save player positions generated with each micro update?
  # Outputs:
  #   fv.epv.list: list of instantaneous transition probabilities/values
  #       for each nmic simulated microtransition
  #   positions: list of simulated player positions if save.positions=TRUE
  
  # dynamic compenents for micro model
  off.x <- as.matrix(tdat[,c("x", "off1_x", "off2_x", "off3_x", "off4_x")])
  off.y <- as.matrix(tdat[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")])
  off.eps.x <- rbind(rep(0, 5), off.x[-1, ] - off.x[-nrow(off.x), ])
  off.eps.y <- rbind(rep(0, 5), off.y[-1, ] - off.y[-nrow(off.y), ])
  off.eps.x[which(abs(off.eps.x) > 2, arr.ind=T)] <- 0
  off.eps.y[which(abs(off.eps.y) > 2, arr.ind=T)] <- 0
  
  def.x <- as.matrix(tdat[,c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")])
  def.y <- as.matrix(tdat[,c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")])
  def.eps.x <- rbind(rep(0, 5), def.x[-1, ] - def.x[-nrow(def.x), ])
  def.eps.y <- rbind(rep(0, 5), def.y[-1, ] - def.y[-nrow(def.y), ])
  def.eps.x[which(abs(def.eps.x) > 2, arr.ind=T)] <- 0
  def.eps.y[which(abs(def.eps.y) > 2, arr.ind=T)] <- 0
  
  fv.epv.list <- vector("list", nmic)
  positions <- vector("list", nmic)
  if(save.positions) {
    positions[[1]] <- list(off.x=off.x, off.y=off.y, def.x=def.x, def.y=def.y)
  }
  
  # get hazards
  tdat.covars <- tdatCovars(tdat)
  fv <- fvMatToDF(fitVals(hyper, tdat, tdat.covars))
  fv.epv.list[[1]] <- fvToEPV(tdat, tdat.covars, fv, ev.out)
  for(i in 2:nmic) {
    print(sprintf("micro %i of %i", i, nmic))
    # update player positions from micro model
    micro.out <- microtdat(micro.def.mod, hyper,
                            tdat, off.eps.x, off.eps.y,
                            def.eps.x, def.eps.y)
    off.eps.x <- micro.out$off.eps.x
    off.eps.y <- micro.out$off.eps.y
    def.eps.x <- micro.out$def.eps.x
    def.eps.y <- micro.out$def.eps.y
    if(save.positions) {
      off.x <- as.matrix(micro.out$tdat[,c("x", "off1_x", "off2_x", "off3_x", "off4_x")])
      off.y <- as.matrix(micro.out$tdat[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")])
      
      def.x <- as.matrix(micro.out$tdat[,c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")])
      def.y <- as.matrix(micro.out$tdat[,c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")])
      positions[[i]] <- list(off.x=off.x, off.y=off.y, def.x=def.x, def.y=def.y)
    }
    
    # get hazards from updates positions
    tdat.covars <- tdatCovars(micro.out$tdat)
    fv <- fvMatToDF(fitVals(hyper, micro.out$tdat, tdat.covars))
    fv.epv.list[[i]] <- fvToEPV(micro.out$tdat, tdat.covars, fv, ev.out)
    tdat <- micro.out$tdat
  }
  return(list(fv.epv.list=fv.epv.list, positions=positions))
}

compressEPV <- function(tdat, fv.epv.list) {
  # comporesses multiresDraw output into a single EPV draw
  # Inputs:
  #   tdat: 'tdat' formatted data frame
  #   fv.epv.list: fv.epv.list component of output of multiresDraw()
  # Outputs:
  #   EPV, and marginal transition probabilities + values
  
  # integrate over time to obtain marginal transition probabilities and values
  nmic <- length(fv.epv.list)
  n <- nrow(fv.epv.list[[1]]$probs)
  p.here <- rep(1, n)
  ppass1 <- ppass2 <- ppass3 <- ppass4 <- pmake <- pmiss <- pTO <- rep(0, n)
  epass1 <- epass2 <- epass3 <- epass4 <- emake <- emiss <- eTO <- rep(0, n)
  for(i in 1:nmic) {
    probs <- fv.epv.list[[i]]$probs
    vals <- fv.epv.list[[i]]$vals
    p.row <- rowSums(probs)
    probs <- probs / p.row
    over1 <- which(p.row > 1)
    if(length(over1) > 0)
      p.row[over1] <- 1
    
    #p.do <- 1 - apply(1 - fv.epv.list[[i]]$probs, 1, prod)
    ppass1 <- ppass1 + probs$ppass1 * p.here * p.row
    ppass2 <- ppass2 + probs$ppass2 * p.here * p.row
    ppass3 <- ppass3 + probs$ppass3 * p.here * p.row
    ppass4 <- ppass4 + probs$ppass4 * p.here * p.row
    pmake <- pmake + probs$pmake * p.here * p.row
    pmiss <- pmiss + probs$pmiss * p.here * p.row
    pTO <- pTO + probs$pTO * p.here * p.row
    epass1 <- epass1 + vals$epass1 * probs$ppass1 * p.here * p.row
    epass2 <- epass2 + vals$epass2 * probs$ppass2 * p.here * p.row
    epass3 <- epass3 + vals$epass3 * probs$ppass3 * p.here * p.row
    epass4 <- epass4 + vals$epass4 * probs$ppass4 * p.here * p.row
    emake <- emake + vals$emake * probs$pmake * p.here * p.row
    emiss <- emiss + vals$emiss * probs$pmiss * p.here * p.row
    eTO <- eTO + vals$eTO * probs$pTO * p.here * p.row
    p.here <- p.here * (1 - p.row)
    eother <- vals$enow
  }
  
  epv <- (epass1 + epass2 + epass3 + epass4 + emake + emiss + eTO + eother * p.here)
  epass1 <- epass1 / ppass1
  epass2 <- epass2 / ppass2
  epass3 <- epass3 / ppass3
  epass4 <- epass4 / ppass4
  emake <- emake / pmake
  emiss <- emiss / pmiss
  eTO <- eTO / pTO
  
  probs.now <- fv.epv.list[[1]]$probs
  vals.now <- fv.epv.list[[1]]$vals
  names(probs.now) <- c("pass1", "pass2", "pass3", "pass4", "make", "miss", "TO")
  names(vals.now) <- c("pass1", "pass2", "pass3", "pass4", "make", "miss", "TO", "state")
  
  # when events happen, EPV is conditioned on these transitions
  epv[which(tdat$event_id %in% c(3, 4))] <- ((fv.epv.list[[1]]$vals$emake * fv.epv.list[[1]]$probs$pmake + 
                                                 fv.epv.list[[1]]$vals$emiss * fv.epv.list[[1]]$probs$pmiss) / (fv.epv.list[[1]]$probs$pmake + fv.epv.list[[1]]$probs$pmiss))[which(tdat$event_id %in% c(3, 4))]
  #get pass recipient:
  pass.idx <- which(tdat$event_id %in% c(22, 25))
  recips <- tdat[pass.idx, c("off1_ent", "off2_ent", "off3_ent", "off4_ent")]
  pass.ent <- sapply(seq(pass.idx), function(j) match(tdat$entity[pass.idx[j] + 1], recips[j, ]))
  epv[pass.idx[which(pass.ent == 1)]] <- fv.epv.list[[1]]$vals$epass1[pass.idx[which(pass.ent == 1)]]
  epv[pass.idx[which(pass.ent == 2)]] <- fv.epv.list[[1]]$vals$epass2[pass.idx[which(pass.ent == 2)]]
  epv[pass.idx[which(pass.ent == 3)]] <- fv.epv.list[[1]]$vals$epass3[pass.idx[which(pass.ent == 3)]]
  epv[pass.idx[which(pass.ent == 4)]] <- fv.epv.list[[1]]$vals$epass4[pass.idx[which(pass.ent == 4)]]
  epv[which(tdat$event_id == 7)] <- fv.epv.list[[1]]$vals$eTO[which(tdat$event_id == 7)]
  return(list(epv=epv, probs=data.frame(pass1=ppass1, pass2=ppass2, pass3=ppass3, pass4=ppass4,
                                        make=pmake, miss=pmiss, TO=pTO, other=p.here),
              vals=data.frame(pass1=epass1, pass2=epass2, pass3=epass3, pass4=epass4,
                              make=emake, miss=emiss, TO=eTO, other=eother),
              probs.now=probs.now,
              vals.now=vals.now))
}

microtdat <- function(micro.def.mod, hyper, 
                       tdat, 
                       off.eps.x, off.eps.y, 
                       def.eps.x, def.eps.y){
  # Updates player positions using micro model
  # Inputs:
  #   tdat: 'tdat' formatted data frame
  #   hyper: list of hyperparameters from getHyperParams()
  #   micro.def.mod: defensive micro model output of microDefModel()
  #   off/def.eps.x/y: differenced positions
  # Ouputs:
  #   tdat: 'tdat' updated with position innovations
  #   off/def.eps.x/y: predicted differenced positions
  
  off.x <- as.matrix(tdat[,c("x", "off1_x", "off2_x", "off3_x", "off4_x")])
  off.y <- as.matrix(tdat[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")])
  
  # match players to parameters of micro update model
  off0.idx <- match(tdat$entity, hyper$player.ids)
  off1.idx <- match(tdat$off1_ent, hyper$player.ids)
  off2.idx <- match(tdat$off2_ent, hyper$player.ids)
  off3.idx <- match(tdat$off3_ent, hyper$player.ids)
  off4.idx <- match(tdat$off4_ent, hyper$player.ids)
  
  intercepts.x.0 <- hyper$intercepts.x.with[off0.idx]
  intercepts.x.1 <- hyper$intercepts.x.without[off1.idx]
  intercepts.x.2 <- hyper$intercepts.x.without[off2.idx]
  intercepts.x.3 <- hyper$intercepts.x.without[off3.idx]
  intercepts.x.4 <- hyper$intercepts.x.without[off4.idx]
  coefs.x.0 <- hyper$coefs.x.with[off0.idx]
  coefs.x.1 <- hyper$coefs.x.without[off1.idx]
  coefs.x.2 <- hyper$coefs.x.without[off2.idx]
  coefs.x.3 <- hyper$coefs.x.without[off3.idx]
  coefs.x.4 <- hyper$coefs.x.without[off4.idx]
  sigmas.x.0 <- hyper$sigmas.x.with[off0.idx]
  sigmas.x.1 <- hyper$sigmas.x.without[off1.idx]
  sigmas.x.2 <- hyper$sigmas.x.without[off2.idx]
  sigmas.x.3 <- hyper$sigmas.x.without[off3.idx]
  sigmas.x.4 <- hyper$sigmas.x.without[off4.idx]
  new.x.mu.0 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$x, tdat$y)))$proj$A * t(hyper$spats.x.with[, off0.idx]))
  new.x.mu.1 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$off1_x, tdat$off1_y)))$proj$A * t(hyper$spats.x.without[, off1.idx]))
  new.x.mu.2 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$off2_x, tdat$off2_y)))$proj$A * t(hyper$spats.x.without[, off2.idx]))
  new.x.mu.3 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$off3_x, tdat$off3_y)))$proj$A * t(hyper$spats.x.without[, off3.idx]))
  new.x.mu.4 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$off4_x, tdat$off4_y)))$proj$A * t(hyper$spats.x.without[, off4.idx]))
  off.new.x.0 <- intercepts.x.0 + tdat$x + new.x.mu.0 + coefs.x.0 * off.eps.x[, 1] + rnorm(length(sigmas.x.0), 0, sigmas.x.0)
  off.new.x.1 <- intercepts.x.1 + tdat$off1_x + new.x.mu.1 + coefs.x.1 * off.eps.x[, 2] + rnorm(length(sigmas.x.1), 0, sigmas.x.1)
  off.new.x.2 <- intercepts.x.2 + tdat$off2_x + new.x.mu.2 + coefs.x.2 * off.eps.x[, 3] + rnorm(length(sigmas.x.2), 0, sigmas.x.2)
  off.new.x.3 <- intercepts.x.3 + tdat$off3_x + new.x.mu.3 + coefs.x.3 * off.eps.x[, 4] + rnorm(length(sigmas.x.3), 0, sigmas.x.3)
  off.new.x.4 <- intercepts.x.4 + tdat$off4_x + new.x.mu.4 + coefs.x.4 * off.eps.x[, 5] + rnorm(length(sigmas.x.4), 0, sigmas.x.4)
  
  intercepts.y.0 <- hyper$intercepts.y.with[off0.idx]
  intercepts.y.1 <- hyper$intercepts.y.without[off1.idx]
  intercepts.y.2 <- hyper$intercepts.y.without[off2.idx]
  intercepts.y.3 <- hyper$intercepts.y.without[off3.idx]
  intercepts.y.4 <- hyper$intercepts.y.without[off4.idx]
  coefs.y.0 <- hyper$coefs.y.with[off0.idx]
  coefs.y.1 <- hyper$coefs.y.without[off1.idx]
  coefs.y.2 <- hyper$coefs.y.without[off2.idx]
  coefs.y.3 <- hyper$coefs.y.without[off3.idx]
  coefs.y.4 <- hyper$coefs.y.without[off4.idx]
  sigmas.y.0 <- hyper$sigmas.y.with[off0.idx]
  sigmas.y.1 <- hyper$sigmas.y.without[off1.idx]
  sigmas.y.2 <- hyper$sigmas.y.without[off2.idx]
  sigmas.y.3 <- hyper$sigmas.y.without[off3.idx]
  sigmas.y.4 <- hyper$sigmas.y.without[off4.idx]
  new.y.mu.0 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$x, tdat$y)))$proj$A * t(hyper$spats.y.with[, off0.idx]))
  new.y.mu.1 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$off1_x, tdat$off1_y)))$proj$A * t(hyper$spats.y.without[, off1.idx]))
  new.y.mu.2 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$off2_x, tdat$off2_y)))$proj$A * t(hyper$spats.y.without[, off2.idx]))
  new.y.mu.3 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$off3_x, tdat$off3_y)))$proj$A * t(hyper$spats.y.without[, off3.idx]))
  new.y.mu.4 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(tdat$off4_x, tdat$off4_y)))$proj$A * t(hyper$spats.y.without[, off4.idx]))
  off.new.y.0 <- intercepts.y.0 + tdat$y + new.y.mu.0 + coefs.y.0 * off.eps.y[, 1] + rnorm(length(sigmas.y.0), 0, sigmas.y.0)
  off.new.y.1 <- intercepts.y.1 + tdat$off1_y + new.y.mu.1 + coefs.y.1 * off.eps.y[, 2] + rnorm(length(sigmas.y.1), 0, sigmas.y.1)
  off.new.y.2 <- intercepts.y.2 + tdat$off2_y + new.y.mu.2 + coefs.y.2 * off.eps.y[, 3] + rnorm(length(sigmas.y.2), 0, sigmas.y.2)
  off.new.y.3 <- intercepts.y.3 + tdat$off3_y + new.y.mu.3 + coefs.y.3 * off.eps.y[, 4] + rnorm(length(sigmas.y.3), 0, sigmas.y.3)
  off.new.y.4 <- intercepts.y.4 + tdat$off4_y + new.y.mu.4 + coefs.y.4 * off.eps.y[, 5] + rnorm(length(sigmas.y.4), 0, sigmas.y.4)
  
  mod.x <- micro.def.mod$mod.x
  mod.y <- micro.def.mod$mod.y
  
  opts.x <- tdat$ball_x*.27 + off.x*.62 + 4.75*.11
  opts.y <- tdat$ball_y*.27 + off.y*.62 + 25*.11
  def.x <- as.matrix(tdat[,c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")])
  def.y <- as.matrix(tdat[,c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")])
  
  # update defensive matchup component
  wts <- array(NA, dim=c(nrow(off.x), 5, 5))
  for(i in 1:nrow(off.x)){
    def.x.n <- def.x[i,]
    def.y.n <- def.y[i,]
    guarding <- sapply(1:5, function(j) {
      #wts <- exp(-.3*(opts.x[i,] - def.x.n[j])^2 - .3*(opts.y[i,] - def.y.n[j])^2)
      wts <- 1/((opts.x[i,] - def.x.n[j])^2 + (opts.y[i,] - def.y.n[j])^2)
      wts <- wts/sum(wts)
      return(wts)
    })
    wts[i,,] <- guarding
  }
  
  # new positions
  off.new.x <- cbind(off.new.x.0, off.new.x.1, off.new.x.2, off.new.x.3, off.new.x.4)
  off.new.y <- cbind(off.new.y.0, off.new.y.1, off.new.y.2, off.new.y.3, off.new.y.4)
  
  new.opts.x <- tdat$ball_x*.27 + off.new.x*.62 + 4.75*.11
  new.opts.y <- tdat$ball_y*.27 + off.new.y*.62 + 25*.11
  
  def.new.x <- mod.x$coef[1] + def.x + 
    mod.x$coef[2]*def.eps.x + 
    mod.x$coef[4] * t(sapply(1:nrow(off.x), 
                             function(i) colSums((new.opts.x[i, ] - opts.x[i, ])*wts[i, , ]))) + 
    mod.x$coef[3] * (def.x - t(sapply(1:nrow(def.x), function(i) colSums(new.opts.x[i, ] * wts[i, , ])))) + 
    summary(mod.x)$sigma * matrix(rnorm(prod(dim(def.x))), nrow=nrow(def.x), ncol=ncol(def.x))
  def.new.y <- mod.y$coef[1] + def.y + 
    mod.y$coef[2]*def.eps.y + 
    mod.y$coef[4] * t(sapply(1:nrow(off.y), 
                             function(i) colSums((new.opts.y[i, ] - opts.y[i, ])*wts[i, , ]))) + 
    mod.y$coef[3] * (def.y - t(sapply(1:nrow(def.y), function(i) colSums(new.opts.y[i, ] * wts[i, , ])))) + 
    summary(mod.y)$sigma * matrix(rnorm(prod(dim(def.y))), nrow=nrow(def.y), ncol=ncol(def.y))
  
  def.eps.x <- def.new.x - def.x
  def.eps.y <- def.new.y - def.y
  off.eps.x <- off.new.x - off.x
  off.eps.y <- off.new.y - off.y
  
  off.new.x[which(off.new.x < 0, arr.ind=T)] <- 0
  off.new.y[which(off.new.y < 0, arr.ind=T)] <- 0
  def.new.x[which(def.new.x < 0, arr.ind=T)] <- 0
  def.new.y[which(def.new.y < 0, arr.ind=T)] <- 0
  off.new.x[which(off.new.x > 94, arr.ind=T)] <- 94
  off.new.y[which(off.new.y > 50, arr.ind=T)] <- 50
  def.new.x[which(def.new.x > 94, arr.ind=T)] <- 94
  def.new.y[which(def.new.y > 50, arr.ind=T)] <- 50
  
  ball.new.x <- tdat$ball_x + off.eps.x[, 1]
  ball.new.y <- tdat$ball_y + off.eps.y[, 1]
  
  
  # update covariate summaries of positional data
  tdat[,c("x", "off1_x", "off2_x", "off3_x", "off4_x")] <- off.new.x
  tdat[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")] <- off.new.y
  tdat[,c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")] <- def.new.x
  tdat[,c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")] <- def.new.y
  tdat[, "ball_x"] <- ball.new.x
  tdat[, "ball_y"] <- ball.new.y
  tdat[,"ndef"] <- getNdef(tdat, "x", "y")
  tdat[, "eP"] <- getPlace(tdat) 
  threept.out <- getThree(tdat)
  tdat[, "bask.dist"] <- threept.out$bask.dist
  tdat[, "threept"] <- threept.out$threept
  return(list(tdat=tdat, off.eps.x=off.eps.x, off.eps.y=off.eps.y, 
              def.eps.x=def.eps.x, def.eps.y=def.eps.y))
}

# is this function below necessary??? pls no
doEPV <- function(tmats, use.leagavg=0, fullres=TRUE){
  evs <- list()
  for(i in 1:nrow(teammates.all)){
    # if(i %% 50 == 0)
    #   print(sprintf("%i of %i", i, nrow(teammates.all)))
    temp <- tryCatch(calcEV(tmats, teammates.all[i,], use.leagavg), error = function(e) e)
    if(!(inherits(temp, "error"))) evs[[i]] <- temp else evs[[i]] <- NA
  }
  
  p.do <- epv.do <- rep(NA, nrow(tdat))
  ppass1 <- epass1 <- rep(NA, nrow(tdat))
  ppass2 <- epass2 <- rep(NA, nrow(tdat))
  ppass3 <- epass3 <- rep(NA, nrow(tdat))
  ppass4 <- epass4 <- rep(NA, nrow(tdat))
  pmake <- emake <- rep(NA, nrow(tdat))
  pmiss <- emiss <- rep(NA, nrow(tdat))
  pTO <- eTO <- rep(NA, nrow(tdat))
  
  if(id %in% use.leagavg) fvs <- fvs.pos else fvs <- fvs.ind
  
  ppass1 <- fvs$fvs_pass1
  ppass2 <- fvs$fvs_pass2
  ppass3 <- fvs$fvs_pass3
  ppass4 <- fvs$fvs_pass4
  
  pmake <- fvs$fvs_take*fvs$fvs_make
  pmiss <- fvs$fvs_take*(1-fvs$fvs_make)
  pTO <- fvs$fvs_TO
  
  for(i in 1:nrow(teammates.all)){
    if(is.na(evs[[i]])[1]) next
    these <- which(tdat$off1_ent == teammates[i,1] & tdat$off2_ent == teammates[i,2] & 
                     tdat$off3_ent == teammates[i,3] & tdat$off4_ent == teammates[i,4])
    place <- which(teammates.all[i,] == id)
    if(place == 1){
      adder <- c(14,28,42,56)
    } else if(place == 2){
      adder <- c(0,28,42,56)
    } else if(place == 3){
      adder <- c(0,14,42,56)
    } else if(place == 4){
      adder <- c(0,14,28,56)
    } else if(place == 5){
      adder <- c(0,14,28,42)
    } 
    statepass1 <- match(paste(eP_pass1[these], cont_pass1[these], sep="-"), state_nms)
    epass1[these] <- evs[[i]][adder[1]+statepass1]
    
    statepass2 <- match(paste(eP_pass2[these], cont_pass2[these], sep="-"), state_nms)
    epass2[these] <- evs[[i]][adder[2]+statepass2]
    
    statepass3 <- match(paste(eP_pass3[these], cont_pass3[these], sep="-"), state_nms)
    epass3[these] <- evs[[i]][adder[3]+statepass3]
    
    statepass4 <- match(paste(eP_pass4[these], cont_pass4[these], sep="-"), state_nms)
    epass4[these] <- evs[[i]][adder[4]+statepass4]
    
    if(!(fullres)){
      statein <- match(paste(tdat$eP[these], tdat$ndef[these] < cont.limit, sep="-"), state_nms)
      epv.do[these] <- evs[[i]][14*(place-1) + statein]
    }
  }
  
  emake <- 2*(tdat$threept == 0) + 3*(tdat$threept == 1)
  emiss <- rep(.15, nrow(tdat))
  eTO <- rep(0, nrow(tdat))
  
  p.do <- ppass1 + ppass2 + ppass3 + ppass4 + pmake + pmiss + pTO
  
  if(fullres){
    epv.do <- (ppass1*epass1 + ppass2*epass2 + ppass3*epass3 + ppass4*epass4 + pmake*emake + pmiss*emiss + pTO*eTO)/p.do
  }
  
  epv.do[which(tdat$event_id %in% c(3,4))] <- ((pmake*emake + pmiss*emiss)/(pmake + pmiss))[which(tdat$event_id %in% c(3,4))]
  epv.do[which(tdat$event_id == 7)] <- 0
  epv.do[which(tdat$event_id == 31)] <- epass1[which(tdat$event_id == 31)]
  epv.do[which(tdat$event_id == 32)] <- epass2[which(tdat$event_id == 32)]
  epv.do[which(tdat$event_id == 33)] <- epass3[which(tdat$event_id == 33)]
  epv.do[which(tdat$event_id == 34)] <- epass4[which(tdat$event_id == 34)]
  
  filter.epv <- epv.do
  thesediff <- diff(tdat$time)
  ends <- c(which(thesediff < 0 | thesediff > 100), nrow(tdat))
  #ends <- ends[which(substr(tdat$game[ends],7,8) %in% c("01","02","03"))]
  starts <- c(1, which(thesediff < 0 | thesediff > 100) + 1)
  #starts <- starts[which(substr(tdat$game[ends],7,8) %in% c("01","02","03"))]
  posses <- lapply(1:length(starts), function(i) starts[i]:ends[i])
  for(i in 1:length(posses)){
    if(length(posses[[i]]) > 1){
      for(j in seq(ends[i]-1, starts[i], -1)){
        if(!(is.na(p.do[j])) & p.do[j] < 1){
          filter.epv[j] <- epv.do[j]*p.do[j] + (1-p.do[j])*filter.epv[j+1]
        }
      }
    }
  }
  
  ppp <- cbind(ppass1, ppass2, ppass3, ppass4)
  epp <- cbind(epass1, epass2, epass3, epass4)
  epass <- rowSums(ppp*epp)/rowSums(ppp)
  shot_sat <- mean((epv.do - epass)[which(tdat$event_id %in% c(3,4))], na.rm=T)
  
  mm <- sapply(1:nrow(ppp), function(i) sum(pmake[i] + pmiss[i] > ppp[i,], na.rm=T))
  eshot <- (pmake*emake + pmiss*emiss)/(pmake + pmiss)
  pass_sat <- mean((epv.do - eshot)[which(mm >= 3 & tdat$event_id %in% c(31,32,33,34))], na.rm=T)
  
  eps <- 1/25
  dpass1 <- diff(c(ppass1,0)) 
  dpass2 <- diff(c(ppass2,0))
  dpass3 <- diff(c(ppass3,0))
  dpass4 <- diff(c(ppass4,0))
  dmake <- diff(c(pmake,0))
  dmiss <- diff(c(pmiss,0))
  dTO <- diff(c(pTO,0))
  pno <- 1 - p.do
  dno <- diff(c(pno,0))
  
  part1 <- (ppass1 + eps*dpass1*pno/(1-eps*dno))*epass1
  part2 <- (ppass2 + eps*dpass2*pno/(1-eps*dno))*epass2
  part3 <- (ppass3 + eps*dpass3*pno/(1-eps*dno))*epass3
  part4 <- (ppass4 + eps*dpass4*pno/(1-eps*dno))*epass4
  part5 <- (pmake + eps*dmake*pno/(1-eps*dno))*emake
  part6 <- (pmiss + eps*dmiss*pno/(1-eps*dno))*emiss
  part7 <- (pTO + eps*dTO*pno/(1-eps*dno))*0
  
  epv.mic <- (part1 + part2 + part3 + part4 + part5 + part6 + part7)
  epv.mic <- epv.mic/(1-pno/(1-eps*dno))
  epv.mic[which(tdat$event_id %in% c(3,4,7,31,32,33,34))] <- epv.do[which(tdat$event_id %in% c(3,4,7,31,32,33,34))]
  
  
  return(data.frame(game=tdat$game, quarter=tdat$quarter, time=tdat$time, game_clock=tdat$game_clock,
                    epv.do=epv.do, epv.mic=epv.mic, epv.filter=filter.epv, ppass1=ppass1,
                    ppass2=ppass2,
                    ppass3=ppass3,
                    ppass4=ppass4,
                    pmake=pmake,
                    pmiss=pmiss,
                    pTO=pTO,
                    epass1=epass1,
                    epass2=epass2,
                    epass3=epass3,
                    epass4=epass4,
                    emake=emake,
                    emiss=emiss,
                    eTO=eTO))
}
