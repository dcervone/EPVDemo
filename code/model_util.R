library(xtable)
library(data.table)
library(matrixStats)
library(fields)
library(INLA)

#Creates functional basis "mesh"
boundary <- matrix(c(0,0, 47,0, 47,50, 0,50), 4, 2, byrow=T)
mesh_boundary <- inla.mesh.segment(boundary)
pts <- as.matrix(rbind(expand.grid(2:32, seq(2, 48, 2))))
mesh <- inla.mesh.create.helper(points=pts,boundary=mesh_boundary,offset=c(1,25),
                                cutoff=2,max.edge=c(10,40))
materncov <- inla.spde.create(mesh, model="matern")
mesh.proj <- inla.mesh.projector(mesh, dims=c(200,200))

# covariates used in hierarchical modeling
n.basis <- 10
spatiotemp.covariates <- c("doff1", "doff2", "doff3", "ddef", "ndef", "ball.lastsec", "dribble")
car.covariates <- c(spatiotemp.covariates,"b0", paste("b", seq(n.basis), sep=""), "wts", "oos")
car.connections <- 8

latent_var <- as.vector(players$position)
latent_var[latent_var == "Point-Guard"] <- 1
latent_var[latent_var == "Guard"] <- 2
latent_var[latent_var == "Shooting-Guard"] <- 3
latent_var[latent_var == "Guard-Forward"] <- 4
latent_var[latent_var == "Small-Forward"] <- 5
latent_var[latent_var == "Forward"] <- 6
latent_var[latent_var == "Power-Forward"] <- 7
latent_var[latent_var == "Forward-Center"] <- 8
latent_var[latent_var == "Center"] <- 9

# formulas for hierarchical models (covariates only)
formul.take <- y ~ dribble + ndef + ball.lastsec
formul.TO <- y ~ dribble + ndef + ball.lastsec
formul.make <- y ~ dribble + ndef
formul.pass <- y ~ dribble + ndef + ball.lastsec + doff1 + doff2 + doff3 + ddef

# COARSENED STATES
# note, entityPlace covariate in covariate_util needs to correspond with these states
states <- expand.grid(c("behind", "per", "rest", "key", "cor3", "cen3", "other"), c(TRUE, FALSE))
state_nms <- paste(states[,1], states[,2], sep="-")
three.pt.states <- c("other", "cor3", "cen3")
cont.limit <- log(1 + 5) #5 feet

# different macrotransition  types, and corresponding
types <- c("TAKE", "MAKE", "PASS1", "PASS2", "PASS3", "PASS4", "TO")
outcomes <- list(c(3,4), 3, 31, 32, 33, 34, 7)
inla.names <- types
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

load(sprintf("%s/spatial_bases.Rdata", data.dir))
bases <- list(take.basis, take.basis, pass1.basis, pass2.basis, pass3.basis, pass4.basis, TO.basis)

microDefModel <- function(datex) {
  # offensive-defensive matchup component
  off.x <- as.matrix(datex[, c("x", "off1_x", "off2_x", "off3_x", "off4_x")])
  opts.x <- datex$ball_x * .27 + off.x * .62 + 4.75 * .11
  off.y <- as.matrix(datex[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")])
  opts.y <- datex$ball_y * .27 + off.y * .62 + 25 * .11
  def.x <- as.matrix(datex[, c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")])
  def.y <- as.matrix(datex[, c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")])
  
  wts <- array(NA, dim=c(nrow(off.x), 5, 5))
  for(i in 1:nrow(off.x)){
    def.x.n <- def.x[i, ]
    def.y.n <- def.y[i, ]
    guarding <- sapply(1:5, function(j) {
      #wts <- exp(-.3*(opts.x[i,] - def.x.n[j])^2 - .3*(opts.y[i,] - def.y.n[j])^2)
      wts <- 1/((opts.x[i, ] - def.x.n[j])^2 + (opts.y[i, ] - def.y.n[j])^2)
      wts <- wts/sum(wts)
      return(wts)
    })
    wts[i, , ] <- guarding
  }
  
  means.x <- t(sapply(1:nrow(off.x), function(j) colSums(opts.x[j, ] * wts[j, , ])))
  means.y <- t(sapply(1:nrow(off.y), function(j) colSums(opts.y[j, ] * wts[j, , ])))
  
  #AR1, difference to optimal, AR1 of optimal
  def.eps.x <- as.vector(diff(def.x))
  residual.x <- as.vector((def.x - means.x)[-1, ])
  opt.eps.x <- as.vector(diff(means.x))
  bads <- which(abs(def.eps.x) > .5)
  residual.x <- residual.x[-bads]
  def.eps.x <- def.eps.x[-bads]
  opt.eps.x <- opt.eps.x[-bads]
  
  mod.x <- lm(def.eps.x[-1] ~ def.eps.x[-length(def.eps.x)] + residual.x[-length(residual.x)] + opt.eps.x[-length(opt.eps.x)])
  
  def.eps.y <- as.vector(diff(def.y))
  residual.y <- as.vector((def.y - means.y)[-1, ])
  opt.eps.y <- as.vector(diff(means.y))
  bads <- which(abs(def.eps.y) > .5)
  residual.y <- residual.y[-bads]
  def.eps.y <- def.eps.y[-bads]
  opt.eps.y <- opt.eps.y[-bads]
  
  mod.y <- lm(def.eps.y[-1] ~ def.eps.y[-length(def.eps.y)] + residual.y[-length(residual.y)] + opt.eps.y[-length(opt.eps.y)])
  return(list(mod.x=mod.x, mod.y=mod.y))
}

getHyperParams <- function(datex) {
  # loads multiresolution transition model parameters for all players in datex
  inlas <- vector("list", length(inla.names))
  for(i in 1:length(inla.names)) {
    # print(sprintf("loading inla output %i of %i", i, length(inla.names)))
    load(sprintf("%s/INLA_%s.Rdata", data.dir, inla.names[i]))
    inlas[[i]] <- inla.out.lite
  }
  player.ids <- unique(unlist(datex[, c("entity", sapply(1:4, function(q) sprintf("off%i_ent", q)),
                                           sapply(1:5, function(q) sprintf("def%i_ent", q)))]))
  player.ids <- player.ids[which(player.ids > 0)]
  micro.inlas <- vector("list",  length(player.ids))
  macro.means <- vector("list",  length(player.ids))
  for(i in 1:length(player.ids)) {
    pl <- player.ids[i]
    load(sprintf("%s/micros/%s.Rdata", data.dir, pl))
    # print(sprintf("loading micro %i out of %i", i, length(player.ids)))
    micro.inlas[[i]] <- list(with=with.ball, without=without.ball)
    macro.means[[i]] <- vector("list", length(inla.names))
    for(j in 1:length(inla.names)) {
      id.n <- match(pl, players$player_id)
      y.fix <- inlas[[j]]$summary.fixed[, "mean"]
      temp <- names(inlas[[j]]$summary.random)
      basis.inds <- c(which(temp == "p.int"), grep("p.b[0-9][0-9]*", temp))
      cov.inds <- setdiff(seq(length(inlas[[j]]$summary.random)), basis.inds)
      y.rand <- c(inlas[[j]]$summary.random$p.int[id.n, "mean"], 
                  sapply(cov.inds, 
                         function(k) inlas[[j]]$summary.random[[k]][id.n, "mean"]),
                  inlas[[j]]$summary.random$p.b1[id.n + nrow(players) * (1:n.basis), "mean"])
      macro.means[[i]][[j]] <- y.fix + y.rand
    }
  }
  
  intercepts.x.with <- sapply(1:length(player.ids), 
                              function(j) micro.inlas[[j]]$with$io.x$summary.fixed["intercept", "mean"])
  intercepts.y.with <- sapply(1:length(player.ids), 
                              function(j) micro.inlas[[j]]$with$io.y$summary.fixed["intercept", "mean"])
  coefs.x.with <- sapply(1:length(player.ids), 
                         function(j) micro.inlas[[j]]$with$io.x$summary.fixed["dif", "mean"])
  coefs.y.with <- sapply(1:length(player.ids), 
                         function(j) micro.inlas[[j]]$with$io.y$summary.fixed["dif", "mean"])
  spats.x.with <- sapply(1:length(player.ids), 
                         function(j) micro.inlas[[j]]$with$io.x$summary.random$spatial[,"mean"])
  spats.y.with <- sapply(1:length(player.ids), 
                         function(j) micro.inlas[[j]]$with$io.y$summary.random$spatial[,"mean"])
  sigmas.x.with <- sapply(1:length(player.ids), 
                          function(j) exp(-.5 * micro.inlas[[j]]$with$io.x$mode$theta[1]))
  sigmas.y.with <- sapply(1:length(player.ids), 
                          function(j) exp(-.5 * micro.inlas[[j]]$with$io.y$mode$theta[1]))
  
  intercepts.x.without <- sapply(1:length(player.ids), 
                                 function(j) micro.inlas[[j]]$without$io.x$summary.fixed["intercept", "mean"])
  intercepts.y.without <- sapply(1:length(player.ids), 
                                 function(j) micro.inlas[[j]]$without$io.y$summary.fixed["intercept", "mean"])
  coefs.x.without <- sapply(1:length(player.ids), 
                            function(j) micro.inlas[[j]]$without$io.x$summary.fixed["dif", "mean"])
  coefs.y.without <- sapply(1:length(player.ids), 
                            function(j) micro.inlas[[j]]$without$io.y$summary.fixed["dif", "mean"])
  spats.x.without <- sapply(1:length(player.ids), 
                            function(j) micro.inlas[[j]]$without$io.x$summary.random$spatial[,"mean"])
  spats.y.without <- sapply(1:length(player.ids), 
                            function(j) micro.inlas[[j]]$without$io.y$summary.random$spatial[,"mean"])
  sigmas.x.without <- sapply(1:length(player.ids), 
                             function(j) exp(-.5 * micro.inlas[[j]]$without$io.x$mode$theta[1]))
  sigmas.y.without <- sapply(1:length(player.ids), 
                             function(j) exp(-.5 * micro.inlas[[j]]$without$io.y$mode$theta[1]))
  
  return(list(player.ids=player.ids, 
              macro.means=macro.means,
              intercepts.x.with=intercepts.x.with,
              intercepts.y.with=intercepts.y.with,
              coefs.x.with=coefs.x.with,
              coefs.y.with=coefs.y.with,
              spats.x.with=spats.x.with,
              spats.y.with=spats.y.with,
              sigmas.x.with=sigmas.x.with,
              sigmas.y.with=sigmas.y.with,
              intercepts.x.without=intercepts.x.without,
              intercepts.y.without=intercepts.y.without,
              coefs.x.without=coefs.x.without,
              coefs.y.without=coefs.y.without,
              spats.x.without=spats.x.without,
              spats.y.without=spats.y.without,
              sigmas.x.without=sigmas.x.without,
              sigmas.y.without=sigmas.y.without))
}

datexCovars <- function(datex) {
  dat.fix <- vector("list", length(inla.names))
  dat.spat <- vector("list", length(inla.names))
  for(j in 1:length(inla.names)) {
    dat.fix[[j]] <- as.matrix(datex[, inlas.covariates[[j]]])
    if (j %in% grep("PASS", inla.names)) {
      temp <- grep("[0-9]", strsplit(inla.names[j], NULL)[[1]], value=TRUE)
      spat <- cBind(inla.mesh.projector(mesh, loc=as.matrix(datex[,c("x", "y")]))$proj$A,
                    inla.mesh.projector(mesh, loc=as.matrix(datex[,c(sprintf("off%s_x", temp), sprintf("off%s_y", temp))]))$proj$A)
    } else {
      spat <- inla.mesh.projector(mesh, loc=as.matrix(datex[,c("x", "y")]))$proj$A      
    }
    dat.spat[[j]] <- spat%*%t(bases[[j]])
  }
  eP_pass1 <- as.vector(apply(datex[,c("off1_x", "off1_y")], 1, function(r) entityPlace(r[1], r[2])))
  eP_pass2 <- as.vector(apply(datex[,c("off2_x", "off2_y")], 1, function(r) entityPlace(r[1], r[2]))) 
  eP_pass3 <- as.vector(apply(datex[,c("off3_x", "off3_y")], 1, function(r) entityPlace(r[1], r[2]))) 
  eP_pass4 <- as.vector(apply(datex[,c("off4_x", "off4_y")], 1, function(r) entityPlace(r[1], r[2])))
  
  this_cont <- datex$ndef < cont.limit
  cont_pass1 <- getNdef(datex, "off1_x", "off1_y") < cont.limit
  cont_pass2 <- getNdef(datex, "off2_x", "off2_y") < cont.limit
  cont_pass3 <- getNdef(datex, "off3_x", "off3_y") < cont.limit
  cont_pass4 <- getNdef(datex, "off4_x", "off4_y") < cont.limit
  
  dat.state <- data.frame(eP_pass1, eP_pass2, eP_pass3, eP_pass4, 
                          this_cont, 
                          cont_pass1, cont_pass2, cont_pass3, cont_pass4)
  return(list(dat.fix=dat.fix, dat.spat=dat.spat,
              dat.state=dat.state))
}

evLineups <- function(datex, use.leagavg=0) {
  # finds coarsened state expected point value for each lineup combination in datex
  # inputs:
  #   datex: possession-formated data frame
  #   use.leagavg: player_id for whom we substitute league-average parameters
  teammates <- datex[,c("entity", "off1_ent", "off2_ent", "off3_ent", "off4_ent")]
  teammates <- unique(teammates)
  teammates.all <- t(apply(teammates, 1, function(r) teammatesRower(r[1], r[-1])))
  tmats <- loadAllTmats(unique(as.vector(teammates.all)))
  evs <- list()
  adder.mat <- matrix(NA, nrow=nrow(teammates), ncol=5)
  theses <- vector("list", nrow(teammates))
  for(i in 1:nrow(teammates.all)){
    # if(i %% 50 == 0)
    #   print(sprintf("%i of %i", i, nrow(teammates.all)))
    temp <- tryCatch(calcEV(tmats, teammates.all[i,], use.leagavg), error = function(e) e)
    if(!(inherits(temp, "error"))) 
      evs[[i]] <- temp 
    else 
      evs[[i]] <- NA
    these <- which(datex$entity == teammates[i, 1] & 
                     datex$off1_ent == teammates[i,2] & datex$off2_ent == teammates[i,3] & 
                     datex$off3_ent == teammates[i,4] & datex$off4_ent == teammates[i,5])
    place <- which(teammates.all[i,] == teammates[i, 1])
    theses[[i]] <- these    
    if(place == 1){
      adder.mat[i, ] <- c(14,28,42,56,0)
    } else if(place == 2){
      adder.mat[i, ] <- c(0,28,42,56,14)
    } else if(place == 3){
      adder.mat[i, ] <- c(0,14,42,56,28)
    } else if(place == 4){
      adder.mat[i, ] <- c(0,14,28,56,42)
    } else if(place == 5){
      adder.mat[i, ] <- c(0,14,28,42,56)
    } 
    
  }
  return(list(teammates=teammates, 
              teammates.all=teammates.all, evs=evs,
              theses=theses, adder.mat=adder.mat))
}

fvToEPV <- function(datex, datex.covars, fv, ev.out) {
  p.do <- epv.do <- rep(NA, nrow(datex))
  ppass1 <- epass1 <- rep(NA, nrow(datex))
  ppass2 <- epass2 <- rep(NA, nrow(datex))
  ppass3 <- epass3 <- rep(NA, nrow(datex))
  ppass4 <- epass4 <- rep(NA, nrow(datex))
  pmake <- emake <- rep(NA, nrow(datex))
  pmiss <- emiss <- rep(NA, nrow(datex))
  pTO <- eTO <- rep(NA, nrow(datex))
  enow <- rep(NA, nrow(datex))
  
  ppass1 <- fv$fvs_pass1
  ppass2 <- fv$fvs_pass2
  ppass3 <- fv$fvs_pass3
  ppass4 <- fv$fvs_pass4
  
  pmake <- fv$fvs_take*fv$fvs_make
  pmiss <- fv$fvs_take*(1-fv$fvs_make)
  pTO <- fv$fvs_TO
  
  evs <- ev.out$evs
  for(i in 1:nrow(ev.out$teammates.all)){
    these <- ev.out$theses[[i]]
    adder <- ev.out$adder.mat[i, ]
    statepass1 <- match(paste(datex.covars$dat.state$eP_pass1[these], datex.covars$dat.state$cont_pass1[these], sep="-"), state_nms)
    epass1[these] <- evs[[i]][adder[1]+statepass1]
    
    statepass2 <- match(paste(datex.covars$dat.state$eP_pass2[these], datex.covars$dat.state$cont_pass2[these], sep="-"), state_nms)
    epass2[these] <- evs[[i]][adder[2]+statepass2]
    
    statepass3 <- match(paste(datex.covars$dat.state$eP_pass3[these], datex.covars$dat.state$cont_pass3[these], sep="-"), state_nms)
    epass3[these] <- evs[[i]][adder[3]+statepass3]
    
    statepass4 <- match(paste(datex.covars$dat.state$eP_pass4[these], datex.covars$dat.state$cont_pass4[these], sep="-"), state_nms)
    epass4[these] <- evs[[i]][adder[4]+statepass4]
    
    state.now <- match(paste(datex$eP[these], datex$ndef[these] < cont.limit, sep="-"), state_nms)
    enow[these] <- evs[[i]][adder[5]+state.now]  
  }
  
  emake <- 2*(datex$threept == 0) + 3*(datex$threept == 1)
  emiss <- rep(.15, nrow(datex))
  eTO <- rep(0, nrow(datex))
  return(list(probs=data.frame(ppass1, ppass2, ppass3, ppass4, 
                               pmake, pmiss, pTO),
              vals=data.frame(epass1, epass2, epass3, epass4,
                              emake, emiss, eTO, enow)))
}

fitVals <- function(hyper, datex, datex.covars) {
  fvs <- matrix(NA, nrow=nrow(datex), ncol=length(inla.names))
  off0.idx <- match(datex$entity, hyper$player.ids)
  for (i in 1:length(inla.names)) {
    for(j in na.omit(unique(off0.idx))) {
      inds <- which(off0.idx == j)
      dat.fix <- datex.covars$dat.fix[[i]][inds, ]
      dat.spat <- datex.covars$dat.spat[[i]][inds, ]
      y <- hyper$macro.means[[j]][[i]]
      fvs[inds, i] <- as.numeric(as.vector(y[1] + dat.fix %*% y[1 + 1:ncol(dat.fix)] + dat.spat %*% tail(y, ncol(dat.spat))))
    }
  }
  return(fvs)
}

fvMatToDF <- function(fv) {
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

allCalcs <- function(datex, hyper, micro.def.mod, ev.out, nmic=50, save.positions=FALSE) {
  off.x <- as.matrix(datex[,c("x", "off1_x", "off2_x", "off3_x", "off4_x")])
  off.y <- as.matrix(datex[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")])
  off.eps.x <- rbind(rep(0, 5), off.x[-1, ] - off.x[-nrow(off.x), ])
  off.eps.y <- rbind(rep(0, 5), off.y[-1, ] - off.y[-nrow(off.y), ])
  off.eps.x[which(abs(off.eps.x) > 2, arr.ind=T)] <- 0
  off.eps.y[which(abs(off.eps.y) > 2, arr.ind=T)] <- 0
  
  def.x <- as.matrix(datex[,c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")])
  def.y <- as.matrix(datex[,c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")])
  def.eps.x <- rbind(rep(0, 5), def.x[-1, ] - def.x[-nrow(def.x), ])
  def.eps.y <- rbind(rep(0, 5), def.y[-1, ] - def.y[-nrow(def.y), ])
  def.eps.x[which(abs(def.eps.x) > 2, arr.ind=T)] <- 0
  def.eps.y[which(abs(def.eps.y) > 2, arr.ind=T)] <- 0
  
  fv.epv.list <- vector("list", nmic)
  positions <- vector("list", nmic)
  if(save.positions) {
    positions[[1]] <- list(off.x=off.x, off.y=off.y, def.x=def.x, def.y=def.y)
  }
  
  datex.covars <- datexCovars(datex)
  fv <- fvMatToDF(fitVals(hyper, datex, datex.covars))
  fv.epv.list[[1]] <- fvToEPV(datex, datex.covars, fv, ev.out)
  for(i in 2:nmic) {
    print(sprintf("micro %i of %i", i, nmic))
    micro.out <- microDatex(micro.def.mod, hyper,
                            datex, off.eps.x, off.eps.y,
                            def.eps.x, def.eps.y)
    off.eps.x <- micro.out$off.eps.x
    off.eps.y <- micro.out$off.eps.y
    def.eps.x <- micro.out$def.eps.x
    def.eps.y <- micro.out$def.eps.y
    if(save.positions) {
      off.x <- as.matrix(micro.out$datex[,c("x", "off1_x", "off2_x", "off3_x", "off4_x")])
      off.y <- as.matrix(micro.out$datex[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")])
      
      def.x <- as.matrix(micro.out$datex[,c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")])
      def.y <- as.matrix(micro.out$datex[,c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")])
      positions[[i]] <- list(off.x=off.x, off.y=off.y, def.x=def.x, def.y=def.y)
    }
    
    datex.covars <- datexCovars(micro.out$datex)
    fv <- fvMatToDF(fitVals(hyper, micro.out$datex, datex.covars))
    fv.epv.list[[i]] <- fvToEPV(micro.out$datex, datex.covars, fv, ev.out)
    datex <- micro.out$datex
  }
  return(list(fv.epv.list=fv.epv.list, positions=positions))
}

compressEPV <- function(datex, fv.epv.list) {
  nmic <- length(fv.epv.list)
  n <- nrow(fv.epv.list[[1]]$probs)
  p.here <- rep(1, n)
  ppass1 <- ppass2 <- ppass3 <- ppass4 <- pmake <- pmiss <- pTO <- rep(0, n)
  epass1 <- epass2 <- epass3 <- epass4 <- emake <- emiss <- eTO <- rep(0, n)
  for(i in 1:nmic) {
    probs <- fv.epv.list[[i]]$probs
    vals <- fv.epv.list[[i]]$vals
    p.row <- rowSums(probs)
    probs <- probs/p.row
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
  epass1 <- epass1/ppass1
  epass2 <- epass2/ppass2
  epass3 <- epass3/ppass3
  epass4 <- epass4/ppass4
  emake <- emake/pmake
  emiss <- emiss/pmiss
  eTO <- eTO/pTO
  
  probs.now <- fv.epv.list[[1]]$probs
  vals.now <- fv.epv.list[[1]]$vals
  names(probs.now) <- c("pass1", "pass2", "pass3", "pass4", "make", "miss", "TO")
  names(vals.now) <- c("pass1", "pass2", "pass3", "pass4", "make", "miss", "TO", "state")
  
  epv[which(datex$event_id %in% c(3, 4))] <- ((fv.epv.list[[1]]$vals$emake * fv.epv.list[[1]]$probs$pmake + 
                                                 fv.epv.list[[1]]$vals$emiss * fv.epv.list[[1]]$probs$pmiss) / (fv.epv.list[[1]]$probs$pmake + fv.epv.list[[1]]$probs$pmiss))[which(datex$event_id %in% c(3, 4))]
  #get pass recipient:
  pass.idx <- which(datex$event_id %in% c(22, 25))
  recips <- datex[pass.idx, c("off1_ent", "off2_ent", "off3_ent", "off4_ent")]
  pass.ent <- sapply(seq(pass.idx), function(j) match(datex$entity[pass.idx[j] + 1], recips[j, ]))
  epv[pass.idx[which(pass.ent == 1)]] <- fv.epv.list[[1]]$vals$epass1[pass.idx[which(pass.ent == 1)]]
  epv[pass.idx[which(pass.ent == 2)]] <- fv.epv.list[[1]]$vals$epass2[pass.idx[which(pass.ent == 2)]]
  epv[pass.idx[which(pass.ent == 3)]] <- fv.epv.list[[1]]$vals$epass3[pass.idx[which(pass.ent == 3)]]
  epv[pass.idx[which(pass.ent == 4)]] <- fv.epv.list[[1]]$vals$epass4[pass.idx[which(pass.ent == 4)]]
  epv[which(datex$event_id == 7)] <- fv.epv.list[[1]]$vals$eTO[which(datex$event_id == 7)]
  return(list(epv=epv, probs=data.frame(pass1=ppass1, pass2=ppass2, pass3=ppass3, pass4=ppass4,
                                        make=pmake, miss=pmiss, TO=pTO, other=p.here),
              vals=data.frame(pass1=epass1, pass2=epass2, pass3=epass3, pass4=epass4,
                              make=emake, miss=emiss, TO=eTO, other=eother),
              probs.now=probs.now,
              vals.now=vals.now))
}

microDatex <- function(micro.def.mod, hyper, 
                       datex, 
                       off.eps.x, off.eps.y, 
                       def.eps.x, def.eps.y){
  off.x <- as.matrix(datex[,c("x", "off1_x", "off2_x", "off3_x", "off4_x")])
  off.y <- as.matrix(datex[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")])
  
  off0.idx <- match(datex$entity, hyper$player.ids)
  off1.idx <- match(datex$off1_ent, hyper$player.ids)
  off2.idx <- match(datex$off2_ent, hyper$player.ids)
  off3.idx <- match(datex$off3_ent, hyper$player.ids)
  off4.idx <- match(datex$off4_ent, hyper$player.ids)
  
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
  new.x.mu.0 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$x, datex$y)))$proj$A * t(hyper$spats.x.with[, off0.idx]))
  new.x.mu.1 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$off1_x, datex$off1_y)))$proj$A * t(hyper$spats.x.without[, off1.idx]))
  new.x.mu.2 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$off2_x, datex$off2_y)))$proj$A * t(hyper$spats.x.without[, off2.idx]))
  new.x.mu.3 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$off3_x, datex$off3_y)))$proj$A * t(hyper$spats.x.without[, off3.idx]))
  new.x.mu.4 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$off4_x, datex$off4_y)))$proj$A * t(hyper$spats.x.without[, off4.idx]))
  off.new.x.0 <- intercepts.x.0 + datex$x + new.x.mu.0 + coefs.x.0 * off.eps.x[, 1] + rnorm(length(sigmas.x.0), 0, sigmas.x.0)
  off.new.x.1 <- intercepts.x.1 + datex$off1_x + new.x.mu.1 + coefs.x.1 * off.eps.x[, 2] + rnorm(length(sigmas.x.1), 0, sigmas.x.1)
  off.new.x.2 <- intercepts.x.2 + datex$off2_x + new.x.mu.2 + coefs.x.2 * off.eps.x[, 3] + rnorm(length(sigmas.x.2), 0, sigmas.x.2)
  off.new.x.3 <- intercepts.x.3 + datex$off3_x + new.x.mu.3 + coefs.x.3 * off.eps.x[, 4] + rnorm(length(sigmas.x.3), 0, sigmas.x.3)
  off.new.x.4 <- intercepts.x.4 + datex$off4_x + new.x.mu.4 + coefs.x.4 * off.eps.x[, 5] + rnorm(length(sigmas.x.4), 0, sigmas.x.4)
  
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
  new.y.mu.0 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$x, datex$y)))$proj$A * t(hyper$spats.y.with[, off0.idx]))
  new.y.mu.1 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$off1_x, datex$off1_y)))$proj$A * t(hyper$spats.y.without[, off1.idx]))
  new.y.mu.2 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$off2_x, datex$off2_y)))$proj$A * t(hyper$spats.y.without[, off2.idx]))
  new.y.mu.3 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$off3_x, datex$off3_y)))$proj$A * t(hyper$spats.y.without[, off3.idx]))
  new.y.mu.4 <- rowSums(inla.mesh.projector(mesh, loc=as.matrix(cbind(datex$off4_x, datex$off4_y)))$proj$A * t(hyper$spats.y.without[, off4.idx]))
  off.new.y.0 <- intercepts.y.0 + datex$y + new.y.mu.0 + coefs.y.0 * off.eps.y[, 1] + rnorm(length(sigmas.y.0), 0, sigmas.y.0)
  off.new.y.1 <- intercepts.y.1 + datex$off1_y + new.y.mu.1 + coefs.y.1 * off.eps.y[, 2] + rnorm(length(sigmas.y.1), 0, sigmas.y.1)
  off.new.y.2 <- intercepts.y.2 + datex$off2_y + new.y.mu.2 + coefs.y.2 * off.eps.y[, 3] + rnorm(length(sigmas.y.2), 0, sigmas.y.2)
  off.new.y.3 <- intercepts.y.3 + datex$off3_y + new.y.mu.3 + coefs.y.3 * off.eps.y[, 4] + rnorm(length(sigmas.y.3), 0, sigmas.y.3)
  off.new.y.4 <- intercepts.y.4 + datex$off4_y + new.y.mu.4 + coefs.y.4 * off.eps.y[, 5] + rnorm(length(sigmas.y.4), 0, sigmas.y.4)
  
  mod.x <- micro.def.mod$mod.x
  mod.y <- micro.def.mod$mod.y
  
  opts.x <- datex$ball_x*.27 + off.x*.62 + 4.75*.11
  opts.y <- datex$ball_y*.27 + off.y*.62 + 25*.11
  def.x <- as.matrix(datex[,c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")])
  def.y <- as.matrix(datex[,c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")])
  
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
  
  off.new.x <- cbind(off.new.x.0, off.new.x.1, off.new.x.2, off.new.x.3, off.new.x.4)
  off.new.y <- cbind(off.new.y.0, off.new.y.1, off.new.y.2, off.new.y.3, off.new.y.4)
  
  new.opts.x <- datex$ball_x*.27 + off.new.x*.62 + 4.75*.11
  new.opts.y <- datex$ball_y*.27 + off.new.y*.62 + 25*.11
  
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
  
  ball.new.x <- datex$ball_x + off.eps.x[, 1]
  ball.new.y <- datex$ball_y + off.eps.y[, 1]
  
  
  datex[,c("x", "off1_x", "off2_x", "off3_x", "off4_x")] <- off.new.x
  datex[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")] <- off.new.y
  datex[,c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")] <- def.new.x
  datex[,c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")] <- def.new.y
  datex[, "ball_x"] <- ball.new.x
  datex[, "ball_y"] <- ball.new.y
  datex[,"ndef"] <- getNdef(datex, "x", "y")
  datex[, "eP"] <- getPlace(datex) 
  threept.out <- getThree(datex)
  datex[, "bask.dist"] <- threept.out$bask.dist
  datex[, "threept"] <- threept.out$threept
  return(list(datex=datex, off.eps.x=off.eps.x, off.eps.y=off.eps.y, 
              def.eps.x=def.eps.x, def.eps.y=def.eps.y))
}

teammatesRower <- function(id, teamrow){
  lv.me <- as.numeric(latent_var[match(id, players$player_id)]) + id/1000000
  lv.team <- as.numeric(latent_var[match(teamrow, players$player_id)]) + as.numeric(teamrow)/1000000
  newrow.ind <- order(c(lv.me, lv.team))
  return(c(id, as.numeric(teamrow))[newrow.ind])  
}

loadAllTmats <- function(ids) {
  tmats <- vector("list", length(ids))
  for(i in 1:length(ids)) {
    # print(sprintf("loading tmat %i of %i", i, length(ids)))
    if(is.na(ids[i])) {
      tmats[[i]] <- list(id=NA, tmat.ind=NULL, tmat.pos=NULL)
      next
    }
    id <- ids[i]
    if(file.exists(sprintf("%s/tmats/%s.Rdata", data.dir, id))) {
      load(sprintf("%s/tmats/%s.Rdata", data.dir, id))
      tmats[[i]] <- list(id=id, tmat.ind=tmat.ind, tmat.pos=tmat.pos)
    } else {
      tmats[[i]] <- list(id=id, tmat.ind=NULL, tmat.pos=NULL)
    }
  }
  return(tmats)
}

tmatRowMaker <- function(tmats, id, place, use.leagavg=0){
  tmat.ids <- sapply(tmats, function(obj) obj$id)
  idx <- which(tmat.ids == id)
  if(length(idx) == 0)
    idx <- which(is.na(tmat.ids))
  #if(length(idx) != 1) {
  #  print("error: id matched more than once")
  #  return()
  #}
  if(is.null(tmats[[idx]]$tmat.ind)) {
    pos <- as.vector(players$position[which(players$player_id == id)])
    these <- which(players$position == pos)
    id <- players$player_id[these[1]]
    load(sprintf("%s/tmats/%s.Rdata", data.dir, id))
    mat <- tmat.pos
  } else {
    if (id %in% use.leagavg) {
      mat <- tmats[[idx]]$tmat.pos
    } else {
      mat <- tmats[[idx]]$tmat.ind
    }
  }
  if(place == 1){
    tmat <- cbind(
      mat$micros, mat$passes1, mat$passes2,
      mat$passes3, mat$passes4, mat$absorbs)
  } else if(place == 2){
    tmat <- cbind(
      mat$passes1, mat$micros, mat$passes2,
      mat$passes3, mat$passes4, mat$absorbs)
  } else if(place == 3){
    tmat <- cbind(
      mat$passes1, mat$passes2, mat$micros, 
      mat$passes3, mat$passes4, mat$absorbs)
  } else if(place == 4){
    tmat <- cbind(
      mat$passes1, mat$passes2,
      mat$passes3, mat$micros, mat$passes4, mat$absorbs)
  } else {
    tmat <- cbind(
      mat$passes1, mat$passes2,
      mat$passes3, mat$passes4, mat$micros, mat$absorbs)
  } 
  return(tmat)
}

calcEV <- function(tmats, teamrow, use.leagavg=0){
  p1_mat <- tmatRowMaker(tmats, teamrow[1], 1, use.leagavg)
  p2_mat <- tmatRowMaker(tmats, teamrow[2], 2, use.leagavg)
  p3_mat <- tmatRowMaker(tmats, teamrow[3], 3, use.leagavg)
  p4_mat <- tmatRowMaker(tmats, teamrow[4], 4, use.leagavg)
  p5_mat <- tmatRowMaker(tmats, teamrow[5], 5, use.leagavg)
  tmat <- rbind(p1_mat, p2_mat, p3_mat, p4_mat, p5_mat)
  tmat[which(is.na(tmat), arr.ind=T)] <- 0
  tmat <- as.matrix(tmat)
  bad.rows <- which(rowSums(tmat) == 0) # never visited
  tmat[cbind(bad.rows, rep(ncol(tmat), length(bad.rows)))] <- 1
  tmat[cbind(bad.rows, bad.rows)] <- 1
  tmat <- tmat/rowSums(tmat, na.rm=T)
  tmatTrans <- tmat[1:nrow(tmat), 1:nrow(tmat)]
  tmatAbsorb <- tmat[, -(1:nrow(tmat))]
  absProb <- solve(diag(1, nrow(tmat)) - tmatTrans)%*%tmatAbsorb 
  absValues <- c(2,.15,3,.15,0)
  vals <- as.numeric(absProb%*%absValues)
  vals[bad.rows] <- NA
  return(vals)
}

doEPV <- function(tmats, use.leagavg=0, fullres=TRUE){
  evs <- list()
  for(i in 1:nrow(teammates.all)){
    # if(i %% 50 == 0)
    #   print(sprintf("%i of %i", i, nrow(teammates.all)))
    temp <- tryCatch(calcEV(tmats, teammates.all[i,], use.leagavg), error = function(e) e)
    if(!(inherits(temp, "error"))) evs[[i]] <- temp else evs[[i]] <- NA
  }
  
  p.do <- epv.do <- rep(NA, nrow(datex))
  ppass1 <- epass1 <- rep(NA, nrow(datex))
  ppass2 <- epass2 <- rep(NA, nrow(datex))
  ppass3 <- epass3 <- rep(NA, nrow(datex))
  ppass4 <- epass4 <- rep(NA, nrow(datex))
  pmake <- emake <- rep(NA, nrow(datex))
  pmiss <- emiss <- rep(NA, nrow(datex))
  pTO <- eTO <- rep(NA, nrow(datex))
  
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
    these <- which(datex$off1_ent == teammates[i,1] & datex$off2_ent == teammates[i,2] & 
                     datex$off3_ent == teammates[i,3] & datex$off4_ent == teammates[i,4])
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
      statein <- match(paste(datex$eP[these], datex$ndef[these] < cont.limit, sep="-"), state_nms)
      epv.do[these] <- evs[[i]][14*(place-1) + statein]
    }
  }
  
  emake <- 2*(datex$threept == 0) + 3*(datex$threept == 1)
  emiss <- rep(.15, nrow(datex))
  eTO <- rep(0, nrow(datex))
  
  p.do <- ppass1 + ppass2 + ppass3 + ppass4 + pmake + pmiss + pTO
  
  if(fullres){
    epv.do <- (ppass1*epass1 + ppass2*epass2 + ppass3*epass3 + ppass4*epass4 + pmake*emake + pmiss*emiss + pTO*eTO)/p.do
  }
  
  epv.do[which(datex$event_id %in% c(3,4))] <- ((pmake*emake + pmiss*emiss)/(pmake + pmiss))[which(datex$event_id %in% c(3,4))]
  epv.do[which(datex$event_id == 7)] <- 0
  epv.do[which(datex$event_id == 31)] <- epass1[which(datex$event_id == 31)]
  epv.do[which(datex$event_id == 32)] <- epass2[which(datex$event_id == 32)]
  epv.do[which(datex$event_id == 33)] <- epass3[which(datex$event_id == 33)]
  epv.do[which(datex$event_id == 34)] <- epass4[which(datex$event_id == 34)]
  
  filter.epv <- epv.do
  thesediff <- diff(datex$time)
  ends <- c(which(thesediff < 0 | thesediff > 100), nrow(datex))
  #ends <- ends[which(substr(datex$game[ends],7,8) %in% c("01","02","03"))]
  starts <- c(1, which(thesediff < 0 | thesediff > 100) + 1)
  #starts <- starts[which(substr(datex$game[ends],7,8) %in% c("01","02","03"))]
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
  shot_sat <- mean((epv.do - epass)[which(datex$event_id %in% c(3,4))], na.rm=T)
  
  mm <- sapply(1:nrow(ppp), function(i) sum(pmake[i] + pmiss[i] > ppp[i,], na.rm=T))
  eshot <- (pmake*emake + pmiss*emiss)/(pmake + pmiss)
  pass_sat <- mean((epv.do - eshot)[which(mm >= 3 & datex$event_id %in% c(31,32,33,34))], na.rm=T)
  
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
  epv.mic[which(datex$event_id %in% c(3,4,7,31,32,33,34))] <- epv.do[which(datex$event_id %in% c(3,4,7,31,32,33,34))]
  
  
  return(data.frame(game=datex$game, quarter=datex$quarter, time=datex$time, game_clock=datex$game_clock,
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
