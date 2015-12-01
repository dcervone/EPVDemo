# **********************
#
# Functions that load, manipulate, or generate parameter
# inferences for multiresolution transition models
#
# **********************

microDefModel <- function(tdat) {
  # Parameter inference for defensive microtransition model
  # Inputs
  #   tdat: 'tdat' formatted data frame
  # Outputs
  #   mod.x, mod.y: lm objects with parameters for defensive microtransitions
  
  # matrices of player-positions
  off.x <- as.matrix(tdat[, c("x", "off1_x", "off2_x", "off3_x", "off4_x")])
  off.y <- as.matrix(tdat[,c("y", "off1_y", "off2_y", "off3_y", "off4_y")])
  def.x <- as.matrix(tdat[, c("def1_x", "def2_x", "def3_x", "def4_x", "def5_x")])
  def.y <- as.matrix(tdat[, c("def1_y", "def2_y", "def3_y", "def4_y", "def5_y")])
  
  # optimal defender points (from Franks et al 2014)
  opts.x <- tdat$ball_x * .27 + off.x * .62 + 4.75 * .11
  opts.y <- tdat$ball_y * .27 + off.y * .62 + 25 * .11
  
  # distances of each defender from optimal points
  wts <- array(NA, dim=c(nrow(off.x), 5, 5))
  for(i in 1:nrow(off.x)){
    def.x.n <- def.x[i, ]
    def.y.n <- def.y[i, ]
    guarding <- sapply(1:5, function(j) {
      wts <- 1 / ((opts.x[i, ] - def.x.n[j])^2 + (opts.y[i, ] - def.y.n[j])^2)
      wts <- wts / sum(wts)
      return(wts)
    })
    wts[i, , ] <- guarding
  }
  
  # defensive "mean" components based on matchups
  means.x <- t(sapply(1:nrow(off.x), function(j) colSums(opts.x[j, ] * wts[j, , ])))
  means.y <- t(sapply(1:nrow(off.y), function(j) colSums(opts.y[j, ] * wts[j, , ])))
  
  # differenced positions for AR model
  def.eps.x <- as.vector(diff(def.x))
  residual.x <- as.vector((def.x - means.x)[-1, ])
  opt.eps.x <- as.vector(diff(means.x))
  bads <- which(abs(def.eps.x) > .5)
  residual.x <- residual.x[-bads]
  def.eps.x <- def.eps.x[-bads]
  opt.eps.x <- opt.eps.x[-bads]
  
  # x component for defensive microtransitions
  mod.x <- lm(def.eps.x[-1] ~ def.eps.x[-length(def.eps.x)] + residual.x[-length(residual.x)] + opt.eps.x[-length(opt.eps.x)])
  
  # repeat above for y component
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

getHyperParams <- function(tdat) {
  # Loads pre-computed multiresolution transition model parameters 
  # for all players in tdat
  # Inputs
  #   tdat: 'tdat' formatted data frame
  # Outputs
  #   list of player-specific model parameters
  
  # load macro entry model output
  inlas <- vector("list", length(inla.names))
  for(i in 1:length(inla.names)) {
    # print(sprintf("loading inla output %i of %i", i, length(inla.names)))
    load(sprintf("%s/INLA_%s.Rdata", data.dir, inla.names[i]))
    inlas[[i]] <- inla.out
  }
  
  # which players in the league are in tdat
  player.ids <- unique(unlist(tdat[, c("entity", sapply(1:4, function(q) sprintf("off%i_ent", q)),
                                        sapply(1:5, function(q) sprintf("def%i_ent", q)))]))
  player.ids <- player.ids[which(player.ids > 0)]
  
  # for these players, organize micro/macro parameter estimates
  micro.inlas <- vector("list",  length(player.ids))
  macro.means <- vector("list",  length(player.ids))
  for(i in 1:length(player.ids)) {
    pl <- player.ids[i]
    # micro parameters
    load(sprintf("%s/micros/%s.Rdata", data.dir, pl))
    micro.inlas[[i]] <- list(with=with.ball, without=without.ball)
    macro.means[[i]] <- vector("list", length(inla.names))
    for(j in 1:length(inla.names)) {
      # macro parameters
      id.n <- match(pl, players$player_id)
      # fixed effects
      y.fix <- inlas[[j]]$summary.fixed[, "mean"]
      temp <- names(inlas[[j]]$summary.random)
      basis.inds <- c(which(temp == "p.int"), grep("p.b[0-9][0-9]*", temp))
      cov.inds <- setdiff(seq(length(inlas[[j]]$summary.random)), basis.inds)
      # random effects
      y.rand <- c(inlas[[j]]$summary.random$p.int[id.n, "mean"], 
                  sapply(cov.inds, 
                         function(k) inlas[[j]]$summary.random[[k]][id.n, "mean"]),
                  inlas[[j]]$summary.random$p.b1[id.n + nrow(players) * (1:n.basis), "mean"])
      macro.means[[i]][[j]] <- y.fix + y.rand
    }
  }
  
  # re-ogranizing/formatting parameters for output
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

evLineups <- function(tdat, use.leagavg=0) {
  # Finds coarsened state expected point value for each lineup 
  # combination in tdat
  # Inputs:
  #   tdat: 'tdat' formated data frame
  #   use.leagavg: player_id for whom we substitute league-average parameters
  # Outputs:
  #   teammates: unordered lineup matrix
  #   teammates.all: position-ordered lineup matrix
  #   evs: list of coarsened state expected values for each lineup
  #   lineup.inds: which points in tdat correpond to which lineups
  #   adder.mat: identifies ballcarrier's rows' position in transition matrix

  # matrix of offensive lineup combinations
  teammates <- tdat[,c("entity", "off1_ent", "off2_ent", "off3_ent", "off4_ent")]
  teammates <- unique(teammates)
  teammates.all <- t(apply(teammates, 1, function(r) teammatesRower(r[1], r[-1])))
  
  # loads coarsened state transition matrices for all offensive players
  tmats <- loadAllTmats(unique(as.vector(teammates.all)))
  
  # for each lineup, organize rows/columns of transition matrix
  evs <- list()
  adder.mat <- matrix(NA, nrow=nrow(teammates), ncol=5)
  lineup.inds <- vector("list", nrow(teammates))
  for(i in 1:nrow(teammates.all)){
    # if(i %% 50 == 0)
    #   print(sprintf("%i of %i", i, nrow(teammates.all)))
    temp <- tryCatch(calcEV(tmats, teammates.all[i,], use.leagavg), error = function(e) e)
    if(!(inherits(temp, "error"))) 
      evs[[i]] <- temp 
    else 
      evs[[i]] <- NA
    these <- which(tdat$entity == teammates[i, 1] & 
                     tdat$off1_ent == teammates[i,2] & tdat$off2_ent == teammates[i,3] & 
                     tdat$off3_ent == teammates[i,4] & tdat$off4_ent == teammates[i,5])
    place <- which(teammates.all[i,] == teammates[i, 1])
    lineup.inds[[i]] <- these    
    n.state <- length(state_nms)
    if(place == 1){
      adder.mat[i, ] <- c(n.state * 1:4, 0)
    } else if(place == 2){
      adder.mat[i, ] <- c(0, n.state * c(2,3,4,1))
    } else if(place == 3){
      adder.mat[i, ] <- c(0, n.state * c(1,3,4,2))
    } else if(place == 4){
      adder.mat[i, ] <- c(0, n.state * c(1,2,4,3))
    } else if(place == 5){
      adder.mat[i, ] <- c(0, n.state * 1:4)
    } 
    
  }
  return(list(teammates=teammates, 
              teammates.all=teammates.all, evs=evs,
              lineup.inds=lineup.inds, adder.mat=adder.mat))
}

loadAllTmats <- function(ids) {
  # loads coarsened state transition matrix components for players
  # Inputs:
  #   ids: vector of player IDs for whom to load components
  # Outputs:
  #   tmats: list of transtition matrix components by player
  
  tmats <- vector("list", length(ids))
  for(i in 1:length(ids)) {
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
  # organizes transition matrices based on ballcarrier and lineup
  # Inputs:
  #   tmats: list of players' trans. mat. components
  #   id: ballcarrier
  #   place: position in the lineup of player 'id'
  #   use.leagavg: player id for which we substitute league average probabilities
  # Outputs:
  #   tmat: coarsened state probability transition matrix
  
  tmat.ids <- sapply(tmats, function(obj) obj$id)
  idx <- which(tmat.ids == id)
  if(length(idx) == 0)
    idx <- which(is.na(tmat.ids))
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
  # order columns based on lineup
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

teammatesRower <- function(id, teamrow){
  # makes a position-ordered lineup out of ballcarrier and teammates
  # Inputs:
  #   id: ballcarrier id
  #   teamrow: vector of teammate player ids
  # Outputs:
  #   id, position-ordered teammates
  
  lv.me <- as.numeric(players$position.number[match(id, players$player_id)]) + id/1000000
  lv.team <- as.numeric(players$position.number[match(teamrow, players$player_id)]) + as.numeric(teamrow)/1000000
  newrow.ind <- order(c(lv.me, lv.team))
  return(c(id, as.numeric(teamrow))[newrow.ind])  
}

calcEV <- function(tmats, teamrow, use.leagavg=0){
  # Calculates coarsened state expected point values given ballcarrier and
  # transition probability matrix
  # Inputs:
  #   tmats: list of transition probability matrices
  #   teamrow: ballcarrier id + ordered teammates
  #   use.leagavg: player id for whom we substitute league average probabilities
  # Outputs:
  #   vals: list of coarsened state expected
  
  # row groups for each player in lineup
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
  
  # divide into absorbing states
  tmatTrans <- tmat[1:nrow(tmat), 1:nrow(tmat)]
  tmatAbsorb <- tmat[, -(1:nrow(tmat))]
  absProb <- solve(diag(1, nrow(tmat)) - tmatTrans) %*% tmatAbsorb 
  absValues <- c(2,.15,3,.15,0)
  vals <- as.numeric(absProb%*%absValues)
  vals[bad.rows] <- NA
  return(vals)
}

