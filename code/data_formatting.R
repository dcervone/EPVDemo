# **********************
#
# Format moments data into offensive possession halfcourt data, better structure for modeling
#
# **********************

# mark events from list "ev" indicating ball possession
poss.foo <- function(dat, poss, ev, next0=FALSE) {
  for(i in 1:5) {
    these <- which(dat[, grep(sprintf("h%i_event", i), names(dat))] %in% ev)
    poss[these] <- as.vector(dat[these, grep(sprintf("h%i_ent", i), names(dat))])
    if(next0)
      poss[these + 1] <- 0
  }
  for(i in 1:5) {
    these <- which(dat[, grep(sprintf("a%i_event", i), names(dat))] %in% ev)
    poss[these] <- as.vector(dat[these, grep(sprintf("a%i_ent", i), names(dat))])
    if(next0)
      poss[these + 1] <- 0
  }
  return(poss)  
}

# mark events from list "ev" indicating no ball possession
zero.foo <- function(dat, poss, ev, next0=FALSE) {
  for(i in 1:5) {
    these <- which(dat[, grep(sprintf("h%i_event", i), names(dat))] %in% ev)
    poss[these + 1] <- 0
  }
  for(i in 1:5) {
    these <- which(dat[, grep(sprintf("a%i_event", i), names(dat))] %in% ev)
    poss[these + 1] <- 0
  }
  return(poss)  
}

# create variable tracking which entity is ballcarrier
possession.indicator <- function(dat) {
  # eliminate 25 and 3 simultaneous tags
  mult.events <- which(rowSums(!is.na(dat[, grep("event", names(dat))])) > 1)
  if(length(mult.events) > 0) {
    assists <- which(dat[mult.events, grep("event", names(dat))] == 25, arr.ind=T)
    if(nrow(assists) > 0)
      assists[, 1] <- mult.events[assists[, 1]]
    dat[assists] <- NA
  }

  # fill in possession indicator at event points
  poss <- rep(NA, nrow(dat))
  poss <- poss.foo(dat, poss, c(22, 25), TRUE) # passes
  poss <- poss.foo(dat, poss, 7, TRUE) # turnovers
  poss <- poss.foo(dat, poss, c(3, 4), TRUE) # shots
  poss <- poss.foo(dat, poss, 23, FALSE) # possession
  poss <- poss.foo(dat, poss, 21, FALSE) # dribbles
  poss <- poss.foo(dat, poss, c(5, 6), FALSE) # rebounds
  poss <- zero.foo(dat, poss, c(1, 2, 8, 9, 11, 13, 15, 16, 17, 18, 19, 20))

  #diagnose missing pass indicator
  event.idx <- which(!(is.na(poss)))
  idx.end <- event.idx[-1]
  idx.start <- event.idx[-length(event.idx)]
  event.start <- poss[event.idx]
  event.end <- event.start[-1]
  event.start <- event.start[-length(event.start)]
  test <- which(event.start != 0 & 
                  event.start != event.end & 
                  event.end != 0 & 
                  idx.end - idx.start > 1)
  if(length(test) > 0) {
    # missing pass event
    for(i in 1:length(test)) {
      inds <- idx.start[test[i]]:idx.end[test[i]]
      temp <- dat[idx.start[test[i]], grep("event", names(dat))]
      pos.pl <- substr(names(temp)[which(!is.na(temp))], 1, 2)
      pos.x <- as.numeric(as.vector(dat[inds, paste0(pos.pl, "_x")]))
      pos.y <- as.numeric(as.vector(dat[inds, paste0(pos.pl, "_y")]))
      ball.x <- as.numeric(as.vector(dat[inds, "x"]))
      ball.y <- as.numeric(as.vector(dat[inds, "y"]))
      dist <- sqrt((pos.x - ball.x)^2 + (pos.y - ball.y)^2)
      cutoff <- which(dist > 4)
      if(length(cutoff) > 0) {
        ind <- min(cutoff)
        if(ind == 1)
          ind <- 2
        if(ind == length(dist))
          ind <- length(dist) - 1
      } else {
        ind <- length(dist) - 1
      }
      poss[inds[ind]] <- 0
      dat[inds[ind], paste0(pos.pl, "_event")] <- 22
    }
  }

  # fill in ballcarrier
  poss[1] <- 0
  for(i in 2:length(poss)) {
    if(is.na(poss[i]))
      poss[i] <- poss[i-1]
  }
  poss <- poss[1:nrow(dat)]

  # checking possession integrity
  pos.x <- dat$x
  pos.y <- dat$y
  print("checking possession integrity")
  for(i in 1:length(poss)) {
    if(i %% 1000 == 0)
      print(sprintf("row %i of %i", i, length(poss)))
    if(poss[i] != 0) {
      temp <- which(dat[i, ] == poss[i])
      if(length(temp) == 1) {
        pos.pl <- substr(names(dat)[temp], 1, 2)
        pos.x[i] <- as.numeric(as.vector(dat[i, paste0(pos.pl, "_x")]))
        pos.y[i] <- as.numeric(as.vector(dat[i, paste0(pos.pl, "_y")]))
      } else {
        pos.x[i] <- NA
        pos.y[i] <- NA
      }
    }
  }
  pos.dist <- sqrt((pos.x - dat$x)^2 + (pos.y - dat$y)^2)
  poss.logic <- pos.dist < 6
  print(sprintf("check for ball and ballcarrier distance within 6 ft: %.3f", mean(poss.logic, na.rm=T)))
  poss[which(!(poss.logic))] <- 0
  gc()
  return(poss)
}

# column for ballcarrier identity, position, and those of his teammates
rearrange.data <- function(dat, poss) {
  pls <- unique(poss)
  pls <- pls[which(pls > 0)]
  entity <- poss
  x <- rep(NA, nrow(dat))
  y <- rep(NA, nrow(dat))
  event_id <- rep(NA, nrow(dat))
  off.ent.mat <- matrix(NA, nrow=nrow(dat), ncol=4)
  off.x.mat <- matrix(NA, nrow=nrow(dat), ncol=4)
  off.y.mat <- matrix(NA, nrow=nrow(dat), ncol=4)
  off.event.mat <- matrix(NA, nrow=nrow(dat), ncol=4)
  def.ent.mat <- matrix(NA, nrow=nrow(dat), ncol=5)
  def.x.mat <- matrix(NA, nrow=nrow(dat), ncol=5)
  def.y.mat <- matrix(NA, nrow=nrow(dat), ncol=5)
  def.event.mat <- matrix(NA, nrow=nrow(dat), ncol=5)
  team <- rep(NA, nrow(dat))
  
  for(i in 1:length(pls)) {
    print(sprintf("player %i of %i", i, length(pls)))
    pl.inds <- which(poss == pls[i])
    pl.team <- substr(names(dat)[which(dat[pl.inds[1], ] == pls[i])], 1, 1)
    pl.num <- sapply(pl.inds, function(j) {
      temp <- which(dat[j, ] == pls[i])
      if(length(temp) == 1)
        temp
      else
        NA
    })
    pl.num <- substr(names(dat)[pl.num], 2, 2)
    team[pl.inds] <- pl.team
    x[pl.inds] <- dat[cbind(pl.inds, 
                                match(paste0(pl.team, paste0(pl.num, "_x")), names(dat)))]
    y[pl.inds] <- dat[cbind(pl.inds, 
                                match(paste0(pl.team, paste0(pl.num, "_y")), names(dat)))]
    event_id[pl.inds] <- dat[cbind(pl.inds, 
                                       match(paste0(pl.team, paste0(pl.num, "_event")), names(dat)))]
    
    for(k in 1:4) {
      off.ent.mat[pl.inds, k] <- dat[cbind(pl.inds, 
                                               match(paste0(pl.team, paste0(as.vector(sapply(pl.num, function(x) setdiff(1:5, x)[k])), 
                                                                            "_ent")), names(dat)))]
      off.x.mat[pl.inds, k] <- dat[cbind(pl.inds, 
                                             match(paste0(pl.team, paste0(as.vector(sapply(pl.num, function(x) setdiff(1:5, x)[k])), 
                                                                          "_x")), names(dat)))]
      off.y.mat[pl.inds, k] <- dat[cbind(pl.inds, 
                                             match(paste0(pl.team, paste0(as.vector(sapply(pl.num, function(x) setdiff(1:5, x)[k])), 
                                                                          "_y")), names(dat)))]
      off.event.mat[pl.inds, k] <- dat[cbind(pl.inds, 
                                                 match(paste0(pl.team, paste0(as.vector(sapply(pl.num, function(x) setdiff(1:5, x)[k])), 
                                                                              "_event")), names(dat)))]
    }
    op.team <- setdiff(c("a", "h"), pl.team)
    def.ent.mat[pl.inds, ] <- cbind(dat[pl.inds, paste0(op.team, "1_ent")],
                                    dat[pl.inds, paste0(op.team, "2_ent")],
                                    dat[pl.inds, paste0(op.team, "3_ent")],
                                    dat[pl.inds, paste0(op.team, "4_ent")],
                                    dat[pl.inds, paste0(op.team, "5_ent")])
    def.x.mat[pl.inds, ] <- cbind(dat[pl.inds, paste0(op.team, "1_x")],
                                  dat[pl.inds, paste0(op.team, "2_x")],
                                  dat[pl.inds, paste0(op.team, "3_x")],
                                  dat[pl.inds, paste0(op.team, "4_x")],
                                  dat[pl.inds, paste0(op.team, "5_x")])
    def.y.mat[pl.inds, ] <- cbind(dat[pl.inds, paste0(op.team, "1_y")],
                                  dat[pl.inds, paste0(op.team, "2_y")],
                                  dat[pl.inds, paste0(op.team, "3_y")],
                                  dat[pl.inds, paste0(op.team, "4_y")],
                                  dat[pl.inds, paste0(op.team, "5_y")])
    def.event.mat[pl.inds, ] <- cbind(dat[pl.inds, paste0(op.team, "1_event")],
                                      dat[pl.inds, paste0(op.team, "2_event")],
                                      dat[pl.inds, paste0(op.team, "3_event")],
                                      dat[pl.inds, paste0(op.team, "4_event")],
                                      dat[pl.inds, paste0(op.team, "5_event")])
  }
  tdat <- data.frame(dat[, c("game", "time", "quarter", "shot_clock", "game_clock",
                                      "x", "y", "z")],
                          entity, team, x, y, event_id, 
                          off.ent.mat, off.x.mat, off.y.mat, off.event.mat,
                          def.ent.mat, def.x.mat, def.y.mat, def.event.mat)
  names(tdat) <- c("game", "time", "quarter", "shot_clock", "game_clock",
                        "ball_x", "ball_y", "ball_z", "entity", "team",
                        "x", "y", "event_id", 
                        sapply(1:4, function(v) sprintf("off%s_ent", v)),
                        sapply(1:4, function(v) sprintf("off%s_x", v)),
                        sapply(1:4, function(v) sprintf("off%s_y", v)),
                        sapply(1:4, function(v) sprintf("off%s_event", v)),
                        sapply(1:5, function(v) sprintf("def%s_ent", v)),
                        sapply(1:5, function(v) sprintf("def%s_x", v)),
                        sapply(1:5, function(v) sprintf("def%s_y", v)),
                        sapply(1:5, function(v) sprintf("def%s_event", v)))

  for(i in 1:4) {
    tdat[, sprintf("off%i_x", i)] <- as.numeric(as.vector(tdat[, sprintf("off%i_x", i)]))
    tdat[, sprintf("off%i_y", i)] <- as.numeric(as.vector(tdat[, sprintf("off%i_y", i)]))
  }
  for(i in 1:5) {
    tdat[, sprintf("def%i_x", i)] <- as.numeric(as.vector(tdat[, sprintf("def%i_x", i)]))
    tdat[, sprintf("def%i_y", i)] <- as.numeric(as.vector(tdat[, sprintf("def%i_y", i)]))
  }
  tdat[, "x"] <- as.numeric(as.vector(tdat[, "x"]))
  tdat[, "y"] <- as.numeric(as.vector(tdat[, "y"]))
  tdat[, "ball_x"] <- as.numeric(as.vector(tdat[, "ball_x"]))
  tdat[, "ball_y"] <- as.numeric(as.vector(tdat[, "ball_y"]))
  tdat[, "ball_z"] <- as.numeric(as.vector(tdat[, "ball_z"]))
  rm(off.ent.mat, off.x.mat, off.y.mat, off.event.mat,
     def.ent.mat, def.x.mat, def.y.mat, def.event.mat)
  gc()
  return(tdat)
}

# flip data to offensive halfcourt
offensive.halfcourt <- function(dat) {
  print("flipping")
  homeQ1x <- dat$x[which(dat$team == "h" & dat$quarter == 1 & dat$event_id %in% c(3, 4))]
  cols.to.flip.x <- c(match("x", names(dat)), 
                      grep("_x", names(dat)))
  cols.to.flip.y <- c(match("y", names(dat)), 
                      grep("_y", names(dat)))
  if(mean(homeQ1x) < 47) {
    rows.to.flip.h <- which(dat$team == "h" & dat$quarter > 2 & dat$quarter <= 4)
    rows.to.flip.a <- which(dat$team == "a" & dat$quarter <= 2)
  } else {
    rows.to.flip.h <- which(dat$team == "h" & dat$quarter <= 2)
    rows.to.flip.a <- which(dat$team == "a" & dat$quarter > 2 & dat$quarter <= 4)
  } 
  rows.to.flip <- c(rows.to.flip.h, rows.to.flip.a)
  dat[rows.to.flip, cols.to.flip.x] <- 94 - dat[rows.to.flip, cols.to.flip.x]
  dat[rows.to.flip, cols.to.flip.y] <- 50 - dat[rows.to.flip, cols.to.flip.y]
  dat[, cols.to.flip.y] <- 50 - dat[, cols.to.flip.y] #flipping data to correct origin
  gc()
  return(dat)
}

# subset to only plays where ballcarrier in offensive halfcourt
# and clock is moving
offensive.ballcarrier <- function(dat) {
  in.play.indicator <- dat$x < 47 & dat$x > 0 & dat$y < 50 & dat$y > 0
  clock.diff <- c(0, diff(dat$game_clock))
  clock.moving.indicator <- clock.diff < -.01 & clock.diff > -.1
  goods <- which(in.play.indicator & clock.moving.indicator)
  return(dat[goods, ])
}

# ID for each player-touch sequence
get.touchID <- function(dat) {
  ent.diff <- diff(dat$entity)
  clock.diff <- diff(dat$game_clock)
  new.entity <- c(1, as.numeric(ent.diff != 0 | abs(clock.diff) > .1))
  touchID <- cumsum(new.entity)
}
