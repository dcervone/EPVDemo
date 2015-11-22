# **********************
#
# Generate spatiotemporal covariates used in multiresolution transition models
#
# **********************

#distance to nearest defender
getNdef <- function(dat, colx, coly) {
  def_dist1 <- sqrt((dat[, colx] - dat$def1_x)^2 + (dat[, coly] - dat$def1_y)^2)
  def_dist2 <- sqrt((dat[, colx] - dat$def2_x)^2 + (dat[, coly] - dat$def2_y)^2)
  def_dist3 <- sqrt((dat[, colx] - dat$def3_x)^2 + (dat[, coly] - dat$def3_y)^2)
  def_dist4 <- sqrt((dat[, colx] - dat$def4_x)^2 + (dat[, coly] - dat$def4_y)^2)
  def_dist5 <- sqrt((dat[, colx] - dat$def5_x)^2 + (dat[, coly] - dat$def5_y)^2)
  def_dist <- data.frame(def_dist1, def_dist2, def_dist3, def_dist4, def_dist5)
  ndef <- apply(def_dist, 1, function(x) min(x, na.rm=T))
  return(log(1 + ndef))
}

#velocity/acceleration
getVeloc <- function(dat) {
  vx <- c(0, 0, 0, 0, diff(dat$x, 4) / diff(dat$time, 4)) * 1000
  vy <- c(0, diff(dat$y) / diff(dat$time)) * 1000
  ax <- c(0, diff(vx) / diff(dat$time)) * 1000
  ay <- c(0, diff(vy) / diff(dat$time)) * 1000
  return(list(vx=vx, vy=vy, ax=ax, ay=ay))
  #indicataor for acceleration in same/opposite direction as velocity
  #acel <- (sign(vx) == sign(ax))*(sign(vy) == sign(ay))
  #decel <- (sign(vx) != sign(ax))*(sign(vy) != sign(ay))
}

#indicator for has a dribble happened on this possession?
getDrib <- function(dat) {
  time.diff <- diff(dat$time)
  ent.diff <- diff(dat$entity)
  ends <- c(which(time.diff < 0 | time.diff > 100 | ent.diff != 0), nrow(dat))
  starts <- c(1, which(time.diff < 0 | time.diff > 100 | ent.diff != 0) + 1)
  posses <- lapply(1:length(starts), function(i) starts[i]:ends[i])
  
  dribble <- rep(0, nrow(dat))
  for(i in 1:length(posses)){
    these <- which(dat$event_id[posses[[i]]] == 21)
    if(length(these) > 0){
      this <- min(these)
      dribble[posses[[i]][this]:tail(posses[[i]], 1)] <- 1
    }
  }
  return(dribble)
}
  
#distance ball has traveled in last 1 second
getBallDist <- function(dat) {
  ball.vx <- c(0,diff(dat$ball_x) / diff(dat$time)) * 1000
  ball.vy <- c(0,diff(dat$ball_y) / diff(dat$time)) * 1000
  ball.vz <- c(0,diff(dat$ball_z) / diff(dat$time)) * 1000
  ball.sp <- sqrt(ball.vx^2 + ball.vy^2 + ball.vz^2)
  bads <- which(ball.sp > 100)
  if(length(bads) > 0)
    ball.sp[bads] <- NA
  ball.lastsec <- c(rep(0, 25), sapply(1:(nrow(dat) - 25), function(i) mean(ball.sp[i + 0:24], na.rm=T)))
  return(ball.lastsec)
}

#ranked distance to each teammate (doff1 = 1 corresponds to teammate 1 being closest)
getDoff <- function(dat) {
  doff1 <- sqrt((dat$x - dat$off1_x)^2 + (dat$y - dat$off1_y)^2)
  doff2 <- sqrt((dat$x - dat$off2_x)^2 + (dat$y - dat$off2_y)^2)
  doff3 <- sqrt((dat$x - dat$off3_x)^2 + (dat$y - dat$off3_y)^2)
  doff4 <- sqrt((dat$x - dat$off4_x)^2 + (dat$y - dat$off4_y)^2)
  
  doff <- data.frame(doff1, doff2, doff3, doff4)
  orank <- t(apply(doff, 1, order))
  doff1 <- orank[, 1]
  doff2 <- orank[, 2]
  doff3 <- orank[, 3]
  doff4 <- orank[, 4]
  return(list(doff1=doff1, doff2=doff2, doff3=doff3, doff4=doff4))
}

#score for how "open" a teammate is for a pass, from 0 (completely open) to 5 (not open at all)
getOpen <- function(dat) {
  foo_defbox <- function(i) {
    this <- dat[i, ]
    
    defs <- t(matrix(c(this$def1_x, this$def1_y,
                       this$def2_x, this$def2_y,
                       this$def3_x, this$def3_y,
                       this$def4_x, this$def4_y,
                       this$def5_x, this$def5_y), nrow=2, ncol=5) - c(this$x, this$y))
    
    offs <- t(matrix(c(this$off1_x, this$off1_y,
                       this$off2_x, this$off2_y,
                       this$off3_x, this$off3_y,
                       this$off4_x, this$off4_y), nrow=2, ncol=4) - c(this$x, this$y))
    
    angles <- apply(offs, 1, function(x) atan(x[1] / x[2])) 
    
    rot_mats <- lapply(angles, function(x) matrix(c(cos(x), -sin(x), sin(x), cos(x)), nrow=2, ncol=2, byrow=T))
    
    ots <- lapply(1:4, function(i) rot_mats[[i]] %*% offs[i,])  
    dts <- lapply(1:4, function(i) rot_mats[[i]] %*% t(defs))
    slopy <- tan(pi / 2 - .15)
    
    dbox <- tryCatch(sapply(1:4, function(j){
      if(ots[[j]][2, 1] < 0) {
        ots[[j]][2, 1] <- -ots[[j]][2, 1]
        dts[[j]][2, ] <- -dts[[j]][2, ]
      }
      sum(
        sapply(1:5, function(k) min(c(1, ((abs(dts[[j]][2, k]) + 4) / slopy) / abs(dts[[j]][1, k])))^4) * (abs(dts[[j]][2, ]) < abs(ots[[j]][2, 1]))
      )
    }), error = function(e) e)
    if(inherits(dbox, "error"))
      dbox <- rep(NA, 4)
    return(dbox)
  }
  
  dboxmat <- matrix(NA, nrow=nrow(dat), ncol=4)
  for(i in 1:nrow(dat)){
    dboxmat[i, ] <- foo_defbox(i)
  }

  ddef1 <- dboxmat[, 1]
  ddef2 <- dboxmat[, 2]
  ddef3 <- dboxmat[, 3]
  ddef4 <- dboxmat[, 4]
  
  return(list(ddef1=ddef1, ddef2=ddef2, ddef3=ddef3, ddef4=ddef4))
}

#indicator for beyond 3 pt line
getThree <- function(dat) {
  bask.dist <- sqrt((dat$x - 4.75)^2 + (dat$y - 25)^2)
  threept <- as.numeric((dat$x < 14 & (dat$y < 3 | dat$y > 47)) | (dat$x > 14 & bask.dist > 23.75))
  return(list(bask.dist=bask.dist, threept=threept))
}

#more court location indicators
entityPlace <- function(x, y){
  if(is.na(x) || is.na(y)) return(NA)
  dist <- sqrt((x - 4.75)^2 + (y - 25)^2)
  if(x < 4 & (y > 3 & y < 47)) {
    place <- "behind"
  } else if( (x < 14)  & ((y > 0 & y < 3) | (y > 47 & y < 50))) {
    place <- "cor3"
  } else if( (x > 14 & x < 32)  & (dist > 23.75 & dist < 27)){
    place <- "cen3"
  } else if(dist <= 4){
    place <- "rest"
  } else if( (x > 4 & x < 19)  & (y > 19 & y < 31)){
    place <- "key"
  } else {
    place <- "other"
  }
  if(dist <= 23.75 & !(place %in% c("key","rest","cor3","behind"))) place <- "per"
  return(place)
}

getPlace <- function(dat) apply(dat[, c("x", "y")], 1, function(r) entityPlace(r[1], r[2])) 

#distance behind 3 pt arc, distance from basket in 2 pt
placeDist <- function(x, y) {
  if(is.na(x) | is.na(y))
    return(list(val_2=NA, val_3=NA, val_behind=NA))
  dist <- sqrt((x - 4.75)^2 + (y - 25)^2)
  if(x < 4 & (y > 3 & y < 47)){
    val_behind <- 4 - x
    val_3 <- 0
    val_2 <- 0
  } else if((x < 14) &  ((y > 0 & y < 3) | (y > 47 & y < 50))){
    val_3 <- max(c(3 - y, y - 47))
    val_behind <- 0
    val_2 <- 0
  } else if(x > 4 & dist > 23.75) {
    val_3 <- dist - 23.75
    val_2 <- 0
    val_behind <- 0    
  } else {
    val_behind <- 0
    val_3 <- 0
    val_2 <- dist
  }
  return(list(val_2=val_2, val_3=val_3, val_behind=val_behind))
}

getFeet <- function(dat) {
  oo <- lapply(1:nrow(dat), function(j) placeDist(dat$x[j], dat$y[j]))
  val_2 <- sapply(1:nrow(dat), function(j) oo[[j]]$val_2)
  val_3 <- sapply(1:nrow(dat), function(j) oo[[j]]$val_3)
  val_behind <- sapply(1:nrow(dat), function(j) oo[[j]]$val_behind)
  return(list(val_2=val_2, val_3=val_3, val_behind=val_behind))
}


#one ring to unite them all
getAllCovars <- function(dat){
  ndef <- getNdef(dat, "x", "y")
  getVeloc.o <- getVeloc(dat)
  dribble <- getDrib(dat)
  ball.lastsec <- getBallDist(dat)
  getDoff.o <- getDoff(dat)
  getOpen.o <- getOpen(dat)
  getThree.o <- getThree(dat)
  eP <- getPlace(dat)
  
  vx <- getVeloc.o$vx
  vy <- getVeloc.o$vy
  ax <- getVeloc.o$ax
  ay <- getVeloc.o$ay
  
  doff1 <- as.numeric(getDoff.o$doff1 == 1)
  doff2 <- as.numeric(getDoff.o$doff2 == 1)
  doff3 <- as.numeric(getDoff.o$doff3 == 1)
  doff4 <- as.numeric(getDoff.o$doff4 == 1)
  
  ddef1 <- getOpen.o$ddef1
  ddef2 <- getOpen.o$ddef2
  ddef3 <- getOpen.o$ddef3
  ddef4 <- getOpen.o$ddef4
  
  bask.dist <- getThree.o$bask.dist
  threept <- getThree.o$threept
  
  fo.o <- getFeet(dat)
  val_2 <- fo.o$val_2
  val_3 <- fo.o$val_3
  val_b <- fo.o$val_behind
  
  extras <- data.frame(ndef, dribble, ball.lastsec, eP, vx, vy, ax, ay, doff1, 
                       doff2, doff3, doff4, ddef1, ddef2, ddef3, ddef4,
                       bask.dist, threept, val_2, val_3, val_b)
  return(extras)
}

