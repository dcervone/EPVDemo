library(RColorBrewer)
library(fields)

data.plotter <- function(dat, ind, poss=F, ...) {
  par(mar=c(1, 1, 1, 1))
  plot(c(0, 94), c(0, 50), type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  draw.fullcourt(...)
  pos.x <- dat[ind, "x"]
  pos.y <- dat[ind, "y"]
  for(i in 1:5) {
    text(dat[ind, c(sprintf("h%i_x", i), sprintf("h%i_y", i))], labels=i)
    points(dat[ind, c(sprintf("h%i_x", i), sprintf("h%i_y", i))], cex=2.5, lwd=2)
    text(dat[ind, c(sprintf("a%i_x", i), sprintf("a%i_y", i))], labels=i, col="purple")
    points(dat[ind, c(sprintf("a%i_x", i), sprintf("a%i_y", i))], col="purple", cex=2.5, lwd=2)
    if(!(poss))
      next
    if(dat[ind, sprintf("h%i_ent", i)] == dat[ind, "poss"]) {
      pos.x <- dat[ind, sprintf("h%i_x", i)]
      pos.y <- dat[ind, sprintf("h%i_y", i)]
    }
    if(dat[ind, sprintf("a%i_ent", i)] == dat[ind, "poss"]) {
      pos.x <- dat[ind, sprintf("a%i_x", i)]
      pos.y <- dat[ind, sprintf("a%i_y", i)]
    }
  }
  points(dat[ind, c("x", "y")], col="orange", pch=20, cex=2)
  if(poss)
    points(x=pos.x, y=pos.y, pch=5, col="red", cex=3, lwd=2)
}

transformed.data.plotter <- function(dat, ind, ...) {
  par(mar=c(1, 1, 1, 1))
  plot(c(0, 47), c(0, 50), type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  draw.halfcourt(...)
  pos.x <- dat[ind, "x"]
  pos.y <- dat[ind, "y"]
  for(i in 1:5) {
    if(i < 5) {
      text(dat[ind, c(sprintf("off%i_x", i), sprintf("off%i_y", i))], labels=i)
      points(dat[ind, c(sprintf("off%i_x", i), sprintf("off%i_y", i))], cex=2.5, lwd=2)
    }
    text(dat[ind, c(sprintf("def%i_x", i), sprintf("def%i_y", i))], labels=i, col="purple")
    points(dat[ind, c(sprintf("def%i_x", i), sprintf("def%i_y", i))], col="purple", cex=2.5, lwd=2)
  }
  points(dat[ind, c("ball_x", "ball_y")], col="orange", pch=20, cex=2)
  text(dat[ind, c("x", "y")], labels=0)
  points(dat[ind, c("x", "y")], pch=5, col="red", cex=3, lwd=2)
}

## Court plotting utilities
circle <- function(x, y, r, from=0, to=2*pi, lines=FALSE, ...) {
  theta <- seq(from, to, length=100)
  if (lines)
    lines(x + r * cos(theta), y + r * sin(theta), ...)
  else polygon(x + r * cos(theta), y + r * sin(theta), ...)
}

draw.halfcourt <- function(...) {
  rect(0, 0, 94/2, 50, ...)
  points(5.25, 25, cex = 2)
  segments(47, 0, 47, 50, ...)
  theta1 <- acos((25 - 35 / 12) / 23.75)
  circle(5.25, 25, 23.75, -pi / 2 + theta1, pi / 2 - theta1, TRUE, ...)
  segments(0, 35 / 12, 5.25 + 23.75 * sin(theta1), 35 / 12, ...)
  segments(0, 50 - 35/12, 5.25 + 23.75 * sin(theta1), 50 - 35 / 12, ...)
  circle(19, 25, 6, -pi / 2, pi / 2, TRUE, ...)
  circle(19, 25, 6, pi / 2, 3 * pi / 2, TRUE, lty = 2, ...)
  circle(5.25, 25, 4, -pi / 2, pi / 2, TRUE, ...)
  rect(0, 17, 19, 33, border = "gray", ...)
}

draw.fullcourt <- function(...) {
  rect(0, 0, 94, 50, ...)
  points(c(5.25, 94 - 5.25), c(25, 25), cex = 2)
  segments(47, 0, 47, 50, ...)
  circle(47, 25, 8, ...)
  circle(47, 25, 2, col = "lightgray", ...)
  theta1 <- acos((25 - 35 / 12) / 23.75)
  circle(5.25, 25, 23.75, -pi / 2 + theta1, pi / 2 - theta1, TRUE, ...)
  circle(94 - 5.25, 25, 23.75, pi / 2 + theta1, 3 * pi / 2 - theta1, TRUE, ...)
  segments(0, 35/12, 5.25 + 23.75 * sin(theta1), 35 / 12, ...)
  segments(0, 50 - 35 / 12, 5.25 + 23.75 * sin(theta1), 50 - 35 / 12, ...)
  segments(94, 35 / 12, 94 - 5.25 - 23.75 * sin(theta1), 35 / 12, ...)
  segments(94, 50 - 35 / 12, 94 - 5.25 - 23.75 * sin(theta1), 50 - 35 / 12, ...)
  circle(19, 25, 6, -pi/2, pi/2, TRUE, ...)
  circle(19, 25, 6, pi/2, 3 * pi/2, TRUE, lty = 2, ...)
  circle(94 - 19, 25, 6, pi/2, 3 * pi/2, TRUE, ...)
  circle(94 - 19, 25, 6, -pi/2, pi/2, TRUE, lty = 2, ...)
  circle(5.25, 25, 4, -pi/2, pi/2, TRUE, ...)
  circle(94 - 5.25, 25, 4, pi/2, 3 * pi/2, TRUE, ...)
  rect(0, 17, 19, 33, border = "gray", ...)
  rect(94, 17, 94 - 19, 33, border = "gray", ...)
}

colo <- colorRampPalette(c("white", "#91CF60", "yellow", "red"))
spatialPlot0 <- function(z, cexy=0.5, legend=TRUE, ...) {
  par(mar=c(0.0,0.0,0.0,0.0))
  if(legend) {
    image.plot(seq(0, 47, l=23), seq(0, 50, l=25), matrix(z, nrow=23, ncol=25),
               xlim=c(0, 47), ylim=c(0, 50), zlim=range(z), xlab="", ylab="", xaxt="n", yaxt="n", col=colo(100), useRaster=T, ...)
  } else {
    image(seq(0, 47, l=23), seq(0, 50, l=25), matrix(z, nrow=23, ncol=25),
          xlim=c(0, 47), ylim=c(0, 50), zlim=range(z), xlab="", ylab="", xaxt="n", yaxt="n", col=colo(100), useRaster=T, ...)
  } 
  draw.halfcourt(lwd=2)
}

spatialPlot1 <- function(z, pts1=NULL, pts0=NULL, cexy=0.5, legend=TRUE, ...) {
  par(mar=c(0.0,0.0,0.0,0.0))
  goods <- which(mesh$loc[, 1] >= 0 & mesh$loc[, 1] <= 47 & mesh$loc[, 2] >= 0 & mesh$loc[, 2] <= 50)
  if(legend) {
    image.plot(mesh.proj$x, mesh.proj$y, inla.mesh.project(mesh.proj, z),
               xlim=c(0, 47), ylim=c(0, 50), zlim=range(z[goods]), xlab="", ylab="", xaxt="n", yaxt="n", col=colo(100), useRaster=T, ...)
  } else {
    image(mesh.proj$x, mesh.proj$y, inla.mesh.project(mesh.proj, z),
          xlim=c(0, 47), ylim=c(0, 50), zlim=range(z[goods]), xlab="", ylab="", xaxt="n", yaxt="n", col=colo(100), useRaster=T, ...)
  } 
  draw.halfcourt(lwd=2)
  if(!(is.null(pts1)))
    points(pts1, pch=20, cex=cexy)
  if(!(is.null(pts0)))
    points(pts0, pch=20, col="gray", cex=cexy)
}

spatialPlot2 <- function(z1, z2, pts1=NULL, pts2=pts1, cexy=0.5, legend=TRUE, ...) {
  par(mar=c(0.1,0.5,0.1,0.5))
  par(mfrow=c(1,2))
  goods <- which(mesh$loc[, 1] >= 0 & mesh$loc[, 1] <= 47 & mesh$loc[, 2] >= 0 & mesh$loc[, 2] <= 50)
  if(legend) {
    image.plot(mesh.proj$x, mesh.proj$y, inla.mesh.project(mesh.proj, z1),
               xlim=c(0, 47), ylim=c(0, 50), zlim=range(z1[goods]), xlab="", ylab="", xaxt="n", yaxt="n", col=colo(100), useRaster=T, ...)
    draw.halfcourt(lwd=2)
    if(!(is.null(pts1)))
      points(pts1, pch=20, cex=cexy)
    image.plot(mesh.proj$x, mesh.proj$y, inla.mesh.project(mesh.proj, z2),
               xlim=c(0, 47), ylim=c(0, 50), zlim=range(z2[goods]), xlab="", ylab="", xaxt="n", yaxt="n", col=colo(100), useRaster=T, ...)
    draw.halfcourt(lwd=2)
    if(!(is.null(pts2)))
      points(pts2, pch=20, cex=cexy)
  } else {
    image(mesh.proj$x, mesh.proj$y, inla.mesh.project(mesh.proj, z1),
          xlim=c(0, 47), ylim=c(0, 50), zlim=range(z1[goods]), xlab="", ylab="", xaxt="n", yaxt="n", col=colo(100), useRaster=T, ...)
    draw.halfcourt(lwd=2)
    if(!(is.null(pts1)))
      points(pts1, pch=20, cex=cexy)
    image(mesh.proj$x, mesh.proj$y, inla.mesh.project(mesh.proj, z2),
          xlim=c(0, 47), ylim=c(0, 50), zlim=range(z2[goods]), xlab="", ylab="", xaxt="n", yaxt="n", col=colo(100), useRaster=T, ...)
    draw.halfcourt(lwd=2)
    if(!(is.null(pts2)))
      points(pts2, pch=20, cex=cexy)
  } 
}


vectorPlot <- function(io, ...){
  mesh.proj <- inla.mesh.projector(mesh, dims=c(40,40))
  v.x <- inla.mesh.project(mesh.proj, io$io.x$summary.random$spatial$mean)
  v.y <- inla.mesh.project(mesh.proj, io$io.y$summary.random$spatial$mean)
  goods <- which(mesh.proj$lattice$loc[,1] >= 0 & mesh.proj$lattice$loc[,1] <= 47 & mesh.proj$lattice$loc[,2] >= 0 & mesh.proj$lattice$loc[,2] <= 50)
  v.x <- as.vector(v.x)
  v.y <- as.vector(v.y)
  v.x[setdiff(1:length(v.x), goods)] <- NA
  v.y[setdiff(1:length(v.y), goods)] <- NA
  
  scal <- 1 / (.04 ^ 2) / 2
  mag <- sqrt(v.x ^ 2 + v.y ^ 2)*scal
  lwds <- 0.5 + 3*mag/max(c(3.5, max(mag, na.rm=T)))
  cols <- colorRamp(c("yellow","red"))(mag/max(c(3.5, max(mag, na.rm=T))))
  plot(c(0,47),c(0,50),xlab="", ylab="", xaxt="n", yaxt="n", type="n") #, ...)
  draw.halfcourt(lwd=2)
  for(i in goods){
    arrows(mesh.proj$lattice$loc[i,1], mesh.proj$lattice$loc[i,2], 
           mesh.proj$lattice$loc[i,1] + scal*as.vector(v.x)[i],
           mesh.proj$lattice$loc[i,2] + scal*as.vector(v.y)[i],
           length=0.065,
           angle=30, lwd=lwds[i], col=rgb(cols[i,1], cols[i,2], cols[i,3], maxColorValue=255))
  }
  
}
