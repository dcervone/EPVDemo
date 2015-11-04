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

