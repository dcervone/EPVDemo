code.dir <- "../code"
data.dir <- "../data"
library(data.table)

fls <- list.files(sprintf("%s/EPVdraws", data.dir), full.names=FALSE)
epv.list <- vector("list", length(fls))

load(file=sprintf("%s/new.dat.Rdata", data.dir))
datex <- new.dat
dat <- read.csv(file=sprintf("%s/2013_11_01_MIA_BKN.csv", data.dir))

for(i in 1:length(epv.list)) {
  print(sprintf("loading %i of %i", i, length(epv.list)))
  load(sprintf("%s/EPVdraws/%s", data.dir, fls[i]))
  epv.list[[i]] <- epv
}

epv <- lapply(epv.list, function(obj) obj$epv)
probs <- lapply(epv.list, function(obj) obj$probs)
vals <- lapply(epv.list, function(obj) obj$vals)
epv <- Reduce("+", epv) / length(epv)
probs <- Reduce("+", probs) / length(probs)
vals <- Reduce("+", vals) / length(vals)
probs.now <- epv.list[[1]]$probs.now
vals.now <- epv.list[[1]]$vals.now

epv.na.diff <- c(0, diff(is.na(epv)))
chunk.starts <- which(epv.na.diff == -1)
chunk.ends <- which(epv.na.diff == 1) - 1
if(!is.na(tail(epv, 1)))
  chunk.ends <- c(chunk.ends, length(epv))
chunks <- lapply(1:length(chunk.starts), function(i) seq(chunk.starts[i], chunk.ends[i]))
epv.smooth <- epv
for(chunk in chunks) {
  if(length(chunk) <= 2)
    next
  spl <- smooth.spline(x=chunk, y=epv[chunk], spar=0.2)
  epv.smooth[chunk] <- spl$y
}

all.epv <- data.frame(time=datex$time, poss=datex$entity, poss_x=datex$x, poss_y=datex$y, off1_ent=datex$off1_ent,
                      off2_ent=datex$off2_ent, off3_ent=datex$off3_ent, off4_ent=datex$off4_ent, 
                      epv=epv, epv.smooth=epv.smooth, probs=probs, vals=vals,
                      probs.now=probs.now, vals.now=vals.now)

all.epv <- data.table(all.epv)
all.dat <- data.table(dat)
setkey(all.epv, time)
setkey(all.dat, time)
dat.epv <- merge(all.dat, all.epv, all.x=TRUE)
epv.temp <- dat.epv$epv.smooth
dt <- c(0, diff(dat.epv$time))
for(i in 2:length(epv.temp)) {
  if(is.na(epv.temp[i]) & dt[i] > 25 & dt[i] < 55)
    epv.temp[i] <- epv.temp[i-1]
}
dat.epv$epv.smooth <- epv.temp
#go through iterate, repeat epv for 3,4,7,22,25 until non-NA event_id
dat.epv[which(dat.epv$poss_x > 47), ] <- NA
# dat.epv[which(dat.epv$poss_x < 0), ] <- NA
# dat.epv[which(dat.epv$poss_y > 50), ] <- NA
# dat.epv[which(dat.epv$poss_y < 0), ] <- NA
epv.na.diff <- c(0, diff(is.na(dat.epv$epv.smooth)))
chunk.starts <- which(epv.na.diff == -1)
chunk.ends <- which(epv.na.diff == 1) - 1
if(!is.na(tail(dat.epv$epv.smooth, 1)))
  chunk.ends <- c(chunk.ends, length(dat.epv$epv.smooth))
chunks <- lapply(1:length(chunk.starts), function(i) seq(chunk.starts[i], chunk.ends[i]))

# rename pbp.seq to possession.id
pbp.seq <- rep(NA, nrow(dat.epv))
for(i in 1:length(chunks)) {
  pbp.seq[chunks[[i]]] <- i
}
dat.epv <- data.frame(dat.epv, pbp.seq=pbp.seq)

save(dat.epv, file=sprintf("%s/full.epv.results.Rdata", data.dir))

if (FALSE) {
  library(animation)
  setwd("~/EPV")
  ani.options(outdir=getwd(), interval=.04, ani.width=1000, ani.height=574)
  makeGIF <- function(inds, name) {
    saveGIF(for(ind in inds) {
      full.plotter(dat.epv, ind, poss=T, legend=T, epv=T, inds=inds)
    }, movie.name=sprintf("%s.gif", name))
    dev.off()
  }
  
  source(sprintf("%s/graphics.R", code.dir))
  makeGIF(which(dat.epv$possession.id == 6), "poss6")
  
  
  
}