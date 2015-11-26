# **********************
#
# Combines independent EPV draws output from epv_draw.R
#
# **********************

code.dir <- "../code"
data.dir <- "../data"
library(data.table)

load(file=sprintf("%s/tdat.Rdata", data.dir))
dat <- read.csv(file=sprintf("%s/2013_11_01_MIA_BKN.csv", data.dir))

fls <- list.files(sprintf("%s/EPVdraws", data.dir), full.names=FALSE)
draw.list <- vector("list", length(fls))

for(i in 1:length(draw.list)) {
  print(sprintf("loading %i of %i", i, length(draw.list)))
  load(sprintf("%s/EPVdraws/%s", data.dir, fls[i]))
  draw.list[[i]] <- draw
}

epv <- lapply(draw.list, function(obj) obj$epv)
probs <- lapply(draw.list, function(obj) obj$probs)
vals <- lapply(draw.list, function(obj) obj$vals)
epv <- Reduce("+", epv) / length(epv)
probs <- Reduce("+", probs) / length(probs)
vals <- Reduce("+", vals) / length(vals)
probs.now <- draw.list[[1]]$probs.now
vals.now <- draw.list[[1]]$vals.now

# slight spline smoothing of EPV
# ALSO, for first moment of every touch, if no event happens, set equal to second moment of touch
# since micro model will be undefined.
epv.smooth <- epv
chunks <- lapply(1:max(tdat$touchID, na.rm=T), function(i) which(tdat$touchID == i))
for(chunk in chunks) {
  if(length(chunk) <= 2)
    next
  epv[chunk[1]] <- epv[chunk[2]]
  probs[chunk[1], ] <- probs[chunk[2], ]
  vals[chunk[1], ] <- vals[chunk[2], ]
  probs.now[chunk[1], ] <- probs.now[chunk[2], ]
  vals.now[chunk[1], ] <- vals.now[chunk[2], ]

  chunk <- chunk[which(!is.na(epv[chunk]))]
  spl <- smooth.spline(x=chunk, y=epv[chunk], spar=0.2)
  epv.smooth[chunk] <- spl$y
}

epv.table <- data.frame(time=tdat$time, poss=tdat$entity, poss_x=tdat$x,
                      epv=epv, epv.smooth=epv.smooth, probs=probs, vals=vals,
                      probs.now=probs.now, vals.now=vals.now)

save(epv.table, file=sprintf("%s/combined.epv.draws.Rdata", data.dir))