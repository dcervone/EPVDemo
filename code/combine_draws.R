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

all.epv <- data.frame(time=datex$time, poss=datex$entity, poss_x=datex$x, off1_ent=datex$off1_ent,
                      off2_ent=datex$off2_ent, off3_ent=datex$off3_ent, off4_ent=datex$off4_ent, epv=epv, probs=probs, vals=vals,
                      probs.now=probs.now, vals.now=vals.now)

all.epv <- data.table(all.epv)
all.dat <- data.table(dat)
setkey(all.epv, time)
setkey(all.dat, time)
dat.epv <- merge(all.dat, all.epv, all.x=TRUE)
epv.temp <- dat.epv$epv
dt <- c(0, diff(dat.epv$time))
for(i in 2:length(epv.temp)) {
  if(is.na(epv.temp[i]) & dt[i] > 35 & dt[i] < 55)
    epv.temp[i] <- epv.temp[i-1]
}
dat.epv$epv <- epv.temp
#go through iterate, repeat epv for 3,4,7,22,25 until non-NA event_id
dat.epv[which(dat.epv$poss_x > 47), ] <- NA

save(dat.epv, file=sprintf("%s/full.epv.results.Rdata", data.dir))

