# Some random subsetting of results

ents <- c("226806", "214152", "296572", "172631", "61849", "2989", "172537", "3501",
"3235", "330047", "292404", "253980", "66858", "3497", "410831", "75344", "3329", 
"3315", "3062", "3257", "265013", "173147")

# 226806,214152,296572,172631,61849,2989,172537,3501,3235,330047,292404,253980,66858,3497,410831,75344,3329,3315,3062,3257,265013,173147

shrinkMicro <- function(obj) {
  obj$io.x <- obj$io.x$inla_out[c("summary.fixed", "summary.random", "mode")]
  obj$io.x$mode <- obj$io.x$mode["theta"]
  obj$io.y <- obj$io.y$inla_out[c("summary.fixed", "summary.random", "mode")]
  obj$io.y$mode <- obj$io.y$mode["theta"]
  return(obj)
}

for(i in 1:length(ents)) {
  print(sprintf("%i of %i", i, length(ents)))
  load(sprintf("/n/gstore/Labs/Bornn_Lab/hazards13_old/%s/micro.Rdata", ents[i]))
  with.ball <- shrinkMicro(with.ball)
  without.ball <- shrinkMicro(without.ball)
  save(with.ball, without.ball, file=sprintf("/n/gstore/Labs/Bornn_Lab/epv_data/2013/models/micros/%s.Rdata", ents[i]))
}

for(i in 1:length(ents)) {
  print(sprintf("%i of %i", i, length(ents)))
  load(sprintf("%s/tmats/%s.Rdata", data.dir, ents[i]))
  #load(sprintf("/n/gstore/Labs/Bornn_Lab/epv_data/2013/%s/rmat.Rdata", ents[i]))
  tmat.ind <- tmat_ind[c("micros", "passes1c", "passes2c", "passes3c", "passes4c", "absorbsc")]
  names(tmat.ind) <- c("micros", "passes1", "passes2", "passes3", "passes4", "absorbs")
  tmat.pos <- tmat_pos[c("micros", "passes1c", "passes2c", "passes3c", "passes4c", "absorbsc")]
  names(tmat.pos) <- c("micros", "passes1", "passes2", "passes3", "passes4", "absorbs")
  save(tmat.ind, tmat.pos, file=sprintf("%s/tmats/%s.Rdata", data.dir, ents[i]))
}

inlas <- c("TAKE", "MAKE", "PASS1", "PASS2", "PASS3", "PASS4", "TO")
for(i in 1:length(inlas)) {
  load(sprintf("%s/INLA_%s.Rdata", data.dir, inlas[i]))
  inla.out <- inla.out.lite
  param.names <- row.names(inla.out$summary.fixed)
  n <- nrow(players)
  player.params <- matrix(NA, nrow=n, ncol=length(param.names))
  y.fix <- inla.out$summary.fixed[, "mean"]
  temp <- names(inla.out$summary.random)
  basis.inds <- c(which(temp == "p.int"), grep("p.b[0-9][0-9]*", temp))
  cov.inds <- setdiff(seq(length(inla.out$summary.random)), basis.inds)
  for(pl in 1:n) {
    y.rand <- c(inla.out$summary.random$p.int[pl, "mean"], 
                sapply(cov.inds, 
                       function(k) inla.out$summary.random[[k]][pl, "mean"]),
                inla.out$summary.random$p.b1[pl + n * (1:n.basis), "mean"])
    player.params[pl, ] <- y.fix + y.rand
  }
  names(player.params) <- param.names
  save(player.params, file=sprintf("%s/%s_params.Rdata", data.dir, inlas[i]))
}

names(player.params) <- param.names
load(sprintf("%s/INLA_TAKE.Rdata", data.dir))
inla.out <- inla.out.lite
# coefficients for time-varying covariates in shot-taking hazard model
inla.out$summary.fixed[, 1:2]

