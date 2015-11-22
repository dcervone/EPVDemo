# **********************
#
# Draws a single EPV curve estimate for a full game of data
# best implemented on a computing cluster
#
# **********************

code.dir <- "../code"
data.dir <- "../data"

# code.dir <- "~/xyhoops/XYHoops/dlc_src/new/demo/code"
# data.dir <- "~/xyhoops/XYHoops/dlc_src/new/demo/data"

args <- as.vector(commandArgs(trailingOnly=TRUE)) # EPV draw number
draw.num <- as.numeric(args[1])

load(sprintf("%s/playerbases.Rdata", data.dir))
players <- read.csv(sprintf("%s/players2013.csv", data.dir))

source(sprintf("%s/constants.R", code.dir))
source(sprintf("%s/data_formatting.R", code.dir))
source(sprintf("%s/covariates.R", code.dir))
source(sprintf("%s/graphics.R", code.dir))
source(sprintf("%s/parameters.R", code.dir))
source(sprintf("%s/EPV_calcs.R", code.dir))


load(sprintf("%s/tdat.Rdata", data.dir))

def.micro <- microDefModel(tdat)

hyper <- getHyperParams(tdat)
ev.out <- evLineups(tdat)
  
draw.raw <- multiresDraw(tdat, hyper, def.micro, ev.out, nmic=50, save.positions=F)
draw <- compressEPV(tdat, draw.raw$fv.epv.list)

dir.create(sprintf("%s/EPVdraws/", data.dir))
save(draw, file=sprintf("%s/EPVdraws/%03.f.Rdata", data.dir, draw.num))
