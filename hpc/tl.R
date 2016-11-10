## Tell R where to look for my packages 

print(.libPaths())
lp = "/scratch/gatto/R_local"

## Load libraries in my home directory
library("MLInterfaces", lib.loc = lp)
library("pRoloc", lib.loc = lp)
library("pRolocdata", lib.loc = lp)
library("BiocParallel", lib.loc = lp)
library("Rmpi")

## Load LOPIT and GOCC data from pRolocdata package
data("hyperLOPIT2015")
data("hyperLOPIT2015goCC")


## keep only 7 overlapping classes for example run for workflow
## rename the following classes as "unknown"
torm  <- c("40S Ribosome", "60S Ribosome",
           "Nucleus - Chromatin", "Nucleus - Non-chromatin")
for (i in seq(torm)) {
  hyperLOPIT2015 <- fDataToUnknown(hyperLOPIT2015,
                                   from = torm[i])
  hyperLOPIT2015goCC <- fDataToUnknown(hyperLOPIT2015goCC,
                                       from = torm[i])
}


## Set params for parallelisation
par <- SnowParam(255L, type = "MPI",
                 stop.on.error = FALSE,
                 log = TRUE,
                 threshold = "DEBUG")


## run opt
res <- knntlOptimisation(hyperLOPIT2015, hyperLOPIT2015goCC, 
                            fcol = "markers", 
                            length.out = 4, 
                            times = 50,
                            xval = 5, k = c(3,3), 
                            BPPARAM = par)

save(res, file = "tl-hyperlopit.rda")

# Tell all slaves to close down, and exit the program
mpi.quit()