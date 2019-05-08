# Parallelize Stuff
#=========================#
require(badmf)
require(MASS)
library(parallel)
require(dplyr)

no_cores = detectCores() - 1

ns=round(2^seq(5, 9, length.out=8)*1.1)  # want to keep track of training set size
niter <- 100  # number of iterations per simulation
K <- 10  # number of folds for cross validation

# the simulations to call themselves
classifiers <- list(badmf.class.fit, rf.class.fit, dec.tree.class.fit)
names(classifiers) <- c("BaD-MF", "RF", "DT")
class.opts=list(list(ntrees=16L, depth.max=5L, size=1L, mc.cores=1L),
                list(ntrees=16L, depth.max=5L, size=1L, mc.cores=1L,
                     train.params=list(ntrees=32L)),
                list(depth.max=5L, size=1L))

d <- 100

# additional arguments for each simulation scenario
sims.algs <- list(badmf.sims.linear, badmf.sims.rtrunk, badmf.sims.rtrunk,
                  badmf.sims.fat_tails, badmf.sims.cross, badmf.sims.xor2,
                  #badmf.sims.radial,
                  badmf.sims.toep, badmf.sims.qdtoep)
opt_args <- list(list(d=d, K=2, signal.lshift=0),
                 list(d=d),
                 list(d=20, K=4),
                 list(d=d, rotate=TRUE),
                 list(d=d),
                 list(d=20),
                 #list(d=2),
                 list(d=d),
                 list(d=d))
names(sims.algs) <- c("No Signal", "Trunk-2", "Trunk-4", "Spread", "Cross", "XOR", #"Sphere",
                      "Toep", "QDToep")

simulations <- do.call(c, lapply(1:length(sims.algs), function(i) {
  do.call(c, lapply(ns, function(n) {
    lapply(1:niter, function(j) {
      list(sim.func=sims.algs[[i]], args=c(list(n=n), opt_args[[i]]), n=n, d=opt_args[[i]]$d, sim=names(sims.algs)[i], iter=j)
    })
  }))
}))

# Setup Algorithms
#=========================#
opath <- '../data/sims'
tmp.opath <- '../data/sims/smol'
dir.create('../data')
dir.create(opath)
dir.create(tmp.opath)
results <- mclapply(simulations, function(sim) {
  sim_dat <- do.call(sim$sim.func, sim$args)
  X <- sim_dat$X; Y <- sim_dat$Y

  #Randomly shuffle the data
  shuf <- sample(length(Y))
  X <- X[shuf,]; Y <- Y[shuf]

  #Create K equally size folds
  folds <- cut(seq(1, length(Y)), breaks=K, labels=FALSE)

  # loop over the folds
  result <- do.call(rbind, lapply(1:K, function(fold) {
    # set up training and testing set
    test.idx <- which(folds == fold, arr.ind=TRUE)
    X.train <- X[-test.idx,,drop=FALSE]; Y.train <- Y[-test.idx]
    X.test <- X[test.idx,,drop=FALSE]; Y.test <- Y[test.idx]
    # loop over the classifiers of interest
    do.call(rbind, lapply(1:length(classifiers), function(i) {
      tryCatch({
        # fit the forest/tree to the training set
        fit <- do.call(classifiers[[i]], c(list(X=X.train, Y=Y.train), class.opts[[i]]))
        # run prediction using S3 class
        Y.pred <- predict(fit, X.test)
        data.frame(Simulation=sim$sim, iter=sim$iter, fold=fold, Classifier=names(classifiers)[i],
                   n=sim$n, misclass.rate=mean(as.character(Y.pred) != Y.test))
      }, error=function(e) return(NULL))
    }))
  }))
  saveRDS(result, file.path(tmp.opath, paste(sim$sim, "_n", sim$n, "_i", sim$iter, '.rds', sep="")))
  return(result)
}, mc.cores=no_cores)

# Aggregate and save
#=================================#
resultso <- results %>%
  bind_rows()
saveRDS(resultso, file.path(opath, paste('badmf_results', '.rds', sep="")))
