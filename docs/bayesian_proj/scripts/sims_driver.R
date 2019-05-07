# Parallelize Stuff
#=========================#
require(lolR)
require(badmf)
require(MASS)
library(parallel)
require(mgc)
require(dplyr)
no_cores = detectCores() - 1

ns=round(2^seq(5, 10, length.out=15)*1.1)  # want to keep track of training set size
niter <- 100  # number of iterations per simulation
K <- 10  # number of folds for cross validation

# the simulations to call themselves
classifiers <- list(badmf.class.fit, rf.class.fit, tree.class.fit)
names(classifiers) <- c("BaD-MF", "RF", "DT")
class.opts=list(list(ntrees=9L, depth.max=4L, size=2L, mc.cores=1L),
                list(ntrees=9L, depth.max=4L, size=2L, mc.cores=1L),
                list(depth.max=10L, size=2L))

# additional arguments for each simulation scenario
sims.algs <- list(badmf.sims.linear, badmf.sims.rtrunk, badmf.sims.rtrunk,
                  badmf.sims.fat_tails, badmf.sims.cross, badmf.sims.xor,
                  badmf.sims.radial)
opt_args <- list(list(K=2, signal.lshift=0),
                 list(),
                 list(K=4),
                 list(rotate=TRUE),
                 list(robust=0.1),
                 list(K=4))
names(sims.algs) <- c("No Signal", "Trunk-2", "Trunk-3", "Fat-Tails (D=1000)", "Cross", "Sphere")

simulations <- do.call(c, lapply(1:length(sims.algs), function(i) {
  do.call(c, lapply(ns, function(n) {
    lapply(1:niter, function(j) {
      list(sim.func=sims.algs[[i]], args=c(list(n=n, d=100), opt_args[[i]]), n=n, d=100, sim=names(sims.algs)[i], iter=j)
    })
  }))
}))

# Setup Algorithms
#=========================#
opath <- '../data/sims'
results <- mclapply(simulations, function(sim) {
  sim_dat <- do.call(sim$sim.func, sim$args)
  X <- sim_dat$X; Y <- sim_dat$Y

  #Randomly shuffle the data
  shuf <- sample(length(Y))
  X <- X[shuf,]; Y <- Y[shuf]

  #Create K equally size folds
  folds <- cut(seq(1, length(Y)), breaks=K, labels=FALSE)

  # loop over the folds
  lapply(1:K, function(fold) {
    # set up training and testing set
    test.idx <- which(folds == fold, arr.ind=TRUE)
    X.train <- X[-test.idx,] Y.train <- Y[-test.idx,]
    X.test <- X[test.idx,]; Y.test <- Y[test.idx,]
    # loop over the classifiers of interest
    lapply(1:length(classifiers), function(i) {
      tryCatch({
        # fit the forest/tree to the training set
        fit <- do.call(classifiers[[i]], c(X.train, Y.train, class.opts[[i]]))
        # run prediction using S3 class
        Y.pred <- predict(fit, X.test)
        data.frame(sim=sim$sim, iter=sim$iter, fold=fold, alg=names(algs)[i],
                   n=sim$n, misclass.rate=mean(Y.pred != Y.test)))
      }, error=function(e) return(NULL))
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()

}, mc.cores=no_cores)

# Aggregate and save
#=================================#
resultso <- results %>%
  bind_rows()
saveRDS(resultso, file.path(opath, paste('badmf_results', '.rds', sep="")))
