##' @title Call CNA using a HMM approach
##' @param ge the input gene expression
##' @param trans.prob the transition probability
##' @param log.mu should the mu be log-scaled?
##' @param max.samps the maximum number of samples accepted.
##' @param perChr run on each chromosome separately?
##' @return a data.frame with info about CNA segments
##' @author Jean Monlong
##' @keywords internal
cnaHMM <- function(ge, trans.prob=.0001, log.mu=TRUE, max.samps=200, perChr=TRUE){
  cells = setdiff(colnames(ge), c("chr","start","end"))
  if(length(cells) > max.samps){
    warning('Too many samples. Downsampling to 200.')
    cells = sample(cells, max.samps)
  }
  ## HMM definition
  M = length(cells) # nb of samples
  N = c("loss","neutral","gain") # states
  Mu <- cbind(rep(1/2,M), rep(1,M), rep(3/2,M)) # mean in each state/samples
  if(log.mu){
    Mu = log(Mu)
  }
  A <- matrix(c(1-trans.prob, trans.prob/2, trans.prob/2,
                trans.prob/2, 1-trans.prob, trans.prob/2,
                trans.prob/2, trans.prob/2, 1-trans.prob),
              ncol = length(N), byrow = TRUE)
  Pi <- c(.2, .6, .2)
  runHMM <- function(ge){
    coord.df = ge[,c("chr","start","end")]
    ge.mat = as.matrix(t(ge[,cells]))
    var.est = stats::var(as.numeric(ge.mat), na.rm=TRUE) # var of samples in state 1
    Sigma <- array(0, dim =c(M,M,length(N)))
    Sigma[,,2] <- Sigma[,,3] <- Sigma[,,1] <- var.est * diag(M) # Independent samples
    HMM.o <- suppressWarnings(RcppHMM::verifyModel(list("Model" = "GHMM", "StateNames" = N,
                                                        "A" = A, "Mu" = Mu,
                                                        "Sigma" = Sigma, "Pi" = Pi)))
    ## Find most likely states
    coord.df$CN = suppressWarnings(RcppHMM::viterbi(HMM.o, ge.mat))
    ## Other metrics
    coord.df$mean = apply(ge.mat,2,mean)
    coord.df
  }
  ge = ge[order(ge$chr, ge$start),]
  if(perChr){
    res = lapply(unique(ge$chr), function(chri){
      runHMM(ge[which(ge$chr == chri),])
    })
    res = do.call(rbind, res)
  } else {
    res = runHMM(ge)
  }
  res
}
