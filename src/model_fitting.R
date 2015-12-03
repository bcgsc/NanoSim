library(stats4)

##### Read in data
prefix <- commandArgs(trailingOnly = TRUE)
mis_file <- paste(prefix, "_mis.hist", sep="")
ins_file <- paste(prefix, "_ins.hist", sep="")
del_file <- paste(prefix, "_del.hist", sep="")

mis <- read.delim(mis_file)
mis.freq <- mis$Mismatches.
mis.os <- rep(mis$number.of.bases - 1, mis$Mismatches.)
mis.cum <- ecdf(mis.os)(0:length(mis.freq))

ins <- read.delim(ins_file)
ins.freq <- ins$Insertions.
ins.os <- rep(ins$number.of.bases, ins$Insertions.)
ins.cum <- ecdf(ins.os)(1:length(ins.freq))

del <- read.delim(del_file)
del.freq <- del$Deletions.
del.os <- rep(del$number.of.bases, del$Deletions.)
del.cum <- ecdf(del.os)(1: length(del.freq))

##### Mix model functions
# Poisson-Geometric distribution
rpoisgeommix <- function(n,lambda, p, prob) {
  ifelse(runif(n)<prob,
         rpois(n, lambda),
         rgeom(n,p))
}

dpoisgeommix <- function(x,lambda,p, prob,log=FALSE) {
  r <- prob*dpois(x,lambda)+(1-prob)*dgeom(x,p)
  if (log) log(r) else r
}

ppoisgeommix <- function(x, lambda, p, prob, log=FALSE){
  r <- prob*ppois(x, lambda) + (1-prob)*pgeom(x, p)
  if (log) log(r) else r
}


# Wei-Geom distribution
rweigeommix <- function(n,shape, p, prob, scale = 1) {
  ifelse(runif(n)<prob,
         rweibull(n, shape, scale),
         rgeom(n,p))
}

dweigeommix <- function(x,shape,p, prob,scale = 1, log=FALSE) {
  r <- prob*dweibull(x, shape, scale)+(1-prob)*dgeom(x,p)
  if (log) log(r) else r
}

pweigeommix <- function(x,shape,p, prob,scale = 1, log=FALSE) {
  r <- prob*pweibull(x, shape, scale)+(1-prob)*pgeom(x,p)
  if (log) log(r) else r
}


##### MLE fitting
LL.mis <- function(lambda, p, prob) {
  R = dpoisgeommix(mis.os, lambda, p, prob)
  -sum(log(R))
}
mis.fit <- mle(LL.mis, start = list(lambda=0.5, p = 0.5, prob = 0.5))

LL.del <- function(shape, p, prob, scale) {
  -sum(pweigeommix(x, shape, p, prob, scale, log=TRUE))
}
del.fit <- mle(LL.del, start = list(shape=0.5, p = 0.5, prob = 0.5, scale = 0.3))

write.table(x.fit.table, "test.tsv", quote = FALSE, sep = "\t", row.names = TRUE)
