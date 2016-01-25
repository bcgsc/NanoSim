library(stats4)

##### Read in data
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
# prefix <- "UCSC_phase1b/ecoli"
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
rpoisgeommix <- function(n, l, p, w) {
  ifelse(runif(n) < w,
         rpois(n, l),
         rgeom(n, p))
}

dpoisgeommix <- function(x, l, p, w,log=FALSE) {
  r <- w * dpois(x, l) + (1 - w) * dgeom(x, p)
  if (log) log(r) else r
}

ppoisgeommix <- function(x, l, p, w, log=FALSE){
  r <- w * ppois(x, l) + (1 - w) * pgeom(x, p)
  if (log) log(r) else r
}


# Wei-Geom distribution
rweigeommix <- function(n, shape, p, w, scale = 1) {
  ifelse(runif(n) < w,
         rweibull(n, shape, scale),
         rgeom(n, p))
}

dweigeommix <- function(x, shape, p, w, scale = 1, log=FALSE) {
  r <- w * dweibull(x, shape, scale) + (1 - w) * dgeom(x, p)
  if (log) log(r) else r
}

pweigeommix <- function(x, shape,p, w, scale = 1, log=FALSE) {
  r <- w * pweibull(x, shape, scale) + (1 - w) * pgeom(x, p)
  if (log) log(r) else r
}

weigeom.fit <- function(p, shape1, shape2, scale1, scale2, p1, p2, prob1, prob2, step){
  d <- 1
  for (i in seq(shape1, shape2, step)){
    for (j in seq(scale1, scale2, step)){
      for (n in seq(p1, p2, step)){
        for (m in seq(prob1, prob2, step)){
          e <- pweigeommix(1:length(p), i, n, m, j)
          d1 <- max(abs(e - p))
          if (d1 < d){
            d = d1
            I <- i
            J <- j
            N <- n
            P <- m
          }
        }
      }
    }
  }
  estimate <- c(I, J, N, P, d)
  estimate
}

##### Model fitting
# Mismatch
LL.mis <- function(lambda, prob, weight) {
  R = dpoisgeommix(mis.os, lambda, prob, weight)
  -sum(log(R))
}

mis_fit_func <- function(ll) {
  s <- 0.1
  while (exists("error_model")==FALSE) {
    tryCatch(error_model <- mle(LL.mis, start = list(lambda=s, prob = s, weight = s)),
             error = function(e){print("Try different initial value of MLE")})
    s <- s + 0.1
  }
  return(error_model)
}

mis.fit <- mis_fit_func(LL.mis)

# Indel
indel_search <- function(model.cum){
  tmp <- weigeom.fit(model.cum, 0.1, 1.5, 0.1, 1.5, 0.1, 1, 0.1, 1, 0.02)
  tmp <- weigeom.fit(model.cum, max(tmp[1]-0.02, 0.001), tmp[1]+0.02, max(tmp[2]-0.02, 0.001), tmp[2]+0.02,
                     max(tmp[3]-0.02, 0.001), min(tmp[3]+0.02, 1), max(tmp[4]-0.02, 0.001), min(tmp[4]+0.02, 1), 0.001)
  start <- c(max(tmp[1]-0.001, 0.0001), tmp[1]+0.001, max(tmp[2]-0.001, 0.0001), tmp[2]+0.001,
             max(tmp[3]-0.001, 0.0001), min(tmp[3]+0.001, 1), max(tmp[4]-0.001, 0.0001), min(tmp[4]+0.001, 1))
  end <- weigeom.fit(model.cum, start[1], start[2], start[3], start[4], start[5], start[6], start[7], start[8], 0.0001)
  i <- 0
  p_value <- end[5]
  while (Reduce("|", end %in% start) | i < 1000){
    last_end <- end
    start <- c(max(end[1]-0.001, 0.0001), end[1]+0.001, max(end[2]-0.001, 0.0001), end[2]+0.001,
               max(end[3]-0.001, 0.0001), min(end[3]+0.001, 1), max(end[4]-0.001, 0.0001), min(end[4]+0.001, 1))
    end <- weigeom.fit(model.cum, start[1], start[2], start[3], start[4], start[5], start[6], start[7], start[8], 0.0001)
    if(p_value <= end[5]) {
      end <- last_end
      break
    }
    i <- i + 1
    p_value <- end[5]
  }
  return(end)
}

ins.fit <- indel_search(ins.cum)
del.fit <- indel_search(del.cum)

##### Write to file
model_fit.table <- data.frame(Type = c("mismatch", "insertion", "deletion"),
                              lambda = c(coef(mis.fit)["lambda"], ins.fit[2], del.fit[2]),
                              k = c(0, ins.fit[1], del.fit[1]),
                              prob = c(coef(mis.fit)["prob"], ins.fit[3], del.fit[3]),
                              weight = c(coef(mis.fit)["weight"], ins.fit[4], del.fit[4]))

out_file <- paste(prefix, "_model_profile", sep="")
write.table(model_fit.table, out_file, row.names = FALSE, quote = FALSE, sep = "\t")
