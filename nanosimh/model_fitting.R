library(stats4)

# Ensure that the results will be reproducible
set.seed(1)


##### Read in data
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
# prefix <- "R9/1D/test"
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

poisgeom.fit <- function(p, lambda1, lambda2, p1, p2, prob1, prob2, step){
  d <- 1
  for (i in seq(lambda1, lambda2, step)){
    for (n in seq(p1, p2, step)){
      for (m in seq(prob1, prob2, step)){
        e <- ppoisgeommix(0:(length(p)-1), i, n, m)
        d1 <- max(abs(e - p))
        if (d1 < d){
          d = d1
          I <- i
          N <- n
          P <- m
        }
      }
    }
  }
  estimate <- c(I, N, P, d)
  estimate
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
  while (exists("error_model")==FALSE || s > 1) {
    tryCatch(error_model <- mle(LL.mis, start = list(lambda=s, prob = s, weight = s),
                                method = "L-BFGS-B",
                                lower = list(prob = 0.01, weight = 0.01),
                                upper = list(prob = 1, weight = 1)),
             error = function(e){print("Try different initial value of MLE")})
    s <- s + 0.1
  }
  return(error_model)
}

mis.fit.tmp1 <- mis_fit_func(LL.mis)

mis_search <- function(model.cum){
  tmp <- poisgeom.fit(model.cum, 0.1, 1, 0.1, 1, 0.1, 1, 0.01)
  tmp <- poisgeom.fit(model.cum, max(tmp[1]-0.02, 0.001), tmp[1]+0.02, max(tmp[2]-0.02, 0.001), min(tmp[2]+0.02, 1), 
                     max(tmp[3]-0.02, 0.001), min(tmp[3]+0.02, 1), 0.001)
  start <- c(max(tmp[1]-0.001, 0.0001), tmp[1]+0.001, max(tmp[2]-0.001, 0.0001), min(tmp[2]+0.001, 1), 
             max(tmp[3]-0.001, 0.0001), min(tmp[3]+0.001, 1))
  end <- poisgeom.fit(model.cum, start[1], start[2], start[3], start[4], start[5], start[6], 0.0001)
  i <- 0
  p_value <- end[4]
  while (Reduce("|", end %in% start) | i < 1000){
    last_end <- end
    start <- c(max(end[1]-0.001, 0.0001), end[1]+0.001, max(end[2]-0.001, 0.0001), end[2]+0.001,
               max(end[3]-0.001, 0.0001), min(end[3]+0.001, 1), max(end[4]-0.001, 0.0001), min(end[4]+0.001, 1))
    end <- poisgeom.fit(model.cum, start[1], start[2], start[3], start[4], start[5], start[6], 0.0001)
    if(p_value <= end[4]) {
      end <- last_end
      break
    }
    i <- i + 1
    p_value <- end[4]
  }
  return(end)
}

mis.fit.tmp2 <- mis_search(mis.cum)

# Choose the best fit between two methods
if (max(abs(ppoisgeommix(0:(length(mis.cum)-1), coef(mis.fit.tmp1)["lambda"], coef(mis.fit.tmp1)["prob"], coef(mis.fit.tmp1)["weight"]) - mis.cum)) <
      mis.fit.tmp2[4]) {
  mis.fit <- c(coef(mis.fit.tmp1)["lambda"], coef(mis.fit.tmp1)["prob"], coef(mis.fit.tmp1)["weight"])
} else {
  mis.fit <- mis.fit.tmp2
}

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
                              lambda = c(mis.fit[1], ins.fit[2], del.fit[2]),
                              k = c(0, ins.fit[1], del.fit[1]),
                              prob = c(mis.fit[2], ins.fit[3], del.fit[3]),
                              weight = c(mis.fit[3], ins.fit[4], del.fit[4]))

out_file <- paste(prefix, "_model_profile", sep="")
write.table(model_fit.table, out_file, row.names = FALSE, quote = FALSE, sep = "\t")
