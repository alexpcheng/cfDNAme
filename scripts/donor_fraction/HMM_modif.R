correctReadcount2 <- function(x, mappability = 0.9, samplesize = 50000,
                             verbose = TRUE) {
  if (length(x$reads) == 0 | length(x$gc) == 0 | length(x$map) == 0) {
    stop("Missing one of required columns: reads, gc, map")
  }
  
  if(verbose) { message("Applying filter on data...") }
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid], prob = c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid], prob = c(doutlier, 1 - doutlier),
                     na.rm = TRUE)
  x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
            x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  
  if (verbose) { message("Correcting for GC bias...") }
  set <- which(x$ideal)
  #select <- sample(set, min(length(set), samplesize))
  select <- set
  rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads / predict(final, x$gc)
  
  if (verbose) { message("Correcting for mappability bias...") }
  coutlier <- 0.01
  range <- quantile(x$cor.gc[which(x$valid)],
                    prob = c(0, 1 - coutlier), na.rm = TRUE)
  set <- which(x$cor.gc < range[2])
  #select <- sample(set, min(length(set), samplesize))
  select <- set
  final = approxfun(lowess(x$map[select], x$cor.gc[select]))
  x$cor.map <- x$cor.gc / final(x$map)
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(x)
}
