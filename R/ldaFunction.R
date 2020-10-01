ldaFunction <- function (data, lfk, rfk, min_cl, ncl) {
  # test 1000 samples for contrast within classes per feature
  # and that there is at least a minimum number of samples per class
  for (j in seq_along(1:1000)) {
    rand_s <- sample(c(1:lfk), rfk, replace = TRUE)
    if (!contastWithinClassesOrFewPerClass(data, rand_s, min_cl, ncl)) {
      break
    }
  }
  # lda with rfk number of samples
  lda.fit <- lda(class ~ ., data = data, subset = rand_s)
  # coefficients that transform observations to discriminants
  w <- lda.fit$scaling[, 1]
  # scaling of lda coefficients
  w.unit <- w / sqrt(sum(w ^ 2))
  sub_d <- data[rand_s, ]
  ss <- sub_d[, -match("class", colnames(sub_d))]
  xy.matrix <- as.matrix(ss)
  # the original matrix is transformed
  LD <- xy.matrix %*% w.unit
  # effect size is calculated as difference between averaged disciminants
  # of two classes
  effect_size <-
    abs(mean(LD[sub_d[, "class"] == 0]) - mean(LD[sub_d[, "class"] == 1]))
  # scaling lda coefficients by the efect size
  scal <- w.unit * effect_size
  # mean count values per fclass per feature
  rres <- lda.fit$means
  rowns <- rownames(rres)
  lenc <- length(colnames(rres))
  
  coeff <- vector("numeric", length(scal))
  for (v in seq_along(scal)) {
    if (is.na(scal[v]) != TRUE) {
      coeff[v] <- abs(scal[v])
    } else{
      coeff[v] <- 0
    }
    
  }
  # count value differences between means of two classes for each feature
  lda.means.diff <- (lda.fit$means[2, ] - lda.fit$means[1, ])
  # difference between a feature's class means and effect size adjusted lda coefficient
  # are averaged for each feature
  (lda.means.diff + coeff) / 2
}