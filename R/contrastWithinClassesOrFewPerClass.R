contastWithinClassesOrFewPerClass <-
  function(expr_sub_t_df, rand_s, min_cl, ncl) {
    cols <- expr_sub_t_df[rand_s, ]
    cls <- expr_sub_t_df$class[rand_s]
    # if the number of classes is less than the actual number (typically two)
    # of classes in the dataframe then return TRUE
    if (length(unique(cls)) < ncl) {
      return (TRUE)
    }
    # detect if for each class there are not fewer than the minimum (min_cl) number of samples
    if (TRUE %in% c(table(cls) < min_cl)) {
      return (TRUE)
    }
    # separate the randomly selected samples (cols) into a list of the two classes
    drops <- c("class")
    by_class <-
      lapply(seq_along(groups), function(x) {
        cols[cols[, "class"] == groups[x], !(names(cols) %in% drops)]
      })
    
    # makes sure that within each class all features have at least min_cl unique count values
    for (i in seq_along(groups)) {
      unique_counts_per_microb = apply(by_class[[i]], 2, function(x) {
        length(unique(x))
      })
      if ((TRUE %in% c(unique_counts_per_microb <= min_cl) &
           min_cl > 1) |
          (min_cl == 1 & (TRUE %in% c(unique_counts_per_microb <= 1)))) {
        return (TRUE)
      }
    }
    return (FALSE)
    
  }