# ensures that more than half of the values in each for each feature are unique
# if that is not the case then a count value is altered by adding it to a small value
# generated via normal distribution with mean=0 and sd=5% of the count value
createUniqueValues <- function(expr_sub_t_df, groups){
  df <- expr_sub_t_df
  for (i in seq_along(c(1:(ncol(df) - 1)))) {
    for (j in seq_along(groups)) {
      equality = df[df[, "class"] == groups[j], i]
      if (length(unique(equality)) > max(length(equality) * 0.5, 4)) {
        next
      } else{
        for (k in seq_along(equality)) {
          equality[k] = abs(equality[k] + rnorm(
            1,
            mean = 0,
            sd = max(equality[k] * 0.05, 0.01)
          ))
        }
        
      }
      df[df[, "class"] == groups[j], i] = equality
    }
  }
  expr_sub_t_df <- df
}
