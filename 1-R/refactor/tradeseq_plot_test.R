
# .plotSmoothers <- function(model, nPoints = 100, lwd = 2, size = 2/3,
#                            xlab = "Pseudotime",
#                            ylab = "Log(expression + 1)",
#                            border = FALSE,
#                            alpha = 1,
#                            sample = 1)
# {
  ## test variable is a single gam structure of the gamList in the bulk timecourse example 
  data <- test$model
  y <- data$y
  
  #construct time variable based on cell assignments.
  nCurves <- length(test$smooth)
  col <- timeAll <- rep(0, nrow(data))
  for (jj in seq_len(nCurves)) {
    for (ii in seq_len(nrow(data))) {
      if (data[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- data[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }
  
  # plot raw data
  df <- data.frame("time" = timeAll,
                   "gene_count" = y,
                   "lineage" = as.character(col))
  rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
  df <- df[rows, ]
  p <- ggplot(df, aes(x = time, y = log1p(gene_count), col = lineage)) +
    geom_point(size = size) +
    labs(x = xlab, y = ylab) +
    theme_classic() +
    scale_color_viridis_d(alpha = alpha)
  
  
  # predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(model$model, jj, nPoints = nPoints)
    yhat <- predict(model, newdata = df, type = "response")
    if (border) {
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd + 1, colour = "white") +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd)
    } else {
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd)
    }
    
  }
  return(p)
}