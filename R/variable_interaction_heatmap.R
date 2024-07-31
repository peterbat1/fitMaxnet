#' Plot Variable Interaction Heatmap
#'
#' @param feature_summary_table Data.frame. Data table produced by function \link{feature_summary}
#' @param num_repl Integer. Number of replicate models used to prepare feature_summary_table
#' @param factor_ordering Character array. Variable names in an order preferred for plotting. Default of NULL uses the
#' @param min_colour Colour name or hexadecimal colour value. Colour for the start of the colour gradient used to shade heatmap grid cells. Default is "white"
#' @param max_colour Colour name or hexadecimal colour value. Colour for the upper level of the colour gradient. Default is "blue"
#'
#' @details
#' The data provided in the parameter feature_summary_table is the number of replicate models in which each feature was retained (ie was given a non-zero coefficient).
#'
#' This does not represent feature importance. Some indirect inference about this may come from considering the importance of each of the two variables contributing to a feature as reported by \link{varImportance}.
#'
#' It might be possible to compute an estimate of feature importance using the same approach implemented in \link{varImportance}.
#'
#' This function is based on code found here: https://r-graph-gallery.com/79-levelplot-with-ggplot2.html
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{}
#'
interaction_heatmap <- function(feature_summary_table = "",
                                num_repl = NULL,
                                factor_ordering = NULL,
                                min_colour = "white",
                                max_colour = "blue")
{
  if (!("data.frame" %in% class(feature_summary_table))) stop("'feature_summary_table' must be a data.frame")

  if ((is.null(num_repl)) | (!is.numeric(num_repl))) stop("Non-numeric value found in parameter 'num_repl'")

  if (!is.integer(num_repl)) stop("Non-integer value found in parameter 'num_repl'")

  # Remove rows with just the variable which are always zero, and recover variable names
  var_rows <- which(!grepl(":", feature_summary_table[, 1]))
  var_names <- feature_summary_table[var_rows, 1]
  feature_summary_table <- feature_summary_table[-var_rows, ]

  # If no factor ordering has been provided, we default to var_names as
  # presented in the feature summary table
  if (is.null(factor_ordering)) factor_ordering<- var_names

  # Split factor names into variable names
  row_ind <- stringr::str_split_i(feature_summary_table[, 1], ":", 1)
  col_ind <- stringr::str_split_i(feature_summary_table[, 1], ":", 2)

  # Make a data.frame with values in ALL cells of the tile matrix. This takes
  # care of the previously present case of cells missing from upper but present
  # in lower half-matrices (caused by one or more variables which only appeared
  # in the second part of all interaction) terms. The full table will be pruned before plotting...
  full_plot_data <- data.frame(X = c(row_ind, col_ind),
                               Y = c(col_ind, row_ind),
                               Z = c(feature_summary_table[, "votes"]/num_repl, feature_summary_table[, "votes"]/num_repl))

  # Generate hashes to allow us to identify elements in the tile matrix to be
  # retained for plotting so that we only have the upper half matrix (including
  # the diagonal elements)
  hash_full_plot_data <- paste(full_plot_data$X, full_plot_data$Y, sep =  "_")

  # Make a template matrix to allow us to tag required matrix elements to be
  # retained by generating their row-column coordinates
  mat_template <- upper.tri(matrix(0, length(factor_ordering), length(factor_ordering)), diag = TRUE)

  # and generate the row-column coordinates of the elements
  low_ind <- which(mat_template, arr.ind = TRUE)

  # Make a hash coding of the elements to be retained but looking-up var names
  # using the row-column coordinates to be retained
  hash_2 <- paste(factor_ordering[low_ind[, 1]], factor_ordering[low_ind[, 2]], sep = "_")

  # Find the rows of the full data table which represent the items to be
  # retained for plotting; they have matching hashes...
  ii <- match(hash_2, hash_full_plot_data)

  # We can now filter the full data table to leave only those elements on and
  # above the diagonal of the tile matrix
  trimmed_plot_data <- full_plot_data[ii, ]

  # Make factors in preparation for plotting heatmap
    trimmed_plot_data$X <- factor(trimmed_plot_data$X, levels = factor_ordering)
    trimmed_plot_data$Y <- factor(trimmed_plot_data$Y, levels = factor_ordering)

  # Now we can create the desired heatmap plot...
  heatmap <- ggplot2::ggplot(trimmed_plot_data, aes(.data$X, .data$Y, fill = .data$Z)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="blue") +
    xlab("") +
    ylab("") +
    labs(fill = "Strength") +
    theme(axis.text.x = element_text(angle = 90, hjust =1 , vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 7))

  return(heatmap)
}
