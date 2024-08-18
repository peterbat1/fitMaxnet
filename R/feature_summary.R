#
# Peter D. Wilson
# Adjunct Fellow
# Dept. of Biological Sciences
# Faculty of Science and Engineering
# Macquarie University, Sydney, Australia
#
# 2021-03-10



#' Generate a feature summary table
#'
#' @param model_files Character array. Array of file names (including full path) to a set of maxnet models
#'
#' @return A data.frame with variables 'min", 'max' and 'votes' and row names representing variable and feature names, 'votes' is the number of replicate model fits which selected each feature.
#' @export
#'
#' @examples
#' \dontrun{}
#'
feature_summary <- function(model_files = "")
{

  # Check that all model files are accessible...
  if (model_files == "") stop("model-files is empty")

  if (length(model_files) < 2) stop("Two or more model filenames must be provided")

  if (!all(file.exists(model_files))) stop("One or more model files cannot be found")

  results_table <- NULL

  for (this_model_file in model_files)
  {

    load(this_model_file)

    vars_used <- sort(names(maxnet_model$varmax))

    feature_betas <- maxnet_model$betas

    tmp_names <- gsub("I(", "", names(feature_betas), fixed = TRUE)
    ii <- grep("^2)", tmp_names, fixed = TRUE)
    tmp_names[ii] <- gsub("^2)", "", tmp_names[ii], fixed = TRUE)
    tmp_names[ii] <- paste0(tmp_names[ii],":", tmp_names[ii])
    names(feature_betas) <- tmp_names

    beta_values <- data.frame(betaNames = tmp_names, beta = feature_betas)

    features_used <- gsub("I(", "", sort(names(maxnet_model$featuremaxs)), fixed = TRUE)
    ii <- grep("^2)", features_used, fixed = TRUE)
    features_used[ii] <- gsub("^2)", "", features_used[ii], fixed = TRUE)
    features_used[ii] <- paste0(features_used[ii],":",features_used[ii])

    features_used_df <- data.frame(betaNames = features_used)

    ans <- left_join(features_used_df, beta_values, "betaNames")
    ans <- ans[order(ans$betaNames), ]

    if (is.null(results_table))
    {
      results_table <- data.frame(ans[, 2])
      rownames(results_table) <- ans[, 1]
    }
    else
      results_table <- cbind(results_table,  ans[, 2])
  }

  colnames(results_table) <- paste0("R_", stringr::str_pad(as.character(1:length(model_files)), side = "left", width = 2, pad = "0"))

  results_table_zero <- results_table
  results_table_zero[which(is.na(results_table_zero), arr.ind = TRUE)] <- 0

  min_and_max <- t(Rfast::rowMinsMaxs(as.matrix(results_table_zero)))
  #rownames(min_and_max) <- rownames(results_table)

  results_table_summary <- data.frame(feature = rownames(results_table),
                                      min_and_max,
                                      votes = apply(results_table, 1, function(x) {sum(!is.na(x))}))

  return(results_table_summary)
}
