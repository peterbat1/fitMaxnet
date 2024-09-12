

#' Make a smoothed calibration plot
#'
#' Make a smoothed calibration plot for a presence-only model using the method by Phillips and Elith
#'
#' @param pred Data frame. Predicted probability of presence at occurrence and background locations
#' @param negrug Numeric. Model prediction of presence at background points.
#' @param posrug Numeric. Model prediction of presence at occurrence points.
#' @param ideal Function.
#' @param ylim Numeric. Two-element vector giving the x-axis for the plot; default is c(0, 1).
#' @param xlim Numeric. Two-element vector giving the y-axis for the plot; default is c(0, 1).
#' @param capuci Logical. Should the upper confidence limit of the smoothed probability of presence by clamped to 1?
#' @param xlabel Character. Label to be applied to the x-axis.
#' @param ylabel Character. Label to be applied to the y-axis.
#' @param filename Character. Full path to the file into which the graph will be plotted. Default of NULL plots to the standard graphics device.
#' @param title Character. Title for the plot. Default is a generic "Calibration plot".
#'
#' @details {
#' This is an adaptation of code published by Phillips and Elith (2010. POC plots: calibrating species distribution models with presence-only data. Ecology 91:2476–2484).
#'
#' Although it is possible to call this function with user-generated parameter values, it is designed to be called by the function \link{POCplot}.
#'
#' This implementation uses ggplot2 to render the plot.
#'
#' }
#'
#' @return
#' A ggplot2 object
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' this_plot <- calibplot(prediction_data)
#' }

calibplot <- function(pred,
                       negrug,
                       posrug,
                       ideal,
                       ylim = c(0, 1),
                       xlim = c(0, 1),
                       capuci = TRUE,
                       xlabel = "Predicted probability of presence",
                       ylabel = "Probability of presence",
                       filename = NULL,
                       title = "Calibration plot")
{
  #if (!is.null(filename)) png(filename)
  ylow <- pred$y - 2 * pred$se
  ylow[ylow < 0] <- 0
  yhigh <- pred$y + 2 * pred$se

  if (capuci) yhigh[yhigh > 1] <- 1

  plotData <- data.frame(x = pred$x,
                         y = pred$y,
                         ylow = ylow,
                         yhigh = yhigh)

  p <- ggplot2::ggplot(plotData, aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(colour = "blue") +
    ggplot2::geom_segment(x = 0, y = 0, xend = 1, yend = 1, colour = "grey30", linetype = 3) +
    ggplot2::annotate("line", x = plotData[, "x"], y = plotData[, "ylow"], linetype = 1, size = 1, colour = "orange") +

    ggplot2::annotate("line", x = plotData[, "x"], y = plotData[, "yhigh"], linetype = 1, size = 1, colour = "orange") +
    ggplot2::ggtitle(title) +
    xlim(xlim) +
    ylim(ylim) +
    xlab(xlabel) +
    ylab(ylabel) +
    annotate("rug", x = negrug, sides = "t", colour = "black") +
    annotate("rug", x = posrug, sides = "b", colour = "orange") +
    theme(plot.title = element_text(face = "bold", size = 14))

  if (!is.null(filename))
  {
    grDevices::png(filename = filename)
    print(p)
    grDevices::dev.off()
  }
  else
    print(p)

}


####################################################################
#' Smoothed density distribution
#'
#' A function to generate a smoothed density function of observed probability of presence
#'
#' @param pred Numeric. Vector of model predictions for combined occurrence and background points.
#' @param res Numeric. Vector expected presence and background values (referred in the function as 'res' = 'response variable'); 0 = expected background, 1 = expected presence.
#' @param smoothing_df Numeric. Degrees of freedom for a natural cubic spline fitted to model predictions at combined occurrence and background locations.
#'
#' @details {
#' This is an adaptation of code published by Phillips and Elith (2010. POC plots: calibrating species distribution models with presence-only data. Ecology 91:2476–2484).
#'
#' It is normally called by the function \link{POCplot}.
#'
#' The original script used a default smoothing degrees of freedom (\emph{smoothing_df}) of 6, but trials during development of this implementation suggest that a better default value is 4.
#' }
#'
#' @return
#' A data.frame with three columns:
#' \describe{
#' \item{x}{The values used to fit the smotthing function}
#' \item{y}{The y-values of the computed smoothing function}
#' \item{se}{The standard error of y-values of the computed smoothing function
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{}
smoothdist <- function(pred,
                       res,
                       smoothing_df = 4)
{

  gam1 <- glm(res ~ splines::ns(pred, df = smoothing_df), weights = rep(1, length(pred)), family = binomial)
  x <- seq(min(pred), max(pred), length = 512)
  y <- predict(gam1, newdata = data.frame(pred = x), se.fit = TRUE,
               type = "response")
  data.frame(x = x, y = y$fit, se = y$se.fit)
}



#' Presence-only smoothed calibration plot
#'
#' Produce a smoothed calibration plot for a presence-only model by the method of Phillips and Elith
#'
#' @param pred Numeric. Vector of model predictions at occurrence points.
#' @param back Numeric. Vector of model predictions at background points.
#' @param linearize Logical. Should the logistic transform be applied to model predictions? The default is FALSE; see Details.
#' @param capUpperValues Logical. Should computed values in smoother fit be clamped to 1?
#' @param title Character. Title to appear in the plot. Default is a generic "Calibration plot".
#' @param filename Character. Full path to a file into which the calibration plot will be written. Default (NULL) will cause the plot to appear on the standard graphics device.
#' @details {
#' This is an adaptation of code published by Phillips and Elith (2010. POC plots: calibrating species distribution models with presence-only data. Ecology 91:2476–2484).
#'
#' The original script set the default value of the parameter \emph{linearize} to TRUE which applies a logistic transformation to the y-values in the POCplot. In development trials, this frequently produced numerical failures when values approaching 1 where transformed. It is perhaps intended to transform 'raw' values supplied by fitted ENMs to the logistic scale. However, the most frequent output from MaxEnt and MaxEnt-like models is already transformed onto a logistic or complementary log-log ('cloglog') scale. Inadvertent double application of a linearising transformation will cause numeric failures and severely distorted calibration plots. So, the default for \emph{linearize} is set to FALSE.
#'
#' }
#'
#' @return A named list with the follow elements:
#' \describe{
#' \item{predd}{Data frame of results from the smoothing process; columns are: x = predicted probability of presence, y = (modelled) probability of presence, se = std error of the y-value at each x-value}
#' \item{mse}{A numeric object giving the Mean Squared Error between the line of equality between probability of presence (y) and predicted probability of presence (x) and the y-values}
#' }
#' @export
#'
#' @examples
#' \dontrun{}
POCplot <- function(pred,
                    back,
                    linearize = FALSE,
                    capUpperValues = TRUE,
                    title = "Calibration plot",
                    filename = NULL)
{
  isPresence <- c(rep(1, length(pred)), rep(0, length(back)))
  predd <- smoothdist(c(pred, back), isPresence)
  c <- mean(back)*length(back)/length(pred)
  if (linearize)
  {
    fun <- function(x,y) c*y / (1 - y)
    predd$y <- mapply(fun, predd$x, predd$y)
    predd$se <- mapply(fun, predd$x, predd$se)
    ideal <- function(x) x
    ylab <- "Relative probability of presence"
  }
  else
  {
    ideal <- function(x) x / (x + c)
    ylab <- "Probability of presence"
  }

  calibplot(predd, negrug = back, posrug = pred, ideal = ideal, ylabel = ylab,
             capuci = capUpperValues, title = title, filename = filename)
  invisible(list(predd = predd, mse = sum((predd$y - predd$x)^2)/nrow(predd)))
}

