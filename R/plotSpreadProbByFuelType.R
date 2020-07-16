#' Plot spread probability by fuel type and MDC
#'
#' @param spreadProbFuelType data.table. Columns as weather, classType and spreadProb.
#' @param typesOfFuel character. Types of fuel in the order of the formula terms.
#' @param coefToUse character. For labels' purpose. Either bestCoef or meanCoef.
#' @param covMinMax data.table. Covariate min a manx. Used to fix weather (MDC) axis.
#'
#' @return list of the original data table and the plot
#'
#' @author Tati Micheletti
#' @export
#' @importFrom ggplot2 ggplot aes geom_line labs
#' @importFrom viridis scale_color_viridis
#'
#' @rdname plotSpreadProbByFuelType

plotSpreadProbByFuelType <- function(spreadProbFuelType,
                                     typesOfFuel,
                                     coefToUse,
                                     covMinMax = NULL){
  # Fix the weather for the plotting!
if (!is.null(covMinMax)){
  correctedWeather <- rescaleKnown(x = spreadProbFuelType[["weather"]], 
                                   minNew = covMinMax[["weather"]][1]/1000, 
                                   maxNew = covMinMax[["weather"]][2]/1000, 
                                   minOrig = 0, maxOrig = 1)
} else {
  correctedWeather <- spreadProbFuelType[["weather"]]
}
  p1 <- ggplot(data = spreadProbFuelType, 
               aes(x = correctedWeather, y = spreadProb, group = classType)) +
    geom_line(aes(color = classType), size = 1.7) +
    scale_color_viridis(discrete = TRUE, option = "D",
                        name = "Fuel Type", 
                        labels = typesOfFuel) +
    labs(y = "spread probability", x = "MDC", 
         title = paste0("Coefficient: ", coefToUse))
  spreadProbFuelType <- list(DT = spreadProbFuelType, 
                                 plot = p1)
  return(spreadProbFuelType)
}