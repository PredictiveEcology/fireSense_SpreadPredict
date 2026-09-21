defineModule(sim, list(
  name = "fireSense_SpreadPredict",
  description = "Predicts a surface of fire spread probilities using a model fitted with fireSense_SpreadFit.",
  keywords = c("fire spread", "fireSense", "predict"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Tati", "Micheletti", email = "tati.micheletti@gmail.com", role = "aut"),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = "aut"),
    person("Alex M.", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = list(fireSense_SpreadPredict = "1.0.0.9001", SpaDES.core = "0.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "Matrix", "methods", "terra", "SpaDES.core (>=3.0.4)", "stats",
                  "ggplot2", "viridis",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.2.3.9029)"),
  parameters = bindrows(
    defineParameter(name = "lowerSpreadProb", class = "numeric", default = 0.13,
                    desc = "Lower asymptote of the 2- and 3-parameter logistic."),
    defineParameter("maxFireSpread", "numeric", default = 0.28,
                    desc = paste("Upper limit on `spreadProb` used when fitting. Here it is only checked",
                                 "to be the same in every module that defines it.")),
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "Time of the first prediction."),
    defineParameter(name = ".runInterval", class = "numeric", default = 1,
                    desc = "Interval between predictions, in years. `NA` predicts once."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA,
                    desc = "Time of the `save` event, which does nothing. `NA` means never."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated?",
                          "This is generally intended for data-type modules, where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "covMinMax_spread", objectClass = "data.table",
                 desc = paste("Minimum and maximum (2 rows) of each covariate in the fitting data,",
                              "used to rescale the covariates as in `fireSense_SpreadFit`.")),
    expectsInput(objectName = "fireSense_SpreadCovariates", objectClass = "data.table",
                 desc = paste("This year's covariates, from `fireSense_dataPrepPredict`.",
                              "`pixelID` is the cell index of `flammableRTM`.")),
    expectsInput(objectName = "flammableRTM", objectClass = "SpatRaster", sourceURL = NA,
                 desc = "Binary raster, 1 where the pixel is flammable. Template for `fireSense_SpreadPredicted`.")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "fireSense_SpreadPredicted", objectClass = "SpatRaster",
                  desc = "Spread probability of each flammable pixel, this year.")
  ))
)

#' Event dispatcher
#'
#' Events: `init`, `run` (predict, repeated every `.runInterval`), `save` (does nothing).
#'
#' @param sim A `simList`.
#' @param eventTime Time of the event.
#' @param eventType Name of the event.
#' @param debug Not used.
#'
#' @return The `simList`, invisibly.
doEvent.fireSense_SpreadPredict <- function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {

      SpaDES.core::paramCheckOtherMods(sim, paramToCheck = "maxFireSpread", moduleToUse = "all")
      SpaDES.core::paramCheckOtherMods(sim, paramToCheck = "maxFireSpread", moduleToUse = "all")
      
      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run",
                           eventPriority = 5.12)

      if (!is.na(P(sim)$.saveInitialTime)) {
        sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
      }
    },
    run = {
      sim <- spreadPredictRun(sim)

      if (!is.na(P(sim)$.runInterval)) {
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "run",
                             eventPriority = 5.12)
      }
    },
    save = {
      message("fireSense_SpreadPredict: the save event does nothing")
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'",
                  sep = ""
    ))
  )

  invisible(sim)
}

#' Predict this year's spread probability
#'
#' Rescales `sim$fireSense_SpreadCovariates` with `sim$covMinMax_spread`, computes the spread
#' probability for each parameter set (row) in `sim$studyAreaWithSpreadParams$params[[1]]`,
#' and writes the mean over parameter sets to `sim$fireSense_SpreadPredicted`.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly.
spreadPredictRun <- function(sim) {
  moduleName <- current(sim)$moduleName

  fireSense_SpreadCovariates <- copy(sim$fireSense_SpreadCovariates)

  ## Fuel biomass arrives logged (fireSenseUtils::logMinB()). A fit made on LINEAR fuel biomass has
  ## fireSenseUtils::fuelLinearRange, c(0, 1e4), as that covariate's covMinMax_spread, and its
  ## coefficients only mean anything for biomass / 1e4: undo the log with the function the fit used.
  ## A fit made on the log scale has the log range there, and its covariates are left as they are.
  for (cn in intersect(names(sim$covMinMax_spread), names(fireSense_SpreadCovariates))) {
    if (fireSenseUtils::isLinearFuelRange(sim$covMinMax_spread[[cn]]))
      fireSense_SpreadCovariates[[cn]] <- fireSenseUtils::fuelLogToLinear(fireSense_SpreadCovariates[[cn]])
  }

  # Load inputs in the data container
  mod_env <- new.env(parent = globalenv())
  list2env(fireSense_SpreadCovariates, envir = mod_env)
  ## In case there is a response in the formula remove it

  terms <- as.formula(sim$fireSense_spreadFormula) %>%
    terms.formula() %>%
    delete.response()

  formula <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula)

  missing <- !allxy %in% ls(mod_env, all.names = TRUE)
  if (s <- sum(missing)) {
    stop(
      moduleName, "> '", allxy[missing][1L], "'",
      if (s > 1) paste0(" (and ", s - 1L, " other", if (s > 2) "s", ")"),
      " not found in data objects."
    )
  }

  # integers x 1000, the form `spreadProbFromIntegerCovs` expects
  shortAnnDTx1000 <- toX1000(list(fireSense_SpreadCovariates))[[1]] |> setDT()
  colsToUse <- setdiff(names(fireSense_SpreadCovariates), "pixelID")

  # Without fitted parameters there is nothing to predict from; say so instead of
  # dying in rowMeans() on an empty matrix (which is what an unfitted ELF produced
  # when fireSense_SpreadFit had not run first). This must come before anything
  # indexes `params[[1]]`: with zero rows that fails first, "subscript out of bounds".
  nPar <- tryCatch(NROW(sim$studyAreaWithSpreadParams$params[[1]]), error = function(e) 0L)
  if (NROW(sim$studyAreaWithSpreadParams) == 0L || is.null(nPar) || nPar == 0L)
    stop("fireSense_SpreadPredict: sim$studyAreaWithSpreadParams holds no fitted spread ",
         "parameters for this run (", if (!is.null(sim$.runName)) sim$.runName else "unknown",
         "). Either fireSense_SpreadFit has not run yet -- its `run` event must precede this ",
         "module's -- or the shared ledger has no row for this polygon.", call. = FALSE)

  logisticPars <- sim$studyAreaWithSpreadParams$params[[1]]
  
  shortAnnDT <-
    spreadProbFromIntegerCovs(shortAnnDTx1000 = shortAnnDTx1000,
                              yr = time(sim),
                              covMinMax = sim$covMinMax_spread,
                              mutuallyExclusive = NULL, # alraedy done in dataPrepPredict
                              colsToUse = colsToUse,
                              doAssertions = FALSE,
                              logisticPars = logisticPars,
                              maxFireSpread = Par$maxFireSpread
                              )

  parsModel <- length(colsToUse)
  mat <- as.matrix(shortAnnDT[, ..colsToUse])

  # for replicate "best" params from DEoptim
  spreadProbList <- purrr::pmap(.l = list(ind = seq(NROW(sim$studyAreaWithSpreadParams$params[[1]]))),
                     sa = sim$studyAreaWithSpreadParams, function(ind, sa) {
                       par <- sa$params[[1]][ind,] |> as.vector() |> unlist()
                       covPars <- intersect(names(par), colsToUse)
                       covPars <- par[covPars]
                       logisticPars <- par[setdiff(names(par), names(covPars))]
                       # Make sure the order is correct in the matrix
                       matching <- intersect(names(covPars), colnames(mat))
                       missingCovs <- setdiff(colnames(mat), names(covPars))
                       if (length(missingCovs))
                         warning("There are covariates in the sim$fireSense_SpreadCovariates: \n",
                              paste0(missingCovs, collapse = ", "),
                              "\n...that are not in the sim$studyAreaWithSpreadParams")
                       mat <- mat[, matching]

                       logisticAll(logisticPars,
                                   mat, covPars, P(sim)$lowerSpreadProb)
                     })
  spreadProbMat <- do.call(cbind, spreadProbList)
  
  set(shortAnnDT, NULL, "spreadProb", rowMeans(spreadProbMat))

  # Return to raster format
  sim$fireSense_SpreadPredicted <- rast(sim$flammableRTM) ## use flammableRTM as template
  sim$fireSense_SpreadPredicted[shortAnnDT$pixelID] <- shortAnnDT$spreadProb

  invisible(sim)
}
