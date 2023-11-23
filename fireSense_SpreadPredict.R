defineModule(sim, list(
  name = "fireSense_SpreadPredict",
  description = "Predicts a surface of fire spread probilities using a model fitted with fireSense_SpreadFit.",
  keywords = c("fire spread", "fireSense", "predict"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Tati", "Michelleti", email = "tati.micheletti@gmail.com", role = "aut"),
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = "aut")
  ),
  childModules = character(),
  version = list(fireSense_SpreadPredict = "0.0.1", SpaDES.core = "0.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "Matrix", "methods", "terra", "SpaDES.core", "stats",
                  "ggplot2", "viridis",
                  "PredictiveEcology/fireSenseUtils@development"),
  parameters = bindrows(
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start time of the simulation."),
    defineParameter(name = ".runInterval", class = "numeric", default = 1,
                    desc = paste("optional. Interval between two runs of this module, expressed in units of",
                                 "simulation time.Defaults to 1 year.")),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA,
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInterval", class = "numeric", default = NA,
                    desc = "optional. Interval between save events."),
    defineParameter(name = "climCol", class = "character", default = "MDC", min = NA, max = NA,
                    desc = "the name of the climate covariate in sim$fireSense_spreadCovariates"),
    defineParameter(name = "coefToUse", class = "character", default = "meanCoef",
                    desc = paste("Which coefficient to use to predict?",
                                 "The best coefficient (bestCoef) from DEOPtim or ",
                                 "the average (meanCoef; default).")),
    defineParameter(name = "lowerSpreadProb", class = "numeric", default = 0.13,
                    desc = "Lower spread probability"),
    defineParameter(name = "mutuallyExclusiveCols", "list", default = list("youngAge" = "vegPC"), NA, NA,
                    desc = paste("a named list, where the name of the list must be a covariate in the data.table.",
                                 "Covariates matching the values in each list element will be set to 0.",
                                 "List content should be a grep regex.")),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated?",
                          "This is generally intended for data-type modules, where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "covMinMax_spread", objectClass = "data.table",
                 desc = "range used to rescale coefficients during spreadFit"),
    expectsInput(objectName = "fireSense_SpreadCovariates", objectClass = "data.table",
                 desc = "data.table of covariates with pixelID column corresponding to flammableRTM index."),
    expectsInput(objectName = "fireSense_SpreadFitted", objectClass = "fireSense_SpreadFit",
                 desc = "An object of class 'fireSense_SpreadFit' created by the fireSense_SpreadFit module."),
    expectsInput(objectName = "flammableRTM", objectClass = "SpatRaster", sourceURL = NA,
                 desc = "RTM with nonflammable pixels coded as 0 and flammable as 1.")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "fireSense_SpreadPredicted", objectClass = "SpatRaster",
                  desc = "A raster layer of spread probabilities")
  ))
)
## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_SpreadPredict <- function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {

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
      sim <- spreadPredictSave(sim)

      if (!is.na(P(sim)$.saveInterval)) {
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, moduleName, "save", .last())
      }
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
      "' in module '", current(sim)[1, "moduleName", with = FALSE], "'",
      sep = ""
    ))
  )

  invisible(sim)
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initialization;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

spreadPredictRun <- function(sim) {
  moduleName <- current(sim)$moduleName

  fireSense_SpreadCovariates <- copy(sim$fireSense_SpreadCovariates)

  if (!is(sim$fireSense_SpreadFitted, "fireSense_SpreadFit")) {
    stop(moduleName, "> '", sim$fireSense_spreadFitted, "' should be of class 'fireSense_SpreadFit")
  }

  # Load inputs in the data container
  mod_env <- new.env(parent = globalenv())
  list2env(fireSense_SpreadCovariates, env = mod_env)
  ## In case there is a response in the formula remove it
  terms <- as.formula(sim$fireSense_SpreadFitted$formula) %>%
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

  ###################################################
  # Convert stacks to lists of data.table objects --> much more compact
  ###################################################
  # First for stacks that are "annual"

  # # Rescale to numerics and /1000
  if (!is.null(sim$covMinMax_spread)) {
    for (cn in names(sim$covMinMax_spread)) {
      set(
        fireSense_SpreadCovariates, NULL, cn,
        rescaleKnown2(x = fireSense_SpreadCovariates[[cn]],
                      minNew = 0,
                      maxNew = 1000,
                      minOrig = sim$covMinMax_spread[[cn]][1],
                      maxOrig = sim$covMinMax_spread[[cn]][2])
      )
    }
  }

  if (!is.null(P(sim)$mutuallyExclusiveCols)) {
    fireSense_SpreadCovariates <- makeMutuallyExclusive(
      dt = fireSense_SpreadCovariates,
      mutuallyExclusiveCols = P(sim)$mutuallyExclusiveCols
    )
  }

  colsToUse <- setdiff(names(fireSense_SpreadCovariates), "pixelID")
  parsModel <- length(colsToUse)

  
  par <- sim$fireSense_SpreadFitted[[P(sim)$coefToUse]]
  mat <- as.matrix(fireSense_SpreadCovariates[, ..colsToUse])/1000 # Divide by 1000 for the model prediction

  # matrix multiplication
  covPars <- tail(x = par, n = parsModel)
  logisticPars <- head(x = par, n = length(par) - parsModel)
  # Make sure the order is correct in the matrix
  matching <- match(names(covPars), colnames(mat))
  mat <- mat[, matching]
  if (length(logisticPars) == 4) {
    set(fireSense_SpreadCovariates, NULL, "spreadProb", logistic4p(mat %*% covPars, logisticPars))
  } else if (length(logisticPars) == 3) {
    set(fireSense_SpreadCovariates, NULL, "spreadProb", logistic3p(mat %*% covPars, logisticPars,
                                                                   par1 = P(sim)$lowerSpreadProb))
  } else if (length(logisticPars) == 2) {
    set(fireSense_SpreadCovariates, NULL, "spreadProb", logistic2p(mat %*% covPars, logisticPars,
                                                                   par1 = P(sim)$lowerSpreadProb))
  }

  # Return to raster format
  sim$fireSense_SpreadPredicted <- rast(sim$flammableRTM) ## use flammableRTM as template
  ## Need to track what is happening with missing pixels
  if (FALSE) {
    nFlam <- sum(getValues(sim$flammableRTM), na.rm = TRUE)
    nLand <- nrow(sim$landcoverDT)
    nSpread <- nrow(sim$fireSense_SpreadCovariates)
    nIg <- nrow(sim$fireSense_IgnitionAndEscapeCovariates)
  }
  sim$fireSense_SpreadPredicted[fireSense_SpreadCovariates$pixelID] <- fireSense_SpreadCovariates$spreadProb

  invisible(sim)
}
