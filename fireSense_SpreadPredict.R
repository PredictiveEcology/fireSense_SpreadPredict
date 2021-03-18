defineModule(sim, list(
  name = "fireSense_SpreadPredict",
  description = "Predicts a surface of fire spread probilities using a model fitted with fireSense_SpreadFit.",
  keywords = c("fire spread", "fireSense", "predict"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut", "cre")),
    person("Tati", "Michelleti", email = "tati.micheletti@gmail.com", role = "aut"),
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = "aut")
  ),
  childModules = character(),
  version = list(fireSense_SpreadPredict = "0.0.1", SpaDES.core = "0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "Matrix", "methods", "raster", "SpaDES.core", "stats",
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
    defineParameter(name = "lowerSpreadProb", class = "numeric", default = 0.13,
                    desc = "Lower spread probability"),
    defineParameter(name = "coefToUse", class = "character",default = "bestCoef", # meanCoef
                    desc = paste0("Which coefficient to use to predict? The best coefficient (bestCoef) from DEOPtim or ",
                                  "the average (meanCoef). default is bestCoef")),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated?",
                          "This is generally intended for data-type modules, where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "covMinMax", objectClass = "data.table",
                 description = "range used to rescale coefficients during spreadFit"),
    expectsInput(objectName = "fireSense_SpreadCovariates", objectClass = "data.table",
                 desc = "data.table of covariates with pixelID column corresponding to flammableRTM index."),
    expectsInput(objectName = "fireSense_SpreadFitted", objectClass = "fireSense_SpreadFit",
                 desc = "An object of class 'fireSense_SpreadFit' created by the fireSense_SpreadFit module."),
    expectsInput(objectName = "flammableRTM", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "RTM with nonflammable pixels coded as 0 and flammable as 1.")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "fireSense_SpreadProbRaster", objectClass = "RasterLayer",
                  desc = "A raster layer of spread probabilities"),
    createsOutput(objectName = 'spreadPredictedProbability' objectClass = "list",
                  desc = "list of annual spread probabilities")
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_SpreadPredict <- function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {

      sim$spreadPredictedProbability <- list()

      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run")

      if (!is.na(P(sim)$.saveInitialTime)) {
        sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
      }
    },
    run = {
      sim <- spreadPredictRun(sim)

      if (!is.na(P(sim)$.runInterval)) {
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "run")
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

  if (!is(sim$fireSense_SpreadFitted, "fireSense_SpreadFit")) {
    stop(moduleName, "> '", sim$fireSense_spreadFitted, "' should be of class 'fireSense_SpreadFit")
  }

  # Load inputs in the data container
  # list2env(as.list(envir(sim)), envir = mod)

  mod_env <- new.env(parent = glo)
  list2env(fireSense_SpreadCovariates, env = mod_env)
  ## In case there is a response in the formula remove it
  terms <- sim$fireSense_SpreadFitted$formula %>%
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
  whNotNA <- which(!is.na(sim$flammableRTM[]))
  hash <- fastdigest(sim$fireSense_SpreadCovariates)
  fireSenseDataDTx1000  <- annualStackToDTx1000(annualStack = sim$fireSense_SpreadCovariates,
                                whNotNA = whNotNA,
                                .fastHash = hash,
                                timeSim = paste0("year", time(sim)),
                                omitArgs = c("annualStack",
                                      "rasterToMatch"))
  # # Rescale to numerics and /1000
  if (!is.null(sim$covMinMax)) {
    for (cn in colnames(sim$covMinMax)) {
      if (cn != "weather"){
        set(fireSenseDataDTx1000, NULL, cn,
            rescaleKnown(x = fireSenseDataDTx1000[[cn]], minNew = 0, maxNew = 1000,
                         minOrig = sim$covMinMax[[cn]][1], maxOrig = sim$covMinMax[[cn]][2]))
      } else {
        set(fireSenseDataDTx1000, NULL, cn,
            rescaleKnown(x = fireSenseDataDTx1000[[cn]], minNew = 0,
                         maxNew = 1000*(max(fireSenseDataDTx1000[[cn]])/sim$covMinMax[[cn]][2]),
                         minOrig = sim$covMinMax[[cn]][1], maxOrig = sim$covMinMax[[cn]][2]))
      }
    }
  } else {
    fireSenseDataDTx1000 <- fireSenseDataDTx1000
  }
  colsToUse <- setdiff(names(sim$fireSense_SpreadCovariates), 'pixelID')
  parsModel <- length(colsToUse)

  par <- sim$fireSense_SpreadFitted_year2011[[P(sim)$coefToUse]]
  mat <- as.matrix(fireSenseDataDTx1000[, ..colsToUse])/1000 # Divide by 1000 for the model prediction

  # matrix multiplication
  covPars <- tail(x = par, n = parsModel)
  logisticPars <- head(x = par, n = length(par) - parsModel)
  # Make sure the order is correct in the matrix
  matching <- match(names(covPars), colnames(mat))
  mat <- mat[, matching]

  if (length(logisticPars) == 4) {
    set(fireSenseDataDTx1000, NULL, "spreadProb", logistic4p(mat %*% covPars, logisticPars))
  } else if (length(logisticPars) == 3) {
    set(fireSenseDataDTx1000, NULL, "spreadProb", logistic3p(mat %*% covPars, logisticPars,
                                                             par1 = P(sim)$lowerSpreadProb))
  } else if (length(logisticPars) == 2) {
    set(fireSenseDataDTx1000, NULL, "spreadProb", logistic2p(mat %*% covPars, logisticPars,
                                                             par1 = P(sim)$lowerSpreadProb))
  }

  if (time(sim) == start(sim)) {
    # We want a full distribution of the spread prob for each fuel type for the
    # whole range of MDC
    weatherValues <- sort(unique(mat[colnames(mat) == "weather"]))
    thinned <- weatherValues[seq.int(1L, length(weatherValues), 10L)] # thin as we have 30k vals
    # Spread probability of each fuel type.
    # I need the whole thinned vector repeated the n times the number of params
    # (length(covPars)-1)
    thinnedExp <- data.table(weather = rep(thinned, times = length(covPars) - 1))
    matExp <- data.table(matrix(rep(as.numeric(Matrix::diag(length(covPars) - 1)),
                                    each = length(thinned)),
                                ncol = length(covPars) - 1))
    names(matExp) <- names(covPars)[names(covPars) != "weather"]
    m <- as.matrix(cbind(thinnedExp, matExp))
    # Now in data.table format so I can add spreadProb and
    # make the plots
    sim$spreadProbFuelType <- data.table(m)
    sim$spreadProbFuelType$classType <- rep(names(sim$spreadProbFuelType)[-1],
                                            each = length(thinned))
    # Now I calculate the spreadProb
    if (length(logisticPars) == 4) {
      set(sim$spreadProbFuelType, NULL, "spreadProb", logistic4p(m %*% covPars, logisticPars))
    } else if (length(logisticPars) == 3) {
      set(sim$spreadProbFuelType, NULL, "spreadProb", logistic3p(m %*% covPars, logisticPars,
                                                                 par1 = P(sim)$lowerSpreadProb))
    } else if (length(logisticPars) == 2) {
      set(sim$spreadProbFuelType, NULL, "spreadProb", logistic2p(m %*% covPars, logisticPars,
                                                                 par1 = P(sim)$lowerSpreadProb))
    }

    coef <- ifelse(P(sim)$coefToUse == "bestCoef", "best coefficients", "averaged coefficients")
    fuelTypes <- setdiff(names(sim$fireSense_ignitionCovariates), c("pixelID")) #maybe others?

    sim$spreadProbFuelType <- plotSpreadProbByFuelType(spreadProbFuelType = sim$spreadProbFuelType,
                                                       typesOfFuel = fuelTypes,
                                                       coefToUse = coef,
                                                       covMinMax = sim$covMinMax)
  }

  # Return to raster format
  sim$fireSense_SpreadPredicted <- raster(sim$flammableRTM)
  sim$fireSense_SpreadPredicted[whNotNA] <- fireSenseDataDTx1000$spreadProb
  sim$spreadPredictedProbability[[paste0("Year", time(sim))]] <- sim$fireSense_SpreadPredicted

  invisible(sim)
}
