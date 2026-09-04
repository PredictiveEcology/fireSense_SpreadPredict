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
  version = list(fireSense_SpreadPredict = "0.0.2", SpaDES.core = "0.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "Matrix", "methods", "terra", "SpaDES.core (>=3.0.4)", "stats",
                  "ggplot2", "viridis",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.1.0)"),
  parameters = bindrows(
    # defineParameter(name = "climCol", class = "character", default = "MDC", min = NA, max = NA,
    #                 desc = "the name of the climate covariate in `sim$fireSense_spreadCovariates`"),
    defineParameter(name = "coefToUse", class = "character", default = "meanCoef",
                    desc = paste("Which coefficient to use to predict?",
                                 "The best coefficient (bestCoef) from DEOPtim or ",
                                 "the average (meanCoef; default).")),
    defineParameter(name = "lowerSpreadProb", class = "numeric", default = 0.13,
                    desc = "Lower spread probability"),
    defineParameter("maxFireSpread", "numeric", default = 0.28,
                    desc = paste0("optional. Maximum fire spread average to be passed to the `.objFun`.",
                                  "This puts an upper limit on `spreadProb` during optimization. Expected to be same ",
                                  "as fireSpread_SpreadFit parameter value")),
    defineParameter(name = "mutuallyExclusiveCols", "list", default = list("youngAge" = "fuels"), NA, NA,
                    desc = paste("a named list, where the name of the list must be a covariate in the data.table.",
                                 "Covariates matching the values in each list element will be set to 0.",
                                 "List content should be a grep regex.")),
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start time of the simulation."),
    defineParameter(name = ".runInterval", class = "numeric", default = 1,
                    desc = paste("optional. Interval between two runs of this module, expressed in units of",
                                 "simulation time.Defaults to 1 year.")),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA,
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInterval", class = "numeric", default = NA,
                    desc = "optional. Interval between save events."),
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

  # if (!is(sim$fireSense_SpreadFitted, "fireSense_SpreadFit")) {
  #   stop(moduleName, "> '", sim$fireSense_spreadFitted, "' should be of class 'fireSense_SpreadFit")
  # }

  # Load inputs in the data container
  mod_env <- new.env(parent = globalenv())
  list2env(fireSense_SpreadCovariates, envir = mod_env)
  ## In case there is a response in the formula remove it

  if (FALSE) {
    # The old way prior to June 2, 2025

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
    # convert to sim$flammableRTMs
    sim$fireSense_SpreadPredicted <- rast(sim$flammableRTM) ## use flammableRTM as template
    ## Need to track what is happening with missing pixels
    if (FALSE) {
      nFlam <- sum(getValues(sim$flammableRTM), na.rm = TRUE)
      nLand <- nrow(sim$landcoverDT)
      nSpread <- nrow(sim$fireSense_SpreadCovariates)
      nIg <- nrow(sim$fireSense_IgnitionAndEscapeCovariates)
    }
    sim$fireSense_SpreadPredicted[fireSense_SpreadCovariates$pixelID] <- fireSense_SpreadCovariates$spreadProb

  } else {
    # sim$studyAreaWithSpreadParams


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

    ###################################################
    # Convert stacks to lists of data.table objects --> much more compact
    ###################################################
    # First for stacks that are "annual"
    shortAnnDTx1000 <- toX1000(list(fireSense_SpreadCovariates))[[1]] |> setDT()
    colsToUse <- setdiff(names(fireSense_SpreadCovariates), "pixelID")

    logisticPars <- sim$studyAreaWithSpreadParams$params[[1]]
    
    shortAnnDT <- # shortAnnDTx1000 <-
      spreadProbFromIntegerCovs(shortAnnDTx1000 = shortAnnDTx1000, # annDTx1000, nonAnnualDTx1000,
                                # indexNonAnnual, 
                                yr = time(sim),
                                covMinMax = sim$covMinMax_spread,
                                mutuallyExclusive = NULL, # alraedy done in dataPrepPredict
                                colsToUse = colsToUse,
                                doAssertions = FALSE,
                                logisticPars = logisticPars, # covPars, 
                                maxFireSpread = Par$maxFireSpread
                                # , lowerSpreadProb # this is set at 0.13
                                )

    # THIS IS INSIDE fireSenseUtils::objFunInner
    # set(shortAnnDTx1000, NULL, "spreadProb",
    #     logisticAll(logisticPars, mat = as.matrix(shortAnnDTx1000[, ..colsToUse]), covPars, lowerSpreadProb))
    
    parsModel <- length(colsToUse)
    # par <- purrr::pmap(.l = list(ind = seq(NROW(sim$studyAreaWithSpreadParams$params[[1]]))),
    #                    sa = sim$studyAreaWithSpreadParams, function(ind, sa) {
    #                      sa$params[[1]][ind,] |> as.vector() |> unlist()
    #                    })
    mat <- as.matrix(shortAnnDT[, ..colsToUse])

    # for replicate "best" params from DEoptim
    spreadProbList <- purrr::pmap(.l = list(ind = seq(NROW(sim$studyAreaWithSpreadParams$params[[1]]))),
                       sa = sim$studyAreaWithSpreadParams, function(ind, sa) {
                         par <- sa$params[[1]][ind,] |> as.vector() |> unlist()
                         # mat <- as.matrix(fireSense_SpreadCovariates[, ..colsToUse])/1000 # Divide by 1000 for the model prediction
                         
                         covPars <- intersect(names(par), colsToUse)
                         covPars <- par[covPars]
                         logisticPars <- par[setdiff(names(par), names(covPars))]
                         # params <- paramsSeparate(par, parsModel)
                         # matrix multiplication
                         # covPars <- params$covPars
                         # logisticPars <- params$logisticPars
                         # Make sure the order is correct in the matrix
                         matching <- intersect(names(covPars), colnames(mat))
                         missingCovs <- setdiff(colnames(mat), names(covPars))
                         if (length(missingCovs))
                           warning("There are covariates in the sim$fireSense_SpreadCovariates: \n",
                                paste0(missingCovs, collapse = ", "),
                                "\n...that are not in the sim$studyAreaWithSpreadParams")
                         mat <- mat[, matching]

                         logisticAll(logisticPars, #fireSense_SpreadCovariates,
                                     mat, covPars, P(sim)$lowerSpreadProb)
                       })
    spreadProbMat <- do.call(cbind, spreadProbList)
    
    # spreadProbMat <- spreadProbMat[, which.min(colMeans(spreadProbMat))]

    set(shortAnnDT, NULL, "spreadProb", rowMeans(spreadProbMat))

    # par <- sim$studyAreaWithSpreadParams$params[[1]][1,] |> as.vector() |> unlist()
    # mat <- as.matrix(fireSense_SpreadCovariates[, ..colsToUse])/1000 # Divide by 1000 for the model prediction
    #
    # # matrix multiplication
    # covPars <- tail(x = par, n = parsModel)
    # logisticPars <- head(x = par, n = length(par) - parsModel)
    # # Make sure the order is correct in the matrix
    # matching <- match(names(covPars), colnames(mat))
    # mat <- mat[, matching]
    #
    # preds <- logisticAll(logisticPars, mat, covPars, P(sim)$lowerSpreadProb)
    # if (length(logisticPars) == 4) {
    #   set(fireSense_SpreadCovariates, NULL, "spreadProb", logistic4p(mat %*% covPars, logisticPars))
    # } else if (length(logisticPars) == 3) {
    #   set(fireSense_SpreadCovariates, NULL, "spreadProb", logistic3p(mat %*% covPars, logisticPars,
    #                                                                  par1 = P(sim)$lowerSpreadProb))
    # } else if (length(logisticPars) == 2) {
    #   set(fireSense_SpreadCovariates, NULL, "spreadProb", logistic2p(mat %*% covPars, logisticPars,
    #                                                                  par1 = P(sim)$lowerSpreadProb))
    # }

    # Return to raster format
    sim$fireSense_SpreadPredicted <- rast(sim$flammableRTM) ## use flammableRTM as template
    ## Need to track what is happening with missing pixels
    if (FALSE) {
      nFlam <- sum(getValues(sim$flammableRTM), na.rm = TRUE)
      nLand <- nrow(sim$landcoverDT)
      nSpread <- nrow(sim$fireSense_SpreadCovariates)
      nIg <- nrow(sim$fireSense_IgnitionAndEscapeCovariates)
    }
    sim$fireSense_SpreadPredicted[shortAnnDT$pixelID] <- shortAnnDT$spreadProb

  }

  invisible(sim)
}



