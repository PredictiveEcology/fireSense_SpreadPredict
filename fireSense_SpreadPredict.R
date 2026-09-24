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
  version = list(fireSense_SpreadPredict = "1.0.0.9005", SpaDES.core = "0.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "Matrix", "methods", "terra", "SpaDES.core (>=3.0.4)", "stats",
                  "ggplot2", "viridis",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.2.3.9038)"),
  parameters = bindrows(
    defineParameter(name = "lowerSpreadProb", class = "numeric", default = 0.13,
                    desc = "Lower asymptote of the 2- and 3-parameter logistic."),
    defineParameter("ELFblendWidth", "numeric", default = 20000,
                    desc = paste("With several ELFs: each ELF's model also predicts this far (m) outside its own",
                                 "pixels, and where predictions overlap they are averaged with weights that fall",
                                 "linearly from 1 inside the ELF to 0 at this distance outside it. 50/50 at a",
                                 "boundary. The default is the buffer fireSenseUtils::makeELFs() puts around ELFs.")),
    defineParameter("maxFireSpread", "numeric", default = 0.28,
                    desc = paste("Upper limit on `spreadProb` used when fitting. Here it is only checked",
                                 "to be the same in every module that defines it.")),
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "Time of the first prediction."),
    defineParameter(name = ".runInterval", class = "numeric", default = 1,
                    desc = "Interval between predictions, in years. `NA` predicts once."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA,
                    desc = "Time of the `save` event, which does nothing. `NA` means never."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used."),
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
    expectsInput(objectName = "rasterToMatchLargeELF", objectClass = "SpatRaster", sourceURL = NA,
                 desc = paste("Only with several fitted ELFs: each pixel's ELF (`ELFind`), on the grid of",
                              "`flammableRTM`, from `fireSense_ELFs` with a `studyAreaLarge`.")),
    expectsInput(objectName = "flammableRTM", objectClass = "SpatRaster", sourceURL = NA,
                 desc = "Binary raster, 1 where the pixel is flammable. Template for `fireSense_SpreadPredicted`.")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "fireSense_SpreadPredicted", objectClass = "SpatRaster",
                  desc = "Spread probability of each flammable pixel, this year."),
    createsOutput(objectName = "fireSense_SpreadSD", objectClass = "SpatRaster|numeric",
                  desc = paste("The fitted sd of the per-year random effect on logit spread probability",
                               "(`yearSpreadSD`; 0 if the fit has none), for `fireSense`. One number with one",
                               "fitted ELF; with several, a raster blended across ELFs with the weights of",
                               "`fireSense_SpreadPredicted`."))
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
#' Rescales `sim$fireSense_SpreadCovariates` with the fit's covariate ranges, computes the spread
#' probability for each parameter set (row) of the fit, and writes the mean over parameter sets to
#' `sim$fireSense_SpreadPredicted`.
#'
#' With one fitted ELF (one row of `sim$studyAreaWithSpreadParams`) every pixel uses its parameters and
#' `sim$covMinMax_spread`. With several, each ELF's model (its parameter sets, its `covMinMax_spread`, and
#' only the covariates it was fitted with) predicts its own pixels, from `sim$rasterToMatchLargeELF`, and
#' those within `ELFblendWidth` of them; overlapping predictions are averaged with the weights of
#' `ELFblendWeights()`.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly.
spreadPredictRun <- function(sim) {
  covs <- copy(sim$fireSense_SpreadCovariates)
  sa <- sim$studyAreaWithSpreadParams

  # Without fitted parameters there is nothing to predict from; say so instead of
  # dying in rowMeans() on an empty matrix (which is what an unfitted ELF produced
  # when fireSense_SpreadFit had not run first). This must come before anything
  # indexes `params[[1]]`: with zero rows that fails first, "subscript out of bounds".
  nPar <- tryCatch(NROW(sa$params[[1]]), error = function(e) 0L)
  if (NROW(sa) == 0L || is.null(nPar) || nPar == 0L)
    stop("fireSense_SpreadPredict: sim$studyAreaWithSpreadParams holds no fitted spread ",
         "parameters for this run (", if (!is.null(sim$.runName)) sim$.runName else "unknown",
         "). Either fireSense_SpreadFit has not run yet -- its `run` event must precede this ",
         "module's -- or the shared ledger has no row for this polygon.", call. = FALSE)

  if (NROW(sa) == 1L) {
    pred <- spreadProbOneELF(covs, params = sa$params[[1]], covMinMax = sim$covMinMax_spread,
                             formula = sim$fireSense_spreadFormula, yr = time(sim),
                             maxFireSpread = P(sim)$maxFireSpread, lowerSpreadProb = P(sim)$lowerSpreadProb)
    sim$fireSense_SpreadSD <- yearSpreadSDOf(sa$params[[1]])
  } else {
    ids <- as.character(sa[[fireSenseUtils::polygonIDTxt]])
    ## each ELF's weight at each pixel: static, so computed once
    if (is.null(mod$ELFweights) || !identical(attr(mod$ELFweights, "key"), list(ids, covs$pixelID)))
      mod$ELFweights <- ELFblendWeights(sim$rasterToMatchLargeELF, sim$flammableRTM, covs$pixelID, ids,
                                        width = P(sim)$ELFblendWidth)
    w <- mod$ELFweights
    none <- rowSums(w) == 0
    if (any(none))
      warning("fireSense_SpreadPredict: ", sum(none), " flammable pixels are more than ",
              P(sim)$ELFblendWidth, " m from every fitted ELF; they get no spread probability", call. = FALSE)
    acc <- numeric(NROW(covs)); wsum <- numeric(NROW(covs)); accSD <- numeric(NROW(covs))
    for (i in seq_along(ids)) {
      these <- which(w[, i] > 0)
      if (!length(these)) next
      p <- spreadProbOneELF(covs[these], params = sa$params[[i]], covMinMax = sa$covMinMax_spread[[i]],
                            formula = NULL, yr = time(sim), maxFireSpread = P(sim)$maxFireSpread,
                            lowerSpreadProb = P(sim)$lowerSpreadProb, byParams = TRUE)
      rows <- these[match(p$pixelID, covs$pixelID[these])]
      acc[rows] <- acc[rows] + w[rows, i] * p$spreadProb
      accSD[rows] <- accSD[rows] + w[rows, i] * yearSpreadSDOf(sa$params[[i]])
      wsum[rows] <- wsum[rows] + w[rows, i]
    }
    ok <- wsum > 0
    pred <- data.table(pixelID = covs$pixelID[ok], spreadProb = acc[ok] / wsum[ok])
    ## each ELF's year effect sd, blended like the probabilities: fireSense scales one z per year by it
    sim$fireSense_SpreadSD <- rast(sim$flammableRTM)
    sim$fireSense_SpreadSD[covs$pixelID[ok]] <- accSD[ok] / wsum[ok]
  }

  # Return to raster format
  sim$fireSense_SpreadPredicted <- rast(sim$flammableRTM) ## use flammableRTM as template
  sim$fireSense_SpreadPredicted[pred$pixelID] <- pred$spreadProb

  invisible(sim)
}

#' Each ELF's weight at each pixel
#'
#' An ELF's raw weight is 1 on its own pixels and falls linearly to 0 at `width` metres from them (distance
#' to the nearest pixel of the ELF). Weights are then normalised to sum to 1 at each pixel, so a pixel on a
#' boundary between two ELFs is 50/50 and one `width` inside an ELF uses that ELF only.
#'
#' @param elfRas `sim$rasterToMatchLargeELF`: `ELFind` per pixel, as a categorical raster (or the id).
#' @param template `sim$flammableRTM`; `pixelID` indexes its cells, so `elfRas` must share its grid.
#' @param pixelID integer; the cells to weight.
#' @param ids character; the ELFs, in the order of the rows of `studyAreaWithSpreadParams`.
#' @param width numeric; metres.
#' @return matrix, one row per `pixelID`, one column per ELF (`ids`), rows summing to 1 (or 0 where no
#'   ELF is within `width`). `attr(, "key")` records `ids` and `pixelID`.
ELFblendWeights <- function(elfRas, template, pixelID, ids, width) {
  if (is.null(elfRas))
    stop("fireSense_SpreadPredict: sim$studyAreaWithSpreadParams has several ELFs, so each pixel's ELF ",
         "must come from sim$rasterToMatchLargeELF (fireSense_ELFs, with a studyAreaLarge); it is missing",
         call. = FALSE)
  if (!isTRUE(terra::compareGeom(elfRas, template, stopOnError = FALSE)))
    stop("fireSense_SpreadPredict: sim$rasterToMatchLargeELF is not on the grid of sim$flammableRTM",
         call. = FALSE)
  r <- elfRas[[1]]
  v <- terra::values(r, mat = FALSE)
  lv <- terra::levels(r)[[1]]
  lab <- if (is.data.frame(lv) && NCOL(lv) >= 2) as.character(lv[[2]][match(v, lv[[1]])]) else as.character(v)
  w <- vapply(ids, function(id) {
    core <- r; terra::values(core) <- ifelse(lab %in% id, 1, NA)
    if (all(is.na(terra::values(core, mat = FALSE)))) return(numeric(length(pixelID)))
    d <- terra::values(terra::distance(core), mat = FALSE)[pixelID]
    pmax(0, 1 - d / width)
  }, numeric(length(pixelID)))
  w <- matrix(w, ncol = length(ids), dimnames = list(NULL, ids))
  tot <- rowSums(w)
  w[tot > 0, ] <- w[tot > 0, , drop = FALSE] / tot[tot > 0]
  attr(w, "key") <- list(ids, pixelID)
  w
}

#' Spread probability for the pixels of one ELF
#'
#' @param covs `data.table` of this ELF's pixels: `pixelID` and covariates (as `fireSense_dataPrepPredict`
#'   makes them; fuel biomass logged).
#' @param params `data.frame` of the ELF's fitted parameter sets, one per row (ledger `params`).
#' @param covMinMax the ELF's `covMinMax_spread`.
#' @param formula the spread formula, to check the covariates are all there.
#' @param byParams logical; `TRUE` (several ELFs) ignores `formula`: the covariates are those named in
#'   `params`, and the table is cut to them, since it holds every ELF's covariates.
#' @param yr,maxFireSpread,lowerSpreadProb as for `fireSenseUtils::spreadProbFromIntegerCovs()` and
#'   `fireSenseUtils::logisticAll()`.
#' @return `data.table` with `pixelID` and `spreadProb`, the mean over parameter sets.
spreadProbOneELF <- function(covs, params, covMinMax, formula, yr, maxFireSpread, lowerSpreadProb,
                             byParams = FALSE) {
  moduleName <- "fireSense_SpreadPredict"
  covs <- copy(covs)
  ## the per-year random effect is not a covariate coefficient: fireSense applies it (fireSense_SpreadSD)
  if (yearSpreadSDTxt %in% names(params)) {
    params <- as.data.frame(params)
    params[[yearSpreadSDTxt]] <- NULL
  }

  ## Fuel biomass arrives logged (fireSenseUtils::logMinB()). A fit made on LINEAR fuel biomass has
  ## fireSenseUtils::fuelLinearRange, c(0, 1e4), as that covariate's covMinMax_spread, and its
  ## coefficients only mean anything for biomass / 1e4: undo the log with the function the fit used.
  ## A fit made on the log scale has the log range there, and its covariates are left as they are.
  for (cn in intersect(names(covMinMax), names(covs))) {
    if (fireSenseUtils::isLinearFuelRange(covMinMax[[cn]]))
      covs[[cn]] <- fireSenseUtils::fuelLogToLinear(covs[[cn]])
  }

  ## the covariates this ELF was fitted with
  needed <- if (isTRUE(byParams)) {
    setdiff(names(params), unlist(fireSenseUtils::logisticParamNames))
  } else {
    terms <- delete.response(terms.formula(as.formula(formula)))
    all.vars(reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept")))
  }
  missing <- !needed %in% names(covs)
  if (s <- sum(missing)) {
    stop(
      moduleName, "> '", needed[missing][1L], "'",
      if (s > 1) paste0(" (and ", s - 1L, " other", if (s > 2) "s", ")"),
      " not found in data objects."
    )
  }
  ## with several ELFs the table holds every ELF's covariates; use this ELF's only
  if (isTRUE(byParams)) covs <- covs[, c("pixelID", needed), with = FALSE]

  # integers x 1000, the form `spreadProbFromIntegerCovs` expects
  shortAnnDTx1000 <- toX1000(list(covs))[[1]] |> setDT()
  colsToUse <- setdiff(names(covs), "pixelID")

  shortAnnDT <-
    spreadProbFromIntegerCovs(shortAnnDTx1000 = shortAnnDTx1000,
                              yr = yr,
                              covMinMax = covMinMax,
                              mutuallyExclusive = NULL, # alraedy done in dataPrepPredict
                              colsToUse = colsToUse,
                              doAssertions = FALSE,
                              logisticPars = params,
                              maxFireSpread = maxFireSpread
                              )

  mat <- as.matrix(shortAnnDT[, ..colsToUse])

  # for replicate "best" params from DEoptim
  spreadProbList <- lapply(seq_len(NROW(params)), function(ind) {
    par <- params[ind, ] |> as.vector() |> unlist()
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
    logisticAll(logisticPars, mat[, matching, drop = FALSE], covPars, lowerSpreadProb)
  })
  spreadProbMat <- do.call(cbind, spreadProbList)

  data.table(pixelID = shortAnnDT$pixelID, spreadProb = rowMeans(spreadProbMat))
}

yearSpreadSDTxt <- "yearSpreadSD"

#' The fitted sd of the per-year random effect
#'
#' `yearSpreadSD` (fireSense_SpreadFit, fireSenseUtils >= 0.2.3.9041) is one eps per year on logit spread
#' probability. With several retained parameter sets, their mean, as the spread probabilities are averaged.
#'
#' @param params `data.frame` of fitted parameters, one row per retained set.
#' @return Numeric; 0 when the fit has no `yearSpreadSD`.
yearSpreadSDOf <- function(params) {
  if (!yearSpreadSDTxt %in% names(params)) return(0)
  mean(params[[yearSpreadSDTxt]])
}
