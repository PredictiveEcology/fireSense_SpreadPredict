# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_SpreadPredict",
  description = "Predict a probability surface of fire spread probilities from a model fitted using fireSense_SpreadFit.",
  keywords = c("fire spread", "fireSense", "predict"),
  authors = c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = "data", class = "character", default = NULL,
      desc = "optional. A character vector indicating the names of objects in the
              simList environment in which to look for variables in the model. 
              Data objects can be named lists of RasterLayers or RasterStacks
              (for time series), but should all be of one unique class, e.g. 
              RasterLayer. If omitted, or if variables are not found in data
              objects, variables are searched in the simList environment."),
    defineParameter(name = "mapping", class = "character, list", default = NULL,
      desc = "optional. Named character vector or list mapping some or all 
              variables in the model to those in data objects."),
    defineParameter(name = "initialRunTime", class = "numeric", default = start(sim),
      desc = "optional. Simulation time at which to start this module. Defaults 
              to simulation start time."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA,
      desc = "optional. Interval in simulation time units between two runs of
              this module.")
  ),
  inputObjects = data.frame(
    objectName = "fireSense_SpreadFitted",
    objectClass = "fireSense_SpreadFit",
    sourceURL = "",
    other = NA_character_,
    stringsAsFactors = FALSE
  ),
  outputObjects = data.frame(
    objectName = "fireSense_SpreadPredict",
    objectClass = "raster",
    other = NA_character_,
    stringsAsFactors = FALSE
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_SpreadPredict = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    sim <- sim$fireSense_SpreadPredictInit(sim)

  } else if (eventType == "run") {
    sim <- sim$fireSense_SpreadPredictRun(sim)

  } else if (eventType == "save") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event
    
    # e.g., call your custom functions/methods here
    # you can define your own methods below this `doEvent` function
    
    # schedule future event(s)
    
    # e.g.,
    # sim <- scheduleEvent(sim, time(sim) + increment, "fireSense_SpreadPredict", "save")
    
    # ! ----- STOP EDITING ----- ! #
    
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  invisible(sim)
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
fireSense_SpreadPredictInit <- function(sim) {

  sim <- scheduleEvent(sim, eventTime = p(sim)$initialRunTime, "fireSense_SpreadPredict", "run")
  sim

}

fireSense_SpreadPredictRun <- function(sim) {
  
  ## Toolbox: set of functions used internally by fireSense_SpreadPredictRun
    ## Raster predict function
      fireSense_SpreadPredictRaster <- function(model, data, par) {
        
        par[3L] + par[1L] / (1 + (model %>%
          model.matrix(data) %>%
          `%*%` (par[5:length(par)]) %>%
          drop) ^ (-par[2L])) ^ par[4L]
        
      }
  
  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))
  list2env(as.list(envir(sim)), envir = envData)

  if (!is.null(p(sim)$data)) ## Handling data arg
    lapply(p(sim)$data, function(x, envData) if (is.list(sim[[x]])) list2env(sim[[x]], envir = envData), envData = envData)

  ## In case there is a response in the formula remove it
  terms <- sim$fireSense_SpreadFitted$formula %>% terms.formula %>% delete.response
  
  ## Mapping variables names to data
  if (!is.null(p(sim)$mapping)) {
    
    for (i in 1:length(p(sim)$mapping)) {
      
      attr(terms, "term.labels") <- gsub(pattern = names(p(sim)$mapping[i]),
                                         replacement = p(sim)$mapping[[i]], x = attr(terms, "term.labels"))
      
    }
    
  }
  
  formula <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula)
  
  if (all(unlist(lapply(allxy, function(x) is(envData[[x]], "RasterStack"))))) {

    sim$fireSense_SpreadPredictProb <- mget(allxy, envir = envData, inherits = FALSE) %>%
      lapply(unstack) %>%
      c(list(FUN = function(...) stack(list(...)), SIMPLIFY = FALSE)) %>%
      do.call("mapply", args = .) %>%
      lapply(function(x) 
        predict(x, model = formula, fun = fireSense_SpreadPredictRaster, na.rm = TRUE, par = sim$fireSense_SpreadFitted$coef)) %>%
      stack

  } else if (all(unlist(lapply(allxy, function(x) is(envData[[x]], "RasterLayer"))))) {

    sim$fireSense_SpreadPredictProb <- mget(allxy, envir = envData, inherits = FALSE) %>%
      stack %>%
      predict(model = formula, fun = fireSense_SpreadPredictRaster, na.rm = TRUE, par = sim$fireSense_SpreadFitted$coef)
    
  } else {
  
    exist <- allxy %in% ls(envData)
    class <- unlist(lapply(allxy, function(x) is(envData[[x]], "RasterLayer") || is(envData[[x]], "RasterStack")))
    
    if (any(!exist)) {
      stop(paste0("fireSense_SpreadPredict> Variable '", allxy[which(!exist)[1L]], "' not found."))
    } else if (any(class)) {
      stop("fireSense_SpreadPredict> At least one of the data objects is not a RasterLayer.")
    } else {
      stop(paste0("fireSense_SpreadPredict> Variable '", allxy[which(!class)[1L]], "' does not match a RasterLayer."))
    }
  }

  if (!is.na(p(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + p(sim)$intervalRunModule, "fireSense_SpreadPredict", "run")
  
  sim

}


### template for save events
fireSense_SpreadPredictSave <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
