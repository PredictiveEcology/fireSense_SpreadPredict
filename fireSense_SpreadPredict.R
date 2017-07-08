# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_SpreadPredict",
  description = "Predict a probability surface of fire spread probilities from a model fitted using fireSense_SpreadFit.",
  keywords = c("fire spread", "fireSense", "predict"),
  authors = c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = "model", class = "character", default = "fireSense_SpreadFitted",
                    desc = "a character vector indicating the name of a model
                            object created with the fireSense_SpreadFit module."),
    defineParameter(name = "data", class = "character", default = "dataFireSense_SpreadPredict",
                    desc = "a character vector indicating the names of objects
                            in the `simList` environment in which to look for
                            variables in the model. `data` objects can be
                            RasterLayers or RasterStacks (for time series), or
                            named list of either several RasterLayers or several
                            RasterStacks, but RasterLayers and RasterStack can not
                            be mixed together. If omitted, or if variables are
                            not found in `data` objects, variables are searched
                            in the `simList` environment."),
    defineParameter(name = "mapping", class = "character, list", default = NULL,
                    desc = "optional named vector or list of character strings
                            mapping one or more variables in the model formula
                            to those in data objects."),
    defineParameter(name = "initialRunTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start
                            time of the simulation."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA,
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time.")
  ),
  inputObjects = rbind(
    expectsInput(
      objectName = "fireSense_SpreadFitted",
      objectClass = "fireSense_SpreadFit",
      sourceURL = NA_character_,
      desc = "An object of class fireSense_SpreadFitted created by the fireSense_SpreadFit module."
    ),
    expectsInput(
      objectName = "dataFireSense_SpreadPredict",
      objectClass = "RasterLayer, RasterStack, list",
      sourceURL = NA_character_,
      desc = "One or more RasterLayers or RasterStacks, or named lists of RasterLayers or RasterStacks in which to look for variables with which to predict."
    )
  ),
  outputObjects = createsOutput(
    objectName = "fireSense_SpreadPredicted",
    objectClass = "raster",
    desc = "An object whose class depends on that of the inputs, could be a RasterLayer or a RasterStack."
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

  sim <- scheduleEvent(sim, eventTime = P(sim)$initialRunTime, current(sim)$moduleName, "run")
  sim

}

fireSense_SpreadPredictRun <- function(sim) {

  moduleName <- current(sim)$moduleName
  
  ## Toolbox: set of functions used internally by fireSense_SpreadPredictRun
    ## Raster predict function
    fireSense_SpreadPredictRaster <- function(model, data, par) {
      
      par[3L] + par[1L] / (1 + (model %>%
        model.matrix(data) %>%
        `%*%` (par[5:length(par)]) %>%
        drop) ^ (-par[2L])) ^ par[4L]
      
    }
    
  # Check if the model object exists in the simList environment
  if (is.null(sim[[P(sim)$model]])) stop(paste0(moduleName, "> Model object '", P(sim)$model,"' does not exists in the simList environment."))
  
  # Create a container to hold the data
  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))
  
  # Load data in the container
  list2env(as.list(envir(sim)), envir = envData)
  
  lapply(P(sim)$data, function(x, envData) {
    
    if (!is.null(sim[[x]])) {
      
      if (is.list(sim[[x]]) && !is.null(names(sim[[x]]))) {
        
        list2env(sim[[x]], envir = envData)  
          
      } else stop(paste0(moduleName, "> '", x, "' is not a named list."))
      
    }
    
  }, envData = envData)

  ## In case there is a response in the formula remove it
  terms <- sim[[P(sim)$model]]$formula %>% terms.formula %>% delete.response
  
  ## Mapping variables names to data
  if (!is.null(P(sim)$mapping)) {
    
    for (i in 1:length(P(sim)$mapping)) {
      
      attr(terms, "term.labels") <- gsub(
        pattern = names(P(sim)$mapping[i]),
        replacement = P(sim)$mapping[[i]],
        x = attr(terms, "term.labels")
      )
      
    }
    
  }
  
  formula <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula)
  
  if (all(unlist(lapply(allxy, function(x) is(envData[[x]], "RasterLayer"))))) {

    sim$fireSense_SpreadPredicted <- mget(allxy, envir = envData, inherits = FALSE) %>%
      stack %>%
      predict(model = formula, fun = fireSense_SpreadPredictRaster, na.rm = TRUE, par = sim[[P(sim)$model]]$coef)
    
  } else {
  
    exist <- allxy %in% ls(envData)
    class <- unlist(lapply(allxy, function(x) is(envData[[x]], "RasterLayer")))
    
    if (any(!exist)) {
      stop(paste0(moduleName, "> Variable '", allxy[which(!exist)[1L]], "' not found."))
    } else if (any(class)) {
      stop(paste0(moduleName, "> At least one of the data objects is not a RasterLayer."))
    } else {
      stop(paste0(moduleName, "> Variable '", allxy[which(!class)[1L]], "' does not match a RasterLayer."))
    }
  }

  if (!is.na(P(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + P(sim)$intervalRunModule, moduleName, "run")
  
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
