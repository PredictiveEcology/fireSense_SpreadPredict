# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_SpreadPredict",
  description = "Predicts a surface of fire spread probilities using a model fitted with fireSense_SpreadFit.",
  keywords = c("fire spread", "fireSense", "predict"),
  authors = c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
  childModules = character(),
  version = list(SpaDES.core = "0.1.0", fireSense_SpreadPredict = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = "modelName", class = "character",
                    default = "fireSense_SpreadFitted",
                    desc = "a character vector indicating the name of a model
                            object created with the fireSense_SpreadFit module."),
    defineParameter(name = "data", class = "character",
                    default = "dataFireSense_SpreadPredict",
                    desc = "a character vector indicating the names of objects
                            in the `simList` environment in which to look for
                            variables present in the model formula. `data`
                            objects can be RasterLayers or RasterStacks. If
                            variables are not found in `data` objects, they are
                            searched in the `simList` environment. For time 
                            series, `data` objects must be named lists of 
                            RasterStacks or RasterLayers, named starting with 
                            `start(simList)` and ending with `end(simList)` such
                            that variables can be matched for every 
                            `timeunit(simList)` of the simulation."),
    defineParameter(name = "mapping", class = "character, list", default = NULL,
                    desc = "optional named vector or list of character strings
                            mapping one or more variables in the model formula
                            to those in data objects."),
    defineParameter(name = "initialRunTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start
                            time of the simulation."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA,
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time."),
    defineParameter(".useCache", "numeric", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = rbind(
    expectsInput(
      objectName = "fireSense_SpreadFitted",
      objectClass = "fireSense_SpreadFit",
      sourceURL = NA_character_,
      desc = "An object of class 'fireSense_SpreadFit' created by the fireSense_SpreadFit module."
    ),
    expectsInput(
      objectName = "dataFireSense_SpreadPredict",
      objectClass = "RasterLayer, RasterStack",
      sourceURL = NA_character_,
      desc = "One or more RasterLayers or RasterStacks in which to look for variables present in the model formula."
    )
  ),
  outputObjects = createsOutput(
    objectName = "fireSense_SpreadPredicted",
    objectClass = "RasterLayer, RasterStack",
    desc = "An object whose class depends on that of the inputs, could be a RasterLayer or a RasterStack."
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_SpreadPredict = function(sim, eventTime, eventType, debug = FALSE) 
{
  switch(
    eventType,
    init = { sim <- sim$fireSense_SpreadPredictInit(sim) }, 
    run = { sim <- sim$fireSense_SpreadPredictRun(sim) },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      
      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function
      
      # schedule future event(s)
      
      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "fireSense_SpreadPredict", "save")
      
      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  
  invisible(sim)
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
fireSense_SpreadPredictInit <- function(sim)
{
  sim <- scheduleEvent(sim, eventTime = P(sim)$initialRunTime, current(sim)$moduleName, "run")
  invisible(sim)
}

fireSense_SpreadPredictRun <- function(sim) 
{
  stopifnot(is(sim[[P(sim)$modelName]], "fireSense_SpreadFit"))
  
  moduleName <- current(sim)$moduleName
  currentTime <- time(sim, timeunit(sim))
  endTime <- end(sim, timeunit(sim))
  
  ## Toolbox: set of functions used internally by fireSense_SpreadPredictRun
    ## Raster predict function
    fireSense_SpreadPredictRaster <- function(model, data, par) 
    {
      par[3L] + par[1L] / (1 + (model %>%
        model.matrix(data) %>%
        `%*%` (par[5:length(par)]) %>%
        drop) ^ (-par[2L])) ^ par[4L]
    }
    
  # Create a container to hold the data
  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))

  for(x in P(sim)$data) 
  {
    if (!is.null(sim[[x]][[as.character(currentTime)]])) 
    {
      if (is(sim[[x]][[as.character(currentTime)]], "RasterStack")) 
      {
        list2env(setNames(unstack(sim[[x]][[as.character(currentTime)]]), names(sim[[x]][[as.character(currentTime)]])), envir = envData)
      } 
      else if (is(sim[[x]][[as.character(currentTime)]], "RasterLayer")) 
      {
        envData[[x]] <- sim[[x]][[as.character(currentTime)]]
      } 
      else stop(paste0(moduleName, "> '", x, "' is not a RasterLayer or a RasterStack."))
    }
  }
  
  ## In case there is a response in the formula remove it
  terms <- sim[[P(sim)$modelName]]$formula %>% terms.formula %>% delete.response
  
  ## Mapping variables names to data
  if (!is.null(P(sim)$mapping))
  {
    for (i in 1:length(P(sim)$mapping)) 
    {
      attr(terms, "term.labels") %<>% gsub(
        pattern = names(P(sim)$mapping[i]),
        replacement = P(sim)$mapping[[i]],
        x = .
      )
    }
  }

  formula <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula)
  
  missing <- !allxy %in% ls(envData, all.names = TRUE)
  
  if (s <- sum(missing))
    stop(paste0(moduleName, "> '", allxy[missing][1L], "'", if (s > 1) paste0(" (and ", s-1L, " other", if (s>2) "s", ")"),
                " not found in data objects."))

  sim$fireSense_SpreadPredicted[as.character(currentTime)] <- list(
    mget(allxy, envir = envData, inherits = FALSE) %>%
      stack %>%
      predict(model = formula, fun = fireSense_SpreadPredictRaster, na.rm = TRUE, par = sim[[P(sim)$model]]$coef)
  )    
  
  if (!is.na(P(sim)$intervalRunModule) && (currentTime + P(sim)$intervalRunModule) <= endTime) # Assumes time only moves forward
    sim <- scheduleEvent(sim, currentTime + P(sim)$intervalRunModule, moduleName, "run")
  
  invisible(sim)
}


### template for save events
fireSense_SpreadPredictSave <- function(sim) 
{
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
