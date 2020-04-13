# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_SpreadPredict",
  description = "Predicts a surface of fire spread probilities using a model fitted with fireSense_SpreadFit.",
  keywords = c("fire spread", "fireSense", "predict"),
  authors = c(
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))
  ),
  childModules = character(),
  version = list(fireSense_SpreadPredict = "0.0.1", SpaDES.core = "0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SpreadPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = "modelObjName", class = "character",
                    default = "fireSense_SpreadFitted",
                    desc = "a character vector indicating the name of a model
                            object created with the fireSense_SpreadFit module."),
    defineParameter(name = "data", class = "character",
                    default = "dataFireSense_SpreadPredict",
                    desc = "a character vector indicating the names of objects
                            in the `simList` environment in which to look for
                            variables present in the model formula. `data`
                            objects can be RasterLayers, RasterStacks or RasterBricks."),
    defineParameter(name = "mapping", class = "character, list", default = NULL,
                    desc = "optional named vector or list of character strings
                            mapping one or more variables in the model formula
                            to those in data objects."),
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start
                            time of the simulation."),
    defineParameter(name = ".runInterval", class = "numeric", default = 1,
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time. By default, 1 year."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA, 
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInterval", class = "numeric", default = NA, 
                    desc = "optional. Interval between save events."),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
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
  moduleName <- current(sim)$moduleName
  
  switch(
    eventType,
    init = { 
      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run")
      
      if (!is.na(P(sim)$.saveInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
    }, 
    run = { 
      sim <- spreadPredictRun(sim)
      
      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "run")
    },
    save = { 
      sim <- spreadPredictSave(sim)
      
      if (!is.na(P(sim)$.saveInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, moduleName, "save", .last())  
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

spreadPredictRun <- function(sim) 
{
  moduleName <- current(sim)$moduleName
  
  if (!is(sim[[P(sim)$modelObjName]], "fireSense_SpreadFit"))
    stop(moduleName, "> '", P(sim)$modelObjName, "' should be of class 'fireSense_SpreadFit")
  
  ## Toolbox: set of functions used internally by spreadPredictRun
    spreadPredictRaster <- function(model, data, par) 
    {
      par[1L] + (par[2L] - par[1L]) / (1 + (model %>%
        model.matrix(data) %>%
        `%*%` (par[5:length(par)]) %>%
        drop) ^ (-par[3L])) ^ par[4L]
    }
  
  # Load inputs in the data container
  # list2env(as.list(envir(sim)), envir = mod)

  mod_env <- new.env()
    
  for(x in P(sim)$data) 
  {
    if (!is.null(sim[[x]])) 
    {
      if (is(sim[[x]], "RasterStack") || is(sim[[x]], "RasterBrick")) 
      {
        list2env(setNames(unstack(sim[[x]]), names(sim[[x]])), envir = mod_env)
      } 
      else if (is(sim[[x]], "RasterLayer")) 
      {
        mod_env[[x]] <- sim[[x]]
      } 
      else stop(moduleName, "> '", x, "' is not a RasterLayer, a RasterStack or a RasterBrick.")
    }
  }
  
  ## In case there is a response in the formula remove it
  terms <- sim[[P(sim)$modelObjName]]$formula %>% terms.formula %>% delete.response
  
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
  
  missing <- !allxy %in% ls(mod_env, all.names = TRUE)
  
  if (s <- sum(missing))
    stop(moduleName, "> '", allxy[missing][1L], "'",
         if (s > 1) paste0(" (and ", s-1L, " other", if (s>2) "s", ")"),
         " not found in data objects.")

  sim$fireSense_SpreadPredicted <- mget(allxy, envir = mod_env, inherits = FALSE) %>%
    stack %>%
    predict(model = formula, fun = spreadPredictRaster, na.rm = TRUE, par = sim[[P(sim)$model]]$coef)
  
  invisible(sim)
}


### template for save events
spreadPredictSave <- function(sim) 
{
  moduleName <- current(sim)$moduleName
  timeUnit <- timeunit(sim)
  currentTime <- time(sim, timeUnit)
  
  saveRDS(
    sim$fireSense_SpreadPredicted, 
    file = file.path(paths(sim)$out, paste0("fireSense_SpreadPredicted_", timeUnit, currentTime, ".rds"))
  )
  
  invisible(sim)
}
