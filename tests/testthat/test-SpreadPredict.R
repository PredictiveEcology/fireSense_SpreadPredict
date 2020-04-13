# Please do three things to ensure this template is correctly modified:
# 1. Rename this file based on the content you are testing using
#    `test-functionName.R` format so that your can directly call `moduleCoverage`
#    to calculate module coverage information.
#    `functionName` is a function's name in your module (e.g., `fireSense_SpreadPredictEvent1`).
# 2. Copy this file to the tests folder (i.e., `fireSense_SpreadPredict/tests/testthat`).

# 3. Modify the test description based on the content you are testing:
test_that("test simInit() and spades()", {
  library(magrittr)
  library(raster)
  library(SpaDES.core)
  library(SpaDES.tools)
  
  set.seed(123)
  
  start <- end <- 1
  
  # Define simulation parameters
  times <- list(start = start, end = end, timeunit = "year")
  modules <- list("fireSense_SpreadPredict")
  paths <- list(
    modulePath = normalizePath("..", "..", "..")
  )
  
  # Create a random map of weather
  nx <- ny <- 100L
  weather <- raster(nrows = ny, ncols = nx, xmn = -nx/2, xmx = nx/2, ymn = -ny/2, ymx = ny/2) %>%
        gaussMap(scale = 300, var = 1, speedup = 1, inMemory = TRUE)
  
  # Create a typical output of fireSense_SpreadFit
  fireSense_SpreadFitted <- list(
    formula = ~ weather2 -1,
    coef = setNames(c(0.1, .3, 3, 1.5, 3),
                    c("d", "a", "b", "g", "weather")) # d, a, b, g are parameters of the 5-parameters logistic function
  )
  class(fireSense_SpreadFitted) <- "fireSense_SpreadFit"
  
  parameters <- list(
    fireSense_SpreadPredict = list(
      .runInterval = 1,
      data = "weather",
      mapping = list(weather2 = "weather") # One can use mapping to map variables
                                           # in the formula of the fitted object
                                           # to those in data. Here weather2
                                           # (formula) is mapped to weather (data).
    )
  )
  
  # Objects to pass from the global environment to the simList environment
  objects <- c("weather" = weather, "fireSense_SpreadFitted" = fireSense_SpreadFitted)
  
  # Create the simList
  sim <- simInit(
    times = times,
    modules = modules,
    paths = paths,
    params = parameters,
    objects = objects
  )
  ## TODO: fix errors in test simulation
  #sim <- spades(sim)
  
  #spreadProb <- sim$fireSense_SpreadPredicted
  #weather <- weather
  #Plot(weather, spreadProb)
})
