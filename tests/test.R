library(SpaDES)

## RasterLayer
  # mySim <- simInit(
  #   times = list(start = 1, end = 1, timeunit = "year"),
  #   modules = list("fireSense_SpreadPredict"),
  #   paths = list(modulePath = " # replace with empty string instead"),
  #   inputs = data.frame(
  #     files = c("Z:/fireSense_SpreadFitted.RData", "Z:/beta.tif", "Z:/theta.tif"),
  #     functions = c("load", "raster", "raster"),
  #     package = c("base", "raster", "raster"),
  #     stringsAsFactors = FALSE)
  # )

## RasterStack
  mySim <- simInit(
    times = list(start = 1, end = 1, timeunit = "year"),
    modules = list("fireSense_SpreadPredict"),
    paths = list(modulePath = " # replace with empty string instead"),
    inputs = data.frame(
      files = c("Z:/fireSense_SpreadFitted.RData", "Z:/beta__STACK.tif", "Z:/theta__STACK.tif"),
      functions = c("load", "stack", "stack"),
      package = c("base", "raster", "raster"),
      objectName = c("fires", "beta", "theta"),
      stringsAsFactors = FALSE)
  )

spades(mySim)

x11(); Plot(mySim$beta, mySim$theta, mySim$fireSense_SpreadPredictProb)
