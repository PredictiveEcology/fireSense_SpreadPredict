library(magrittr)
library(raster)
library(SpaDES)

set.seed(123)

# Define simulation parameters
times <- list(start = 1, end = 1, timeunit = "year")
modules <- list("fireSense_SpreadPredict")
paths <- list(
  modulePath = file.path("~/Documents/GitHub/McIntire-lab/modulesPrivate/")
)

# Create a random map of weather
nx <- ny <- 100L
weather <- raster(nrows = ny, ncols = nx, xmn = -nx/2, xmx = nx/2, ymn = -ny/2, ymx = ny/2) %>%
  gaussMap(scale = 300, var = 1, speedup = 1, inMemory = TRUE)

# Create a typical output of fireSense_SpreadFit
fireSense_SpreadFitted <- list(
  formula = ~ weather2 -1,
  coef = setNames(c(0.3, 3, 0.1, 1.5, 3),
                  c("A", "B", "D", "G", "weather")) # A, B, D, G are parameters of the 5-parameters logistic function
)
class(fireSense_SpreadFitted) <- "fireSense_SpreadFit"

parameters <- list(
  fireSense_SpreadPredict = list(
    mapping = list(weather2 = "weather") # One can use mapping to map variables
                                         # in the formula of the fitted object
                                         # to those in data. Here weather2
                                         # (formula) is mapped to weather (data).
  )
)

# Objects to pass from the global environment to the simList environment
objects <- c("weather", "fireSense_SpreadFitted")

# Create the simList
sim <- simInit(
  times = times,
  modules = modules,
  paths = paths,
  params = parameters,
  objects = objects
)

sim <- spades(sim)

x11(); Plot(weather, sim$fireSense_SpreadPredicted)
