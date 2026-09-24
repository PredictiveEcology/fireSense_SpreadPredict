## Several fitted ELFs in one study area (the 2-ELF Mackenzie forecast, 2026-09). Every pixel gets a spread
## probability. Each ELF's model predicts its own pixels and those within ELFblendWidth of them; where two
## overlap they are averaged with weights falling linearly from 1 inside an ELF to 0 at ELFblendWidth
## outside it, normalised to sum to 1 (Eliot, 2026-09-22).
##
## Toy: one row of 10 pixels, 5 km wide; ELF "A" is columns 1-5, "B" columns 6-10. All coefficients are 0,
## so each ELF predicts a constant: lower + (maxAsymptote - lower) / 2 = 0.19 for A (max 0.25) and
## 0.20 for B (max 0.27). A is fitted on fuelA, B on fuelB; the covariate table holds both.

multiInputs <- function() {
  crs <- "EPSG:3005"
  flam <- terra::rast(nrows = 1, ncols = 10, xmin = 0, xmax = 50000, ymin = 0, ymax = 5000, crs = crs, vals = 1)
  elf <- terra::rast(flam); terra::values(elf) <- rep(1:2, each = 5)
  levels(elf) <- data.frame(id = 1:2, ELFind = c("A", "B"))
  covs <- data.table::data.table(pixelID = 1:10, MDC = 100, fuelA = 2, fuelB = 3)
  pars <- function(max, fuel) {
    d <- data.frame(maxAsymptote = max, hillSlope1 = 2, inflectionPoint1 = 1, MDC = 0, x = 0)
    names(d)[5] <- fuel
    d
  }
  cmm <- function(fuel) { d <- data.table::data.table(MDC = c(0, 200), f = c(0, 8)); data.table::setnames(d, "f", fuel); d }
  sa <- data.frame(polygonID = c("A", "B"))
  sa$params <- list(pars(0.25, "fuelA"), pars(0.27, "fuelB"))
  sa$covMinMax_spread <- list(cmm("fuelA"), cmm("fuelB"))
  list(flammableRTM = flam, rasterToMatchLargeELF = elf, fireSense_SpreadCovariates = covs,
       studyAreaWithSpreadParams = sa, fireSense_spreadFormula = "~ MDC + fuelA - 1",
       covMinMax_spread = cmm("fuelA"))
}

test_that("each ELF predicts its own pixels, and the boundary is a distance-weighted blend", {
  v <- predVals(toyRun(multiInputs(), params = list(ELFblendWidth = 20000)))
  ## raw weights at pixel centres (5 km apart): A = 1 on 1-5, then 0.75, 0.5, 0.25, 0, 0 on 6-10
  rA <- c(1, 1, 1, 1, 1, 0.75, 0.5, 0.25, 0, 0); rB <- rev(rA)
  expect_equal(v, (rA * 0.19 + rB * 0.20) / (rA + rB), tolerance = 1e-9)
  expect_equal(v[1], 0.19, tolerance = 1e-9)             # 25 km from B: A only
  expect_equal(v[10], 0.20, tolerance = 1e-9)
  expect_equal(v[5], (0.19 + 0.75 * 0.20) / 1.75, tolerance = 1e-9)
  expect_true(all(is.finite(v)))                         # every pixel has a probability
})

test_that("a narrow blend width leaves only the boundary pixels mixed", {
  v <- predVals(toyRun(multiInputs(), params = list(ELFblendWidth = 6000)))
  ## 5 km from the other ELF: raw weight 1 - 5/6
  expect_equal(v[c(1:4, 7:10)], rep(c(0.19, 0.20), each = 4), tolerance = 1e-9)
  expect_equal(v[5], (0.19 + (1/6) * 0.20) / (1 + 1/6), tolerance = 1e-9)
})

test_that("with several ELFs, a missing ELF raster stops with a message that says what is needed", {
  ins <- multiInputs(); ins$rasterToMatchLargeELF <- NULL
  expect_error(toyRun(ins), "rasterToMatchLargeELF")
})
