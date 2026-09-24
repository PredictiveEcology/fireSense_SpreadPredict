# fireSense_SpreadPredict (development version)

- The fitted per-year random effect (`yearSpreadSD`, fireSense_SpreadFit / fireSenseUtils >= 0.2.3.9041) is no
  longer read as a covariate coefficient: with it, a single ELF treated it as a fourth logistic parameter and several
  ELFs stopped with "'yearSpreadSD' not found". It becomes the new output `fireSense_SpreadSD`, which fireSense
  scales one draw per year by: one number with one ELF (the mean over retained parameter sets), a raster blended
  with the spread probabilities' weights with several. 0 when the fit has none.

- Several fitted ELFs in one study area. Every pixel gets a spread probability: each ELF's model (its parameter sets, its `covMinMax_spread` and only the covariates it was fitted with, from its ledger row) predicts its own pixels and those within `ELFblendWidth` (default 20 km) of them, and overlapping predictions are averaged with weights falling linearly from 1 inside an ELF to 0 at `ELFblendWidth` outside it. Each pixel's ELF comes from the new input `rasterToMatchLargeELF` (fireSense_ELFs with a `studyAreaLarge`). One ELF works as before.
- A fit made with fireSense_SpreadFit's `link = "logistic3pUpper"` stores `upperTail1`; prediction uses the upper-tail link for it, chosen by the parameter's name (`fireSenseUtils::logisticAll()`). Needs fireSenseUtils >= 0.2.3.9038.

# fireSense_SpreadPredict 1.0.0

First release from `development` since `master` was last updated (2021-01-27). Full history: https://github.com/PredictiveEcology/fireSense_SpreadPredict/compare/a5b41f9...v1.0.0

## Breaking changes

- Removed input `dataFireSense_SpreadPredict` (RasterLayer, RasterStack).
- Removed output `spreadPredictedProbability` (list).
- Output `fireSense_SpreadPredicted` is now `SpatRaster` (was `RasterLayer, RasterStack`).
- Removed parameters: `data`, `mapping`, `modelObjName`, `typesOfFuel`.

## New features

- New inputs: `covMinMax_spread`, `fireSense_SpreadCovariates`, `flammableRTM`.
- New parameters: `maxFireSpread`, `mutuallyExclusiveCols`.

## Dependencies

- No longer depends on `raster`.
- Now depends on `terra`.

## Testing

- testthat suite and CI (`testthat-module`), including a snapshot of the module's inputs, outputs and parameters in `tests/testthat/test-metadata.R`.
