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
