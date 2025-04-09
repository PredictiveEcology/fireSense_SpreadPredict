---
title: "fireSense_SpreadPredict Manual"
subtitle: "v.0.0.1"
date: "Last updated: 2025-04-08"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
bibliography: citations/references_fireSense_SpreadPredict.bib
link-citations: true
always_allow_html: true
---

# fireSense_SpreadPredict Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:fireSense-SpreadPredict) *fireSense_SpreadPredict*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Jean Marchal <jean.d.marchal@gmail.com> [aut], Eliot McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut, cre], Tati Michelleti <tati.micheletti@gmail.com> [aut], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Alex M. Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Predicts a surface (raster) of fire spread probabilities using a model previously fit with `fireSense_SpreadFit`.
fireSense [@Marchal:2017a; @Marchal:2017b; @Marchal:2019] <!-- TODO -->

### Module inputs and parameters

Describe input data required by the module and how to obtain it (e.g., directly from online sources or supplied by other modules)
If `sourceURL` is specified, `downloadData("fireSense_SpreadPredict", "..")` may be sufficient.

Table \@ref(tab:moduleInputs-fireSense-SpreadPredict) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-fireSense-SpreadPredict)(\#tab:moduleInputs-fireSense-SpreadPredict)List of (ref:fireSense-SpreadPredict) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> covMinMax_spread </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> range used to rescale coefficients during spreadFit </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_SpreadCovariates </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> data.table of covariates with pixelID column corresponding to flammableRTM index. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_SpreadFitted </td>
   <td style="text-align:left;"> fireSense_SpreadFit </td>
   <td style="text-align:left;"> An object of class 'fireSense_SpreadFit' created by the fireSense_SpreadFit module. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableRTM </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> RTM with nonflammable pixels coded as 0 and flammable as 1. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Summary of user-visible parameters (Table \@ref(tab:moduleParams-fireSense-SpreadPredict))


<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-fireSense-SpreadPredict)(\#tab:moduleParams-fireSense-SpreadPredict)List of (ref:fireSense-SpreadPredict) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> climCol </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> MDC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> the name of the climate covariate in `sim$fireSense_spreadCovariates` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> coefToUse </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> meanCoef </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Which coefficient to use to predict? The best coefficient (bestCoef) from DEOPtim or the average (meanCoef; default). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lowerSpreadProb </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Lower spread probability </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mutuallyExclusiveCols </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> vegPC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> a named list, where the name of the list must be a covariate in the data.table. Covariates matching the values in each list element will be set to 0. List content should be a grep regex. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> when to start this module? By default, the start time of the simulation. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between two runs of this module, expressed in units of simulation time.Defaults to 1 year. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. When to start saving output to a file. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between save events. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant </td>
  </tr>
</tbody>
</table>

### Events
<!-- TODO -->
- Module initiation;
- Make predictions.

### Plotting

Write what is plotted. <!-- TODO -->

### Saving

There is currently nothing saved, but this may change in the future. <!-- TODO -->

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-fireSense-SpreadPredict)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-fireSense-SpreadPredict)(\#tab:moduleOutputs-fireSense-SpreadPredict)List of (ref:fireSense-SpreadPredict) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireSense_SpreadPredicted </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A raster layer of spread probabilities </td>
  </tr>
</tbody>
</table>

# Links to other modules

Predictions made with this module can be used with the fire spread component of landscape fire models (e.g., `fireSense`). <!-- TODO -->

### Getting help

- <https://github.com/PredictiveEcology/fireSense_SpreadPredict/issues>

## References

<!-- autogenerated from bibligraphy -->
