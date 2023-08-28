# amazon_droughts

This repository provides the R scripts used in the analyses for the paper 'Critical slowing down of the Amazon forest due to increased drought occurrences':

- script_paper_calc_tac.R: R script to detrend EVI and kNDVI time series, calculate lag-1 autocorrelation (TAC), and calculate linear trend in TAC over time period.
- script_paper_droughtlegacy_models.R: R script to calculate droughtlegacy variables per pixel, and to run spatial simultaneous autoregressive lag model over different drought time periods, different moving window lenghts, and different regions within the Amazon.
- script_paper_vod_comparison.R: R script to compare TAC results from EVI and VOD time series.
- script_paper_kndvi_comparison.R: R script to calculate kNDVI time series and to compare TAC results from EVI and kNDVI time series.
- script_paper_plot_figures.R: R script to plot figures used in paper.
