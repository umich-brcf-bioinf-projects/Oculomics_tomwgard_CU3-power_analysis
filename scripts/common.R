# Oculomics power analysis visualizations and data
# jinquma / cgates 
# 3/22/2021

library(RColorBrewer)

log_start <- function(log_file) {
  sink(log_file, append=FALSE, split=TRUE) 
  Sys.Date()
  getwd()
}

log_stop <- function() {
  sessionInfo()
  sink()
}

output_data_dir = "outputs/data/"
output_figures_dir = "outputs/figures/"
output_logs_dir = "outputs/logs/"

colors_sequential_greens = brewer.pal(n = 9, name = 'Greens')[c(9,7,5,3)]
colors_subphenotype = brewer.pal(n = 8, name = 'Paired')[c(1,2,7,8,6)]
