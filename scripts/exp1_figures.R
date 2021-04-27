# Oculomics analysis, visualizations for proteomics experiment 1
# jinquma/ cgates 
# 2/10/2021

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis')
source('scripts/common.R')
script_name = "exp1_figures" 
log_start(paste0(output_logs_dir, script_name, ".log")) 


library(GGally)
library(RColorBrewer)
library(tidyverse)

###############################################################################
## Load normalized data

#setwd("/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis")
raw_df = readxl::read_xlsx("inputs/UM_F_50cm_2019_0657_0664.xlsx", sheet = 2)
normalized_data_df=raw_df[grepl("*Normalized*", colnames(raw_df))]
colnames(normalized_data_df)=sapply(strsplit(colnames(normalized_data_df), split = ","), function(x){x[[3]]})
normalized_data_df <- normalized_data_df %>% rename(
  `Pool1.1`=" Pool 1",
  `Pool1.2`=" Pool 2",
  `Pool1.3`=" Pool 3",
  `Pool1.4`=" Pool 4",
  `Pool1.5`=" Pool 5")
normalized_data_df=na.omit(normalized_data_df)
normalized_samples_df = normalized_data_df[1:5]
normalized_pools_df = normalized_data_df[6:10]

###############################################################################
## Pairwise Spearman correlation
# Considered Pearson and Spearman correlations; they are quite similar but 
# there are a few high-leverage outliers, so Spearman seemed a more conservative
# option.

subpheno_colors = colors_subphenotype

p_ <- GGally::print_if_interactive
grid.draw.gg <- function(x){
  print(x)
}


# A custom upper function that simplifies the correlation matrix to suppress
# the "Corr:" label
upper_fn <- function(
  data,
  mapping,
  ...,
  stars = TRUE,
  method = "pearson",
  use = "complete.obs",
  display_grid = FALSE,
  digits = 2,
  title_args = list(...),
  group_args = list(...),
  justify_labels = "right",
  align_percent = 0.5,
  title = "") {

  na.rm <-
    if (missing(use)) {
      # display warnings
      NA
    } else {
      (use %in% c("complete.obs", "pairwise.complete.obs", "na.or.complete"))
    }
  
  ggally_statistic(
    data = data,
    mapping = mapping,
    na.rm = na.rm,
    align_percent = align_percent,
    display_grid = FALSE,
    title_args = title_args,
    group_args = group_args,
    justify_labels = justify_labels,
    justify_text = "left",
    sep = "", 
    title = title,
    text_fn = function(x, y) {

      # set exact=FALSE to avoid warning "Cannot compute exact p-value with ties"
      corObj <- stats::cor.test(x, y, method = method, use = use, exact=FALSE)
      
      # make sure all values have X-many decimal places
      cor_est <- as.numeric(corObj$estimate)
      cor_txt <- formatC(cor_est, digits = digits, format = "f")
      
      cor_txt
    }
  )
}


# Ensures the scales are consistent across subplots
lower_fun <- function(data,mapping, min_all, max_all){
  ggplot(data = data, mapping = mapping)+
    geom_point(alpha = 0.3, color=subpheno_colors[2])+
    scale_x_continuous(limits = c(min_all, max_all))+
    scale_y_continuous(limits = c(min_all, max_all))
}


# Just print the sample names on the diag
diag_fn <- function (data, mapping) {
  label <- gsub('`', '', mapping_string(mapping$x))
  p <- ggally_text(label = label, color=subpheno_colors[2], size=8) + theme_void()
  p
}


max_all = max(normalized_samples_df)
min_all = min(normalized_samples_df)
p = ggpairs(normalized_samples_df,
        lower = list(continuous = wrap(lower_fun, min_all=min_all, max_all=max_all)),
        upper = list(continuous = wrap(upper_fn, method= "spearman", size=8, color='black')),
        diag = list(continuous = diag_fn),
        axisLabels='none') +
  theme(text=element_text(size=20),
        panel.grid.minor = element_blank(),
        title=element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "grey", fill=NA, size=1))

p_(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-pairwise_scatter_samples.pdf'), 
       plot=p, 
       width=10, height=10)

max_all = max(normalized_pools_df)
min_all = min(normalized_pools_df)
p = ggpairs(normalized_pools_df,
        lower = list(continuous = wrap(lower_fun, min_all=min_all, max_all=max_all)),
        upper = list(continuous = wrap(upper_fn, method= "spearman", size=8, color='black')),
        diag = list(continuous = diag_fn),
        axisLabels='none') +
  theme(text=element_text(size=20),
        panel.grid.minor = element_blank(),
        title=element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "grey", fill=NA, size=1))

p_(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-pairwise_scatter_pools.pdf'), 
       plot=p, 
       width=10, height=10)

sessionInfo()

print('Done')
###############################################################################

###############################################################################
## Technical vs Biological variation
# A pretty visualization, but ultimately not used in publication

# a=apply(norm_data[1:5], 1, var)
# b=apply(norm_data[6:10], 1, var)
# a=data.frame(a)
# a$sample="BioTechVar"
# colnames(a)[1]="Var"
# 
# b=data.frame(b)
# b$sample="TechVar"
# colnames(b)[1]="Var"
# temp=rbind(a, b)
# mu=temp %>%
#   group_by(sample) %>%
#   summarise(mean=mean(log(Var)))
# ggplot(temp, aes(log(Var), fill=sample))+geom_histogram(alpha=0.8)
# p=ggplot(temp, aes(log(Var), color=sample))+geom_histogram(aes(fill=sample), alpha=0.3, position = "identity")
# p+geom_vline(data=mu, aes(xintercept=mean, color=sample), linetype="dashed")
# 
# sessionInfo()
# print('done')
                                                                                 