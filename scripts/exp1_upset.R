# Oculomics analysis, visualizations for proteomics experiment 2
# jinquma/ cgates 
# 2/22/2021

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis')
source('scripts/common.R')
script_name = "exp1_upset"
log_start(paste0(output_logs_dir, script_name, ".log")) 

library(cowplot)
library(GGally)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
if (!require(devtools, quietly=T)) install.packages("devtools")
if (!require(ComplexUpset, quietly=T)) devtools::install_github("krassowski/complex-upset")
library(ComplexUpset) #https://krassowski.github.io/complex-upset/articles/Examples_Python.html#11-advanced-usage-examples


###############################################################################
## Setup data structures
###############################################################################

###############################################################################
raw_df = readxl::read_xlsx("inputs/UM_F_50cm_2019_0657_0664.xlsx", sheet = 2)
normalized_data_df=raw_df[grepl("*Normalized*", colnames(raw_df))]
colnames(normalized_data_df)=sapply(sapply(strsplit(colnames(normalized_data_df), split = ","), function(x){x[[3]]}), trimws)
normalized_data_df <- normalized_data_df %>% rename(
  `Pool1.1`="Pool 1",
  `Pool1.2`="Pool 2",
  `Pool1.3`="Pool 3",
  `Pool1.4`="Pool 4",
  `Pool1.5`="Pool 5")


###############################################################################
## Load sample phenotype design
# Focus to experiment 2 (batch B)
pheno=read.csv("inputs/design.csv")

pheno$sample <- recode(pheno$sample, 
                       `Pool1-1`="Pool1.1",
                       `Pool1-2`="Pool1.2",
                       `Pool1-3`="Pool1.3",
                       `Pool1-4`="Pool1.4",
                       `Pool1-5`="Pool1.5")
#pheno$sample

batch="A"
myplex=as.character(pheno$sample[pheno$batch==batch])
#view(myplex)
all_abundances = normalized_data_df[names(normalized_data_df) %in% myplex]
#names(normalized_data_df)

design=pheno[pheno$sample %in% colnames(all_abundances),]
design=design[match(colnames(all_abundances), design$sample ),]


temp_df = all_abundances %>%
  as_tibble(rownames=NA) %>%
  mutate_all(negate(is.na))

# Exclude proteins which were not assigned an abundance in any sample
temp_df = temp_df[rowSums(temp_df)>0,]
  
temp_sample_df = design[,c('sample', 'phenotype', 'plex')]
temp_sample_df$phenotype_f = recode_factor(temp_sample_df$phenotype,
                                           `Pool1` = 'technical replicates',
                                           `Control` = 'biological replicates')
#head(temp_sample_df)
temp_sample_df$plex_f = factor(temp_sample_df$plex)
temp_sample_df$sample_f = factor(temp_sample_df$sample,
                                 levels=temp_sample_df$sample[order(temp_sample_df$phenotype_f, temp_sample_df$sample)])

temp_sample_df = temp_sample_df[order(temp_sample_df$sample_f),]
#temp_sample_df



###############################################################################
## Build plots helper functions
###############################################################################


## A function that accepts a dataframe, vector of columns, and text size and 
# returns a list of plots (composite plot along with all the plot components 
# and legend)
#
make_upset_plot_fn = function (df, cols, intersection_breaks, set_size_limits_y, text_size) {
  
  # There is a subtle bug in the upset library (intersection matrix) which
  # assumes the dataframe column order matches the order cols and also assumes
  # the cols are not factors but a simple vector of colnames. The bug causes
  # the intersection matrix rows, set counts, and labels to be mismatched.
  # These next few lines evade that bug by explicitly reordering the columns.
  
  ordered_df = data.frame(df)
  ordered_df = temp_df[ , rev(cols)]
  ordered_cols = colnames(ordered_df)
  temp_extra_cols_df = df[, !names(df) %in% ordered_cols]
  ordered_df = cbind(ordered_df, temp_extra_cols_df)
  
  
  # It's hard to extract and remove the legend from the composite plot, so
  # this helper function allows me to generate it with the legend (so I can 
  # extract the legend) and without (so I can get at the component plots sans
  # legend).
  base_upset_plot_fn = function(show_intersection_legend=TRUE) {
    upset(ordered_df, 
          ordered_cols,
          name='', 
          base_annotations=list(
            'Intersection size'=intersection_size(
                counts=TRUE,
                text=list(size=7),
                show.legend=show_intersection_legend) + 
              scale_y_continuous(name ="proteins in sample\nintersection", breaks=intersection_breaks) +
              guides(color=FALSE)
          ),
          set_sizes = upset_set_size(geom=geom_point(stat='count', size=5, shape=18)) +
            ylab('measured proteins\nfor each sample') +
            expand_limits(y=set_size_limits_y) +
            scale_y_continuous(breaks=c(980, 989, 1000)), # + 
          height_ratio=2,
          width_ratio=.25,
          wrap=FALSE,
          sort_sets=FALSE,
          stripes=c(alpha('grey90', 0.45), alpha('white', 0.3)),
          themes = upset_modify_themes(list(
            'Intersection size'=theme(text=element_text(size=text_size), 
                                      axis.title.y = element_text(vjust=-1)),
            'intersections_matrix'=theme(text=element_text(size=text_size)),
            'overall_sizes'=theme(text=element_text(size=text_size), 
                                  axis.text.x=element_text(angle=90, vjust=0.5),
                                  axis.title.x = element_text(margin=margin(t=10)),
                                  panel.grid.major.x = element_line(colour = 'gray',size=0.75),
                                  panel.grid.major.y = element_blank(),
                                  plot.background=element_rect(fill='transparent', color=NA))
          ))
    )
  }
  
  
  legend = get_legend(base_upset_plot_fn(TRUE))
  upset_plot = base_upset_plot_fn(FALSE)
  return (list(plot=upset_plot,
               intersection_barplot = upset_plot[[2]],
               measured_proteins_barplot = upset_plot[[3]],
               intersection_matrix = upset_plot[[4]],
               legend=legend))
}



## Builds the sample annotation strips on side (ala pheatmap)
tation_strip_fn = function (df, col, fill, fill_palette, legend_name, text_size) {
  strip = ggplot(df) +
    geom_bar(mapping = aes(y=reorder(!! sym(col), desc(!! sym(col))), x=1, fill=!! sym(fill)),
             stat = "identity",
             width = 1) +
    scale_fill_manual(values = fill_palette, name = legend_name) +
    theme_void() +
    theme(text=element_text(size=text_size),
          plot.margin=margin(0,0,0,0))  
    
  legend = get_legend(strip)
  strip =  strip + theme(legend.position = "none")
  return (list(plot=strip, legend=legend))
}


###############################################################################
## Build the plots
###############################################################################

text_size = 20

###############################################################################
## Exp2: Make upset plot sorted by phenotype (exp2_upset_by_phenotype.pdf)

intersection_breaks = c(0,1000)
set_size_limits_y = c(980,1000)
upset_plot = make_upset_plot_fn(temp_df, 
                                temp_sample_df$sample_f,
                                intersection_breaks,
                                set_size_limits_y,
                                text_size)

#plot(upset_plot$plot)

pheno_annotation_strip = make_annotation_strip_fn(temp_sample_df, 'sample_f', 'phenotype_f', colors_subphenotype, '', text_size)
plot(pheno_annotation_strip$legend)
###############################################################################
## Assemble plot components

# This dumb repetition gives me a smidgen more control on where the legend is placed vertically.
# It's cheap and dumb, so if you see a better way, please make it suck less.
layout <- '
#A##
#A##
#A##
#A##
#A##
#A##
#A##
#A##
#A##
#A##
BCD#
BCD#
BCD#
BCD#
BCDE
BCD#
BCD#
BCD#
BCD#
BCD#
'

final_plot = wrap_plots(A = upset_plot$intersection_barplot,
           B = upset_plot$measured_proteins_barplot,
           C = upset_plot$intersection_matrix,
           D = pheno_annotation_strip$plot,
           E = pheno_annotation_strip$legend,
           design = layout,
           widths=c(2, 2, 0.125,2.5),
           heights=c(1,2))

plot(final_plot)

ggsave(filename=paste0(output_figures_dir,script_name,'-by_pool_vs_sample.pdf'), 
       plot=final_plot, 
       width=12, height=10, scale=0.75)


###############################################################################

sessionInfo()

print('Done')
###############################################################################
sink()
