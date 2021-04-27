# Oculomics analysis, visualizations for proteomics experiment 2
# jinquma/ cgates 
# 2/22/2021

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis')
source('scripts/common.R')
script_name = "exp2_upset" 
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
## Load normalized abundance and design
dat_gene<-read.table("inputs/Report_ratio_groupby=gene_protNorm=MD_gu=0.tsv", stringsAsFactors = F, header = T, row.names =1)

dat_gene <- dat_gene %>% rename(
  Pool1.6=Pool12,
  Pool1.7=Pool13,
  Pool1.8=Pool14,
  Pool1.9=Pool15
)

#view(dat_gene)

###############################################################################
## Load sample phenotype design
# Focus to experiment 2 (batch B)
pheno=read.csv("inputs/design.csv")

pheno$sample <- recode(pheno$sample, 
                       `Pool12`="Pool1.6",
                       `Pool13`="Pool1.7",
                       `Pool14`="Pool1.8",
                       `Pool15`="Pool1.9")
#pheno$sample

pheno$plex <- recode(pheno$plex, 
                     `0`="1",
                     `1`="2.1",
                     `2`="2.2",
                     `3`="2.3",
                     `4`="2.4")
#pheno$plex

batch="B"
myplex=as.character(pheno$sample[pheno$batch==batch])
#view(myplex)
all_abundances = dat_gene[names(dat_gene) %in% myplex]

design=pheno[pheno$sample %in% colnames(all_abundances),]
#Reset phenotype names
design$phenotype <- recode(design$phenotype, 
                           `PDR-Clear` = "PDR-L",
                           `PDR-Yellow` = "PDR-M",
                           `PDR-Red` = "PDR-H")
design=design[match(colnames(all_abundances), design$sample ),]


###############################################################################
## Setup dfs to drive plots

###
temp_plex_df = all_abundances %>% 
  as_tibble(rownames=NA) %>%
  rownames_to_column(var='protein') %>%
  pivot_longer(cols=!protein, names_to='sample', values_to='protein_abundance') %>%
  inner_join(design[,c('sample', 'plex')], by='sample') %>%
  group_by(protein, plex) %>%
  summarize(measured_samples=sum(!is.na(protein_abundance))) %>%
  mutate(plex_dropout_full = measured_samples == 0,
         plex_dropout_partial = measured_samples > 0 & measured_samples < 9)

#temp_plex_df

temp_plex_dropout_full_df = temp_plex_df[,c('protein', 'plex', 'plex_dropout_full')] %>% 
  pivot_wider(names_from='plex', names_prefix='plex_dropout_full_', values_from='plex_dropout_full') %>%
  mutate(plex_dropout_full=rowSums(across(where(is.logical)))>0)

temp_plex_dropout_partial_df = temp_plex_df[,c('protein', 'plex', 'plex_dropout_partial')] %>% 
  pivot_wider(names_from='plex', names_prefix='plex_dropout_partial_', values_from='plex_dropout_partial') %>%
  mutate(plex_dropout_partial=rowSums(across(where(is.logical)))>0)

temp_plex_dropout_df = temp_plex_dropout_full_df %>% 
  inner_join(temp_plex_dropout_partial_df) %>%
  select(plex_dropout_full, plex_dropout_partial) %>%
  mutate(dropout=case_when(plex_dropout_partial==TRUE ~ 'intraplex dropout',
                           plex_dropout_full==TRUE ~ 'interplex dropout',
                           TRUE ~ 'no dropout')) %>%
  select(protein, dropout)

temp_plex_dropout_df$dropout = factor(temp_plex_dropout_df$dropout, 
                                      levels=c('no dropout', 'interplex dropout', 'intraplex dropout'))

temp_df = all_abundances %>%
  as_tibble(rownames=NA) %>%
  rownames_to_column(var='protein') %>%
  mutate_at(vars(-protein), negate(is.na)) %>%
  group_by(vars(-protein)) %>%
  inner_join(temp_plex_dropout_df) %>%
  select(-protein)
#temp_df


temp_sample_df = design[,c('sample', 'phenotype', 'plex')]
temp_sample_df$phenotype_f = factor(temp_sample_df$phenotype,
                                    levels=c('Pool1', 'Control', 'PDR-L', 'PDR-M', 'PDR-H'))
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
make_upset_plot_fn = function (df, cols, intersection_breaks, set_size_limits_y, set_size_breaks, text_size) {
  
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
                counts=FALSE,
                mapping=aes(fill=dropout),
                show.legend=show_intersection_legend) + 
              scale_fill_manual(values=c(
                'no dropout'='purple3',
                'interplex dropout'='black',
                'intraplex dropout'='gray'), 'protein dropout') + 
              scale_y_continuous(name ="proteins in sample\nintersection", breaks=intersection_breaks) +
              geom_point(aes(y=exclusive_intersection_size, color=dropout), size=3, shape=18, show.legend=FALSE) +
              scale_color_manual(values=c(
                'no dropout'='purple3',
                'interplex dropout'='black',
                'intraplex dropout'='gray')) + 
              guides(color=FALSE)
          ),
          set_sizes = upset_set_size(geom=geom_point(stat='count', size=5, shape=18)) +
            ylab('measured proteins\nfor each sample') +
            expand_limits(y=set_size_limits_y) + 
            scale_y_continuous(breaks=set_size_breaks),
          height_ratio=2,
          width_ratio=.25,
          wrap=FALSE,
          sort_sets=FALSE,
          stripes=c(alpha('grey90', 0.45), alpha('white', 0.3)),
          themes = upset_modify_themes(list(
            'Intersection size'=theme(text=element_text(size=text_size), 
                                      axis.title.y = element_text(vjust=-10)),
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
make_annotation_strip_fn = function (df, col, fill, fill_palette, legend_name, text_size) {
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

text_size = 12

###############################################################################
## Exp2: Make upset plot sorted by phenotype (exp2_upset_by_phenotype.pdf)

intersection_breaks = c(0,100,727)
set_size_limits_y = c(875,1000)
set_size_breaks = c(875, 900, 925, 950, 975, 1000)
upset_plot = make_upset_plot_fn(temp_df, 
                                temp_sample_df$sample_f,
                                intersection_breaks,
                                set_size_limits_y,
                                set_size_breaks,
                                text_size)

#plot(upset_plot$plot)

# Add a line at break between pool|control and control|PDR
upset_plot$intersection_matrix = upset_plot$intersection_matrix + 
  geom_hline(yintercept = 32.5, color=alpha('cyan', 0.3), size=3)

upset_plot$intersection_matrix = upset_plot$intersection_matrix + 
  geom_hline(yintercept = 22.5, color=alpha('red', 0.3), size=3)

plex_annotation_strip = make_annotation_strip_fn(temp_sample_df, 'sample_f', 'plex_f', colors_sequential_greens, 'TMT plex', text_size)
pheno_annotation_strip = make_annotation_strip_fn(temp_sample_df, 'sample_f', 'phenotype_f', colors_subphenotype, 'phenotype', text_size)


###############################################################################
## Assemble plot components

layout <- '
#A##X
BCDEY
'

plot_elements_by_pheno = list(intersection_barplot = upset_plot$intersection_barplot,
                             measured_proteins_barplot = upset_plot$measured_proteins_barplot,
                             intersection_matrix = upset_plot$intersection_matrix,
                             pheno_annotation_strip = pheno_annotation_strip$plot,
                             plex_annotation_strip = plex_annotation_strip$plot,
                             intersection_barplot_legend = upset_plot$legend)

final_plot = wrap_plots(A = upset_plot$intersection_barplot,
           B = upset_plot$measured_proteins_barplot,
           C = upset_plot$intersection_matrix,
           D = pheno_annotation_strip$plot,
           E = plex_annotation_strip$plot,
           X = upset_plot$legend,
           Y = plot_grid(pheno_annotation_strip$legend, plex_annotation_strip$legend, ncol = 1, align='v'),
           design = layout,
           widths=c(3, 10, 0.3, 0.3, 3),
           heights=c(5,10))

plot(final_plot)

ggsave(filename=paste0(output_figures_dir, script_name, '-by_phenotype.pdf'), 
       plot=final_plot, 
       width=17, height=15, scale=0.75)



###############################################################################
## Make upset plot as above but sorted by plex

sample_plex_df = data.frame(temp_sample_df)
sample_plex_df$sample_f = factor(sample_plex_df$sample,
                                 levels=sample_plex_df$sample[order(sample_plex_df$plex_f, sample_plex_df$phenotype_f, sample_plex_df$sample)])
sample_plex_df = sample_plex_df[order(sample_plex_df$sample_f),]


intersection_breaks = c(0,100,727)
set_size_limits_y = c(875,1000)
set_size_breaks = c(875, 900, 925, 950, 975, 1000)
upset_plot = make_upset_plot_fn(temp_df, 
                                sample_plex_df$sample_f,
                                intersection_breaks,
                                set_size_limits_y,
                                set_size_breaks,
                                text_size)

#plot(upset_plot$plot)

# Highlihght the singleton plex sets
df <- as.data.frame(
      matrix(c(3, 3, 1, 9, 
               6, 6, 19, 27,
               7, 7, 28, 36,
               8, 8, 10, 18),
             nrow = 4, byrow=TRUE))
colnames(df) = c('x1', 'x2', 'y1', 'y2')
#df
upset_plot$intersection_matrix =   upset_plot$intersection_matrix + 
  geom_hline(yintercept = 9.5, color=alpha('green', 0.3), size=3) +
  geom_hline(yintercept = 18.5, color=alpha('green', 0.3), size=3) +
  geom_hline(yintercept = 27.5, color=alpha('green', 0.3), size=3) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, size=3), color=alpha('cyan', 0.3), data = df, show.legend=FALSE)

  
#plot(upset_plot$intersection_matrix)

plex_annotation_strip = make_annotation_strip_fn(sample_plex_df, 'sample_f', 'plex_f', colors_sequential_greens, 'TMT plex', text_size)
pheno_annotation_strip = make_annotation_strip_fn(sample_plex_df, 'sample_f', 'phenotype_f', colors_subphenotype, 'phenotype', text_size)

layout <- '
#A##X
BCDEY
'
plot_elements_by_plex = list(intersection_barplot = upset_plot$intersection_barplot,
                             measured_proteins_barplot = upset_plot$measured_proteins_barplot,
                             intersection_matrix = upset_plot$intersection_matrix,
                             pheno_annotation_strip = pheno_annotation_strip$plot,
                             plex_annotation_strip = plex_annotation_strip$plot,
                             intersection_barplot_legend = upset_plot$legend)

final_plot = wrap_plots(A = upset_plot$intersection_barplot,
                        B = upset_plot$measured_proteins_barplot,
                        C = upset_plot$intersection_matrix,
                        D = plex_annotation_strip$plot,
                        E = pheno_annotation_strip$plot,
                        X = upset_plot$legend,
                        Y = plot_grid(pheno_annotation_strip$legend, plex_annotation_strip$legend, ncol = 1, align='v'),
                        design = layout,
                        widths=c(3, 10, 0.3, 0.3, 3),
                        heights=c(5,10))

plot(final_plot)

ggsave(filename=paste0(output_figures_dir, script_name, '-by_plex.pdf'),
       plot=final_plot,
       width=17, height=15, scale=0.75)



###############################################################################
## Make composite upset plot combining both plots above

layout <- '
#A##X
BCDEY
FGHI#
'

final_plot = wrap_plots(A = plot_elements_by_plex$intersection_barplot,
                        B = plot_elements_by_plex$measured_proteins_barplot + theme(axis.title.x=element_blank()),
                        C = plot_elements_by_plex$intersection_matrix,
                        D = plot_elements_by_plex$plex_annotation_strip,
                        E = plot_elements_by_plex$pheno_annotation_strip,
                        F = plot_elements_by_pheno$measured_proteins_barplot,
                        G = plot_elements_by_pheno$intersection_matrix,
                        H = plot_elements_by_pheno$plex_annotation_strip,
                        I = plot_elements_by_pheno$pheno_annotation_strip,
                        X = plot_elements_by_plex$intersection_barplot_legend,
                        Y = plot_grid(pheno_annotation_strip$legend, plex_annotation_strip$legend, ncol = 1, align='v'),
                        design = layout,
                        widths=c(3, 10, 0.3, 0.3, 3),
                        heights=c(5,10,10))

plot(final_plot)

ggsave(filename=paste0(output_figures_dir, script_name, '-composite.pdf'),
       plot=final_plot,
       width=17, height=25, scale=0.75)


###############################################################################
## Summary stats for figures
###############################################################################

temp_df
# How many proteins measured by each sample
summary(colSums(temp_df[,temp_sample_df$sample_f]))

#Count if different protein dropout categories
summary(temp_df$dropout)

# Proteins not measured in pools
all_abundances[rowSums(temp_df[,c('Pool1.6', 'Pool1.7', 'Pool1.8', 'Pool1.9')])==0, ]

# Proteins measured consistently by exactly one plex 
singleton_plex_proteins = colSums(temp_plex_dropout_full_df[(rowSums(temp_plex_dropout_full_df[,2:5]) == 1), 2:5])
singleton_plex_proteins
sum(singleton_plex_proteins)


print('Done')
###############################################################################
log_stop()

