# Oculomics power analysis visualizations and data
# jinquma/ cgates 
# 2/1/2021

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis')
source('scripts/common.R')
script_name = "exp2_power" 
log_start(paste0(output_logs_dir, script_name, ".log")) 

library(effsize) #https://www.rdocumentation.org/packages/effsize/versions/0.8.1/topics/cohen.d
library(RColorBrewer)
library(ggrepel)
library(ssize.fdr)
library(tidyverse)


###############################################################################
## Setup data structures
###############################################################################

###############################################################################
## Load sample phenotype design
# Focus to experiment 2, non-pool samples
batch="B"
design_df = read.csv("inputs/design.csv") 
design_df = design_df[design_df$batch==batch,] %>% 
  filter(!str_detect(sample, 'Pool'))
design_df$phenotype <- recode(design_df$phenotype, 
                           `PDR-Clear` = "PDR-L",
                           `PDR-Yellow` = "PDR-M",
                           `PDR-Red` = "PDR-H")


###############################################################################
## Load normalized abundance and design
# Drop NAs to create dense matrix
normalized_abundance_df = read.table("inputs/Report_ratio_groupby=gene_protNorm=MD_gu=0.tsv", stringsAsFactors = F, header = T, row.names =1)
normalized_abundance_df = normalized_abundance_df %>%
  drop_na() %>% 
  select(design_df$sample)

# Ensure abundance column names and phenotype row names are consistent
design_df= design_df[match(colnames(normalized_abundance_df), design_df$sample ),]

###############################################################################
## Model diffex across groups
#design_matrix = model.matrix(~as.factor(as.character(design_df$pheno)))
#fit = lmFit(normalized_abundance_df, design_matrix)
#ebayes_results = eBayes(fit)
#diffex_results = topTable(ebayes_results, number =Inf)

diffex_df<-read.table(paste0(output_data_dir, "exp2_diffex-kammers.tsv"), 
                      stringsAsFactors = F,
                      header = T,
                      row.names =14)
head(diffex_df)

###############################################################################
## Calculate the Cohen's d (Hedge-adjusted) effect size, pooled sd, and delta
phenotypes = factor(design_df$pheno)
cohen_hedges_fn = function (data_row) {
  cohen_d = cohen.d(as.matrix(data_row), phenotypes, hedges.correction=TRUE, conf.level=0.95)
  list(delta_mean=abs(cohen_d$estimate * cohen_d$sd),
       pooled_sd=cohen_d$sd,
       effect_size=cohen_d$estimate,
       confidence_interval_95_lower=cohen_d$conf.int[1],
       confidence_interval_95_upper=cohen_d$conf.int[2]
       )
}


effect_size_df = data.frame(bind_rows(apply(normalized_abundance_df, 1, cohen_hedges_fn)), row.names=rownames(normalized_abundance_df))
effect_size_df$abs_effect_ntile = ntile(abs(effect_size_df$effect_size), 100)

###############################################################################
## Export supplemental data table
supplemental_data_export_df = rownames_to_column(effect_size_df, 'gene') %>% inner_join(rownames_to_column(diffex_df, 'gene'), by='gene')
write.table(supplemental_data_export_df, 
            paste0(output_data_dir, script_name, '-gene_effect_size_diffex.tsv'),
            row.names = FALSE,
            col.names = TRUE)


###############################################################################
## Plotting
###############################################################################

## Setup standard colors
color_paired_palette = brewer.pal(n = 8, name = 'Paired')
color_pink = color_paired_palette[5]
color_blue = color_paired_palette[2]

color_control = colors_subphenotype[2]
color_pdr_l = colors_subphenotype[3]
color_pdr_m = colors_subphenotype[4]
color_pdr_h = colors_subphenotype[5]

## Setup a basic df to drive plotting
plot_df = data.frame(effect_size_df)
plot_df$abs_effect_size = abs(plot_df$effect_size)
plot_df$abs_confidence_interval_95_lower = ifelse(plot_df$effect_size<0, -1 * plot_df$confidence_interval_95_upper, plot_df$confidence_interval_95_lower)
plot_df$abs_confidence_interval_95_upper = ifelse(plot_df$effect_size<0, -1 * plot_df$confidence_interval_95_lower, plot_df$confidence_interval_95_upper)


plot_df$gene = rownames(plot_df)
plot_df = plot_df %>% inner_join(rownames_to_column(diffex_df[,c('p.mod','q.mod')], 'gene'), by='gene')
plot_df$inverse_pooled_sd = 1/plot_df$pooled_sd
plot_df = plot_df %>% 
   mutate(adjusted_p_value_f =  case_when(
     `q.mod` <= 0.01 ~ '<= 0.01',
     `q.mod` <= 0.05 ~ '<= 0.05',
     `q.mod` <= 0.10 ~ '<= 0.10',
      TRUE ~ '> 0.10'))
plot_df$adjusted_p_value_f = factor(plot_df$adjusted_p_value_f, 
                              levels=c("<= 0.01", "<= 0.05", "<= 0.10", "> 0.10"))
plot_df=plot_df[order(plot_df$effect_size, plot_df$gene),]

plot_df$selected_gene = plot_df$gene %in% c("CA2", "SERPINA1", "SPINK1", "NEO1")

plot_df$gene_label = ""
plot_df[plot_df$selected_gene,]$gene_label = plot_df[plot_df$selected_gene,]$gene


# HEY: be careful because while these factors are useful for plotting, you have 
# to be deliberate and add them if you want to talk about proportions.
apv_t = table(plot_df$adjusted_p_value_f)
addmargins(apv_t)
prop.table(apv_t)

###############################################################################
## exp2_pvalue_histogram.pdf
p = ggplot(plot_df, aes(x=p.mod)) + 
  geom_histogram(bins=100) +
  scale_x_continuous(breaks=c(0.01,0.05,0.1,0.25, 0.5, 1)) +
  theme_classic() +
  theme(text = element_text(size=30), 
        axis.text.x = element_text(angle = 90, vjust=0.5)) +
  xlab("moderated p-value")

plot(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-pvalue_histogram.pdf'),
       plot=p, width=10, height=10)

###############################################################################
## exp2_adjusted_pvalue_histogram.pdf
colors <- c(rep(colors_sequential_greens[1],1),
            rep(colors_sequential_greens[2],4), 
            rep(colors_sequential_greens[3],5),
            rep(colors_sequential_greens[4], 90))

p = ggplot(plot_df, aes(x=q.mod)) + 
  geom_histogram(bins=100, aes(fill=adjusted_p_value_f)) +
  scale_fill_manual(values = colors_sequential_greens) +
  geom_histogram(bins=100, fill=colors) +
  scale_x_continuous(breaks=c(0.01,0.05,0.1,0.25, 0.5, 1)) +
  theme_classic() +
  theme(text = element_text(size=30), 
        axis.text.x = element_text(angle = 90, vjust=0.5),
        legend.position = c(0.70, 0.50),) +
  labs(fill = "adjusted p-value") +
  xlab("Storey adjusted moderated p-value")

plot(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-adjusted_pvalue_histogram.pdf'),
       plot=p, 
       width=10, height=10)


###############################################################################
## exp2_effect_size_scatter.pdf
p = ggplot(plot_df, aes(x=delta_mean, y=inverse_pooled_sd, color=selected_gene, label=gene_label)) + 
  geom_point(data=subset(plot_df, !selected_gene), shape=19, size=4, color=alpha("grey", 0.5)) +
  geom_point(data=subset(plot_df, selected_gene), shape=19, size=5, color=color_blue) +
  geom_text_repel(color=color_blue, size=10, seed=1, point.padding=.75, min.segment.length = Inf) +
  theme_classic() +
  theme(text = element_text(size=30)) + 
  ylab("1/pooled sd") +
  xlab("delta PDR vs. control")

plot(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-effect_size_scatter.pdf'),
       plot=p,
       width=10, height=10)

###############################################################################
## exp2_power_estimated_effect_size_lineplot.pdf

temp_df = plot_df[order(plot_df$abs_effect_size, plot_df$gene),] 
temp_df$x = seq.int(nrow(temp_df))

tail(temp_df,1)

p = ggplot(temp_df, 
           aes(x=x, y=abs_effect_size, color=adjusted_p_value_f)) + 
  scale_color_manual(breaks=c('<= 0.01', '<= 0.05', '<= 0.10', '> 0.10'),
                    values = colors_sequential_greens) +
  geom_ribbon(aes(ymin=abs_confidence_interval_95_lower, ymax=abs_confidence_interval_95_upper), fill='gray20', alpha=0.2, colour=NA) +
  geom_point(size=5) + 
  geom_point(data=subset(temp_df, selected_gene), aes(y=abs_effect_size+0.20), color=color_blue, fill=color_blue, shape=25, size=5) + 
  theme_classic() +
  theme(legend.position = c(0.40, 0.90),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text=element_text(size=30)) + 
  geom_hline(yintercept=0) +
  labs(color='adjusted p value') +
  ylab("estimated effect size") +
  xlab("all proteins") +
  geom_text(data=subset(temp_df, selected_gene),
            aes(y=abs_confidence_interval_95_upper, label=gene),
            nudge_y=.2,
            color=color_blue, 
            size=10,
            angle=90,
            hjust=0)

plot(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-estimated_effect_size_lineplot.pdf'),
       plot=p,
       width=10, height=10)

###############################################################################
## exp2_select_genes_phenotype_boxplots.pdf
gene_sample_values = as_tibble(t(normalized_abundance_df[plot_df$gene[plot_df$selected_gene],]), rownames='sample') %>% 
  inner_join(design_df[,c("sample", "pheno")])

tall_gene_sample_values = cbind(gene_sample_values[c("sample", "pheno")], stack(gene_sample_values[plot_df$gene[plot_df$selected_gene]]))
colnames(tall_gene_sample_values) = c('sample', 'pheno', 'values', 'gene')

tall_gene_sample_values$gene_f = factor(tall_gene_sample_values$gene, levels=c('SERPINA1', 'NEO1', 'SPINK1', 'CA2'))
tall_gene_sample_values$phenotype = factor(tall_gene_sample_values$pheno, levels=c('PDR', 'Control'))

p = ggplot(tall_gene_sample_values, aes(x=phenotype, y=values, fill = phenotype)) +
  geom_violin(linetype='blank') +
  stat_summary(fun.data=mean_sdl, geom="pointrange") +
  facet_grid( ~ gene_f) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position="bottom",
        legend.title=element_blank()) +
  ylab('normalized protein abundance') +
  labs(fill='sample phenotype')

plot(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-select_genes_phenotype_boxplots.pdf'), 
       plot=p,
       width=10, height=10)

###############################################################################
## exp2_power_select_genes_qqplots.pdf
p = ggplot(tall_gene_sample_values, aes(sample=values, color=phenotype)) + 
  stat_qq_line(color='blue', alpha=0.5, size=2) +
  stat_qq(size=4, alpha=0.7) +
  facet_grid(rows=vars(gene_f), cols=vars(phenotype), scales='fixed') +
  theme_classic() +
  theme(text = element_text(size=30), 
        legend.position="none",
        panel.spacing = unit(2, "lines")) + 
  labs(x = "normal theoretical quantiles", y = "empirical quantiles")

plot(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-select_genes_qqplots.pdf'), 
       plot=p,
       width=10, height=10)


temp_df= tall_gene_sample_values[tall_gene_sample_values$gene=='CA2' & tall_gene_sample_values$phenotype=='PDR',]
shapiro.test(temp_df$values)

###############################################################################
## exp2_select_genes_subphenotype_boxplots.pdf
gene_sample_values = as_tibble(t(normalized_abundance_df[plot_df$gene[plot_df$selected_gene],]), rownames='sample') %>% 
  inner_join(design_df[,c("sample", "phenotype")])

tall_gene_sample_values = cbind(gene_sample_values[c("sample", "phenotype")], stack(gene_sample_values[plot_df$gene[plot_df$selected_gene]]))
colnames(tall_gene_sample_values) = c('sample', 'phenotype', 'values', 'gene')

tall_gene_sample_values$gene_f = factor(tall_gene_sample_values$gene, levels=c('SERPINA1', 'NEO1', 'SPINK1', 'CA2'))
tall_gene_sample_values$phenotype = factor(tall_gene_sample_values$phenotype, levels=c('PDR-H', 'PDR-M', 'PDR-L', 'Control'))
pheno_colors = c(color_pdr_h, color_pdr_m, color_pdr_l, color_control)

p = ggplot(tall_gene_sample_values, aes(x=phenotype, y=values, fill = phenotype)) +
  geom_violin(linetype='blank') +
  stat_summary(fun.data=mean_sdl, geom="pointrange") +
  scale_fill_manual(values=pheno_colors) +
  facet_grid( ~ gene_f) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position="bottom",
        legend.title=element_blank()) +
  ylab('normalized protein abundance') +
  labs(fill='sample phenotype')

plot(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-select_genes_subphenotype_boxplots.pdf'),
       plot=p,
       width=10, height=10)



###############################################################################
## exp2_prospective_power_by_effect_size.pdf
named_genes = c('CA2', 'SPINK1', 'NEO1')
power_fn <- function(x) {
  fdr = 0.05
  desired_power = 0.8
  pi0 = 0.7
  max_samples = 25
  power=ssize.twoSamp(delta = as.double(x[2]), 
                      sigma = as.double(x[3]),
                      fdr,
                      desired_power,
                      pi0,
                      max_samples,
                      "two-sided")$power
  power=as.data.frame(power)
  power$gene = x[1]
  colnames(power) = c('n', 'power', 'gene')
  power
}

temp = plot_df[plot_df$gene %in% named_genes, c('gene', 'delta_mean', 'pooled_sd', 'abs_effect_ntile')]
rownames(temp) = temp$gene
temp$gene_factor = paste0(temp$gene, " (", temp$abs_effect_ntile, "%)")

power_df = do.call(rbind, apply(temp[, c('gene', 'delta_mean', 'pooled_sd')], 1, power_fn)) 
power_df = power_df %>% inner_join(temp[,c('gene', 'gene_factor')])
power_df$gene_factor=factor(power_df$gene_factor, levels = temp[named_genes, 'gene_factor'])

p = ggplot(power_df, aes(x=n, y=power, color=gene_factor)) +
  theme_classic() +
  theme(text = element_text(size=30), 
        legend.position = c(0.75, 0.30), 
        legend.title=element_text(size=20),
        legend.text=element_text(size=20, margin=margin(5,0,5,0, 'pt'))) +
  geom_hline(yintercept=0.8, colour = "gray", size=2, alpha=1) +
  geom_vline(xintercept=8, color="#F8766D", size=2, alpha=0.5) +
  geom_vline(xintercept=11, color="#00BA38", size=2, alpha=.5) +
  geom_vline(xintercept=14, color="#619CFF", size=2, alpha=.5) +
  geom_point(size=3) +
  geom_line(aes(group=gene), size=2) +
  labs(color="Protein (effect size percentile)") +
  scale_y_continuous(breaks=seq(0,1, 0.2)) +
  scale_x_continuous(name="minimum samples per group")


plot(p)

ggsave(filename=paste0(output_figures_dir, script_name, '-prospective_power_by_effect_size.pdf'), 
       plot=p, width=10, height=10)

print('Done')
###############################################################################
log_stop()

