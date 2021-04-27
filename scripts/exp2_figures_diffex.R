# Oculomics analysis, visualizations for proteomics experiment 2
# jinquma/ cgates 
# 3/5/2021

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis')
source('scripts/common.R')
script_name = "exp2_figures_diffex" 
log_start(paste0(output_logs_dir, script_name, ".log")) 

library(tidyverse)
library(RColorBrewer)

###############################################################################
## Setup data structures
###############################################################################

###############################################################################
## Load normalized abundance and design
diffex_df<-read.table(paste0(output_data_dir, "exp2_diffex-kammers.tsv"),
                      stringsAsFactors = F,
                      header = T,
                      row.names = 14)
head(diffex_df)


q_value_cutoff = 0.05
p_value_cutoff = max(diffex_df[diffex_df$q.mod <= q_value_cutoff,]$p.mod)

diffex_df[diffex_df$p.mod>0.154 & diffex_df$p.mod<0.155,]


erythrocyte_gene_symbols = c("ACTB", "HBG2", "CA2", "HBA1", "HBB", "HBD", "BLVRB", "BPGM", "PRDX2", "CA1", "CAT")

protein_accession = read.table("inputs/Report_ratio_groupby=gene_protNorm=MD_gu=0.tsv", stringsAsFactors = F, header = T, row.names =1) %>%
  select("Proteins") %>% 
  separate(Proteins, c(NA, "uniprot", NA),  sep="\\|", remove = TRUE)
plasma = readxl::read_xlsx("inputs/List_20 Plasma Proteins_Uniprot Accession.xlsx", skip = 2) %>%
  select("Accession (Uniprot)", "Depleted by Pierce Top12 Column (Yes/No)") %>%
  rename("uniprot"="Accession (Uniprot)", depleted="Depleted by Pierce Top12 Column (Yes/No)")
plasma_proteins = protein_accession %>% rownames_to_column(var='gene') %>% left_join(plasma) %>% drop_na()

###############################################################################
## Plotting
###############################################################################

logfc_cutoff=0.585
sig_cutoff=p_value_cutoff

# Add direction column
diffex_df$direction = 'NS'
diffex_df$direction[diffex_df$p.mod<=sig_cutoff & diffex_df$logFC <= 0] = 'Down'
diffex_df$direction[diffex_df$p.mod<=sig_cutoff & diffex_df$logFC > 0] = 'Up'

# Determine direction labels (with number per)
direction_table = table(diffex_df$direction)

# up_label = sprintf('Up: %s', direction_table['Up'])
# down_label = sprintf('Down: %s', direction_table['Down'])
# ns_label = sprintf('adj p-val > %0.2f: %s', q_value_cutoff, direction_table['NS'])

diffex_df$direction_count = factor(
  diffex_df$direction,
  levels = c('Up', 'Down', 'NS'),
  labels = c('Upregulated', 'Downregulated', sprintf('adj p-val > %0.2f', q_value_cutoff)))

# if(all( !( c('Up', 'Downl') %in% names(direction_table) ) )) {
#   warning(sprintf('No genes were DE at fdr < %s and |log2fc| > %s. Consider a different threshold.', fdr_cutoff, logfc_cutoff))
# }

# Transform qval to -log10 scale
diffex_df$log10_p_mod = -log10(diffex_df$p.mod)

diffex_df$protein_type = 'other'
diffex_df$protein_type[rownames(diffex_df) %in% erythrocyte_gene_symbols] = 'RBC'
diffex_df$protein_type[rownames(diffex_df) %in% plasma_proteins$gene] = 'plasma'
#table(diffex_df$protein_type)

# Add top 10 Up and 10 Down gene labels
# de_list is assumed to be ordered by Call/diff_exp and then qvalue from deseq2_diffex.R and ballgown_diffex.R
# top = rbind(
#   head(subset(de_list, direction == 'Up'), 10),
#   head(subset(de_list, direction == 'Down'), 10))
# top$label = top$external_gene_name
# de_list = merge(x = de_list, y = top[, c(id,'label')], by = id, all.x = TRUE, sort = FALSE)

# Volcano Plot
volcano_plot = ggplot(diffex_df, aes(x = logFC, y = log10_p_mod, color = direction_count, alpha = 0.5, shape=protein_type)) +
  geom_vline(xintercept = c(0, -1*logfc_cutoff, logfc_cutoff), linetype = c(1, 2, 2), color = c('darkgray', 'gray', 'gray'), size=1) +
  geom_hline(yintercept = -log10(sig_cutoff), linetype = 1, color = 'gray', size=1) +
  scale_alpha(guide = 'none') +
  geom_point(size = 3) +
  scale_shape_manual(name='protein type', values=c(16,9,7)) +
  scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray')) +
  labs(
    x = 'log2 fold-change',
    y = '-log10 p-value') +
  theme_classic() +
  theme(text = element_text(size=20))

plot(volcano_plot)

ggsave(filename = paste0(output_figures_dir, script_name, '-volcano.pdf'),
       plot = volcano_plot,
       height = 6, width = 8)

###############################################################################

diffex_df$logfc_cutoff = ifelse(abs(diffex_df$logFC) >= logfc_cutoff, '>=1.5 FC', '<1.5FC')
diffex_df$adj_pval_cutoff = ifelse(diffex_df$log10_p_mod > -log10(sig_cutoff), 'apv <= 0.05', 'apv > 0.05')

table(diffex_df$adj_pval_cutoff, diffex_df$logfc_cutoff)


cross_tab_df = as.data.frame.matrix(addmargins(table(diffex_df$adj_pval_cutoff, diffex_df$logfc_cutoff)))
cross_tab_df
write.table(cross_tab_df, 
            file = paste0(output_data_dir, script_name, "-diffex_crosstable.tsv"),
            sep = "\t", row.names = TRUE)
###############################################################################

sessionInfo()

print('Done')
###############################################################################
log_stop()