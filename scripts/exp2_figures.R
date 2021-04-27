# Oculomics analysis, visualizations for proteomics experiment 2
# jinquma/ cgates 
# 2/10/2021

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis')
source('scripts/common.R')
script_name = "exp2_figures" 
log_start(paste0(output_logs_dir, script_name, ".log")) 

library(GGally)
library(tidyverse)
library(RColorBrewer)


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
mydat = dat_gene[names(dat_gene) %in% myplex]
mydat = mydat %>% drop_na() # Drop NAs to create dense matrix
#colnames(mydat)

design=pheno[pheno$sample %in% colnames(mydat),]
#Reset phenotype names
design$phenotype <- recode(design$phenotype, 
       `PDR-Clear` = "PDR-L",
       `PDR-Yellow` = "PDR-M",
       `PDR-Red` = "PDR-H")
design=design[match(colnames(mydat), design$sample ),]

#rownames_to_column(head(mydat), var='protein')


###############################################################################
## Plotting
###############################################################################

###############################################################################
## normalized_abundance_boxplot.pdf

temp_df = rownames_to_column(mydat, var='protein') %>% 
  pivot_longer(cols=!protein, names_to='sample', values_to='protein_abundance') %>%
  inner_join(design, by='sample')
temp_df$sample_f = factor(temp_df$sample, 
                          levels=unique(temp_df$sample[order(temp_df$plex, temp_df$sample)])
                  )
#temp_df$sample_f = reorder(temp_df$plex, temp_df$sample_f)
levels(temp_df$sample_f)

temp_df$plex_f = factor(temp_df$plex, levels=c('2.1','2.2','2.3','2.4'), ordered=TRUE)
temp_df = temp_df[order(temp_df$plex_f, temp_df$sample_f),]
#temp_df = na.omit(temp_df)
plex_colors = colors_sequential_greens

p = ggplot(temp_df, aes(x=sample_f, y=protein_abundance, fill=plex_f)) +
  geom_boxplot(color=alpha('black', 0.7)) +
  scale_fill_manual(values=plex_colors) + 
  theme_classic() +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 90, vjust=0.5),
        legend.position="bottom") +
  xlab("sample") +
  ylab("normalized protein abundance") + 
  labs(fill="TMT plex")

plot(p)

ggsave(filename = paste0(output_figures_dir, script_name, '-normalized_abundance_boxplot.pdf'), 
       plot = p, 
       width = 10, height = 10)



###############################################################################
## exp2_pairwise_scatter_pools.pdf

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
    geom_point(alpha = 0.3, color=colors_subphenotype[2])+
    scale_x_continuous(limits = c(min_all, max_all))+
    scale_y_continuous(limits = c(min_all, max_all))
}


# Just print the sample names on the diag
diag_fn <- function (data, mapping) {
  label <- gsub('`', '', mapping_string(mapping$x))
  p <- ggally_text(label = label, color=colors_subphenotype[2], size=8) + theme_void()
  p
}

temp_normalized_pools_df = mydat %>% select(starts_with('Pool'))
#view(temp_normalized_pools_df)

max_all = max(temp_normalized_pools_df)
min_all = min(temp_normalized_pools_df)
p = ggpairs(temp_normalized_pools_df,
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

ggsave(filename = paste0(output_figures_dir, script_name, '-pairwise_scatter_pools.pdf'), 
       plot = p, 
       width = 10, height = 10)


###############################################################################
## pca.pdf

pc = prcomp(t(mydat))
temp = summary(pc)
plex_shapes = c(19, 17, 15, 18)

proportion_variance = 100 * round(summary(pc)$importance[2,][1:2], 2)

datt=data.frame(pheno=design$pheno,
                plex=factor(design$plex),
                PCA=pc$x[,1:2],
                sample=design$sample,
                phenotype=factor(design$phenotype, levels = c("Pool1", "Control", "PDR-L", "PDR-M", "PDR-H")))
#view(datt)
p=ggplot(datt, aes(x=PCA.PC1, y=PCA.PC2, shape=plex, color=phenotype)) + 
  geom_point(size=5) +
  scale_color_manual(values = colors_subphenotype) +
  scale_shape_manual(values = plex_shapes) +
  theme_classic() +
  theme(text=element_text(size=30),
        legend.spacing.y = unit(2, 'line'),
        legend.key.height=unit(2,"line"),
        legend.box.margin=margin(0,0,0,25)) +
  labs(x=paste0("PC1 (",proportion_variance[1],"%)"), 
       y=paste0("PC2 (",proportion_variance[2],"%)"),
       shape="TMT plex")
  
plot(p)

ggsave(filename = paste0(output_figures_dir, script_name, '-pca.pdf'), 
       plot = p, 
       width = 10, height = 6.6)


###############################################################################
## Heatmap annotated with plasma/erythrocyte: all proteins
## exp2_heatmap_all_annotated.pdf

protein = dat_gene[rownames(mydat),] %>%select("Proteins")
protein_accession <- protein %>% separate(Proteins, c(NA, "uniprot", NA),  sep="\\|", remove = TRUE, extra = 'drop')
#head(protein_accession)

plasma=readxl::read_xlsx("inputs/List_20 Plasma Proteins_Uniprot Accession.xlsx", skip = 2) %>%select("Accession (Uniprot)", "Depleted by Pierce Top12 Column (Yes/No)")
plasma$`Depleted by Pierce Top12 Column (Yes/No)`[1:14]="Yes"
plasma$`Depleted by Pierce Top12 Column (Yes/No)`[15:23]="No"
plasma = plasma %>% rename("uniprot"="Accession (Uniprot)", plasma="Depleted by Pierce Top12 Column (Yes/No)")
annotated_protein_accession = protein_accession %>% left_join(plasma)
rownames(annotated_protein_accession) = rownames(protein_accession)
annotated_protein_accession$plasma = recode(annotated_protein_accession$plasma,
                                            "Yes"= "depleted but present",
                                            "No" = "present")

erythrocyte_gene_symbols = c("ACTB", "HBG2", "CA2", "HBA1", "HBB", "HBD", "BLVRB", "BPGM", "PRDX2", "CA1", "CAT")
annotated_protein_accession$erythrocyte = rownames(annotated_protein_accession) %in% erythrocyte_gene_symbols
annotated_protein_accession$erythrocyte = ifelse(annotated_protein_accession$erythrocyte=="TRUE", "present", NA)

annotation_row=data.frame(
  plasma=as.factor(annotated_protein_accession$plasma),
  RBC=as.factor(annotated_protein_accession$erythrocyte),
  row.names=rownames(annotated_protein_accession)
)

annotation_col = data.frame(phenotype=design$phenotype, plex=design$plex)
rownames(annotation_col)=colnames(mydat)

annotation_colors=list(phenotype=c("Pool1"=colors_subphenotype[1], 
                                   "Control"=colors_subphenotype[2], 
                                   "PDR-L"=colors_subphenotype[3], 
                                   "PDR-M"=colors_subphenotype[4],
                                   "PDR-H"=colors_subphenotype[5]),
                       plasma=c("depleted but present"="purple", "present"="black"),
                       RBC=c("present"="black"),
                       plex=c("2.1" = plex_colors[1],
                              "2.2" = plex_colors[2],
                              "2.3" = plex_colors[3],
                              "2.4" = plex_colors[4])
                       )
pdf(file=paste0(output_figures_dir, script_name, '-heatmap_all_annotated.pdf'),
    width=14,
    height=10)
pheatmap::pheatmap(
  mydat,
  clustering_method="ward.D",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  border_color = NA,
  show_colnames = FALSE,
  show_rownames = FALSE,
  fontsize = 20,
  width = 14, height=10,
  cellwidth=15, cellheight=0.65)
dev.off()


###############################################################################
## Heatmap annotated with plasma/erythrocyte: Exclude plasma and erythrocyte proteins
## exp2_heatmap_exclude_erythrocyte_plasma_annotated.pdf

protein = dat_gene[rownames(mydat),] %>%select("Proteins")
protein_accession <- protein %>% separate(Proteins, c(NA, "uniprot", NA),  sep="\\|", remove = TRUE, extra = "drop")
#head(protein_accession)

plasma=readxl::read_xlsx("inputs/List_20 Plasma Proteins_Uniprot Accession.xlsx", skip = 2) %>%select("Accession (Uniprot)", "Depleted by Pierce Top12 Column (Yes/No)")
plasma$`Depleted by Pierce Top12 Column (Yes/No)`[1:14]="Yes"
plasma$`Depleted by Pierce Top12 Column (Yes/No)`[15:23]="No"
plasma = plasma %>% rename("uniprot"="Accession (Uniprot)", plasma="Depleted by Pierce Top12 Column (Yes/No)")
annotated_protein_accession = protein_accession %>% left_join(plasma)
rownames(annotated_protein_accession) = rownames(protein_accession)
annotated_protein_accession$plasma = recode(annotated_protein_accession$plasma,
                                            "Yes"= "depleted but present",
                                            "No" = "present")

erythrocyte_gene_symbols = c("ACTB", "HBG2", "CA2", "HBA1", "HBB", "HBD", "BLVRB", "BPGM", "PRDX2", "CA1", "CAT")
annotated_protein_accession$erythrocyte = rownames(annotated_protein_accession) %in% erythrocyte_gene_symbols
annotated_protein_accession$erythrocyte = ifelse(annotated_protein_accession$erythrocyte=="TRUE", "present", NA)

#head(annotated_protein_accession)

annotation_row=data.frame(
  plasma=as.factor(annotated_protein_accession$plasma),
  RBC=as.factor(annotated_protein_accession$erythrocyte),
  row.names=rownames(annotated_protein_accession)
)

annotation_col = data.frame(phenotype=design$phenotype, plex=design$plex)
rownames(annotation_col)=colnames(mydat)

annotation_colors=list(phenotype=c("Pool1"=colors_subphenotype[1], 
                                   "Control"=colors_subphenotype[2], 
                                   "PDR-L"=colors_subphenotype[3], 
                                   "PDR-M"=colors_subphenotype[4],
                                   "PDR-H"=colors_subphenotype[5]),
                       plasma=c("depleted but present"="purple", "present"="black"),
                       RBC=c("present"="black"),
                       plex=c("2.1" = plex_colors[1],
                              "2.2" = plex_colors[2],
                              "2.3" = plex_colors[3],
                              "2.4" = plex_colors[4])
                      )

subset_genes = mydat[!(rowSums(is.na(annotation_row)) != ncol(annotation_row)),]

pdf(file=paste0(output_figures_dir, script_name, '-heatmap_exclude_erythrocyte_plasma_annotated.pdf'),
    width=15, height=10)
pheatmap::pheatmap(
  subset_genes,
  clustering_method="ward.D",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 20,
  width = 15, height=10,
  cellwidth=15, cellheight=0.65)
dev.off()

###############################################################################
## Heatmap annotated with plasma/erythrocyte: Only plasma and erythrocyte proteins
## heatmap_only_erythrocyte_plasma_annotated.pdf

protein = dat_gene[rownames(mydat),] %>%select("Proteins")
protein_accession <- protein %>% separate(Proteins, c(NA, "uniprot", NA),  sep="\\|", remove = TRUE, extra = 'drop')
#head(protein_accession)

plasma=readxl::read_xlsx("inputs/List_20 Plasma Proteins_Uniprot Accession.xlsx", skip = 2) %>%select("Accession (Uniprot)", "Depleted by Pierce Top12 Column (Yes/No)")
plasma$`Depleted by Pierce Top12 Column (Yes/No)`[1:14]="Yes"
plasma$`Depleted by Pierce Top12 Column (Yes/No)`[15:23]="No"
plasma = plasma %>% rename("uniprot"="Accession (Uniprot)", plasma="Depleted by Pierce Top12 Column (Yes/No)")
annotated_protein_accession = protein_accession %>% left_join(plasma)
rownames(annotated_protein_accession) = rownames(protein_accession)
annotated_protein_accession$plasma = recode(annotated_protein_accession$plasma,
                                            "Yes"= "depleted but present",
                                            "No" = "present")

erythrocyte_gene_symbols = c("ACTB", "HBG2", "CA2", "HBA1", "HBB", "HBD", "BLVRB", "BPGM", "PRDX2", "CA1", "CAT")
annotated_protein_accession$erythrocyte = rownames(annotated_protein_accession) %in% erythrocyte_gene_symbols
annotated_protein_accession$erythrocyte = ifelse(annotated_protein_accession$erythrocyte=="TRUE", "present", NA)
#head(annotated_protein_accession)

annotation_row=data.frame(
  plasma=as.factor(annotated_protein_accession$plasma),
  RBC=as.factor(annotated_protein_accession$erythrocyte),
  row.names=rownames(annotated_protein_accession)
)

annotation_col = data.frame(phenotype=design$phenotype, plex=design$plex)
rownames(annotation_col)=colnames(mydat)

annotation_colors=list(phenotype=c("Pool1"=colors_subphenotype[1], 
                                   "Control"=colors_subphenotype[2], 
                                   "PDR-L"=colors_subphenotype[3], 
                                   "PDR-M"=colors_subphenotype[4],
                                   "PDR-H"=colors_subphenotype[5]),
                       plasma=c("depleted but present"="purple", "present"="black"),
                       RBC=c("present"="black"),
                       plex=c("2.1" = plex_colors[1],
                              "2.2" = plex_colors[2],
                              "2.3" = plex_colors[3],
                              "2.4" = plex_colors[4])
                      )

subset_genes = mydat[rowSums(is.na(annotation_row)) <2,]

pdf(file=paste0(output_figures_dir, script_name, '-heatmap_only_erythrocyte_plasma_annotated.pdf'),
    width=14, height=10)
pheatmap::pheatmap(
  subset_genes,
  clustering_method="ward.D",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 18,
  width = 14, height=10,
  cellwidth=12, cellheight=15)
dev.off()


print('Done')
###############################################################################
log_stop()