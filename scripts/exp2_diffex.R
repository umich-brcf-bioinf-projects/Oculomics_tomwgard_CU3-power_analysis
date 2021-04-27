# Oculomics analysis, visualizations for proteomics experiment 2
# jinquma/ cgates 
# 3/4/2021

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis')
source('scripts/common.R')
script_name = "exp2_diffex" 
log_start(paste0(output_logs_dir, script_name, ".log")) 


library(limma)
library(qvalue)
library(tidyverse)

## Kammers K, Cole RN, Tiengwe C, Ruczinski I. Detecting Significant Changes in Protein Abundance. EuPA Open Proteom. 2015 Jun;7:11-19. doi: 10.1016/j.euprot.2015.02.002. PMID: 25821719; PMCID: PMC4373093.
## http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
## http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r
eb.fit <- function(dat, design){ 
  n <- dim(dat)[1]
  fit <- lmFit(dat, design) 
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coefficients[, 2] 
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2] 
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual) 
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb <- results.eb[order(results.eb$p.mod), ]
  return(results.eb) 
}

dat_gene<-read.table("inputs/Report_ratio_groupby=gene_protNorm=MD_gu=0.tsv", stringsAsFactors = F, header = T, row.names =1)
pheno=read.csv("inputs/design.csv")

batch="B"
myplex=as.character(pheno$sample[pheno$batch==batch])
mydat=dat_gene[names(dat_gene) %in% myplex]
mydat =mydat %>%
  select(-contains("Pool")) %>%
  drop_na()

design=pheno[pheno$sample %in% colnames(mydat),]
design=design[match(colnames(mydat), design$sample ),]
tr<-as.factor(as.character(design$pheno))
design=model.matrix(~tr)

res.eb <- eb.fit(mydat, design)
temp=dat_gene[match(rownames(res.eb), rownames(dat_gene)),]
res.eb$proteins=temp$Proteins
res.eb$Gene=rownames(temp)
res.eb1 <- res.eb %>% separate(proteins, c(NA, "protein", NA),  remove = TRUE, extra = "drop")
sig=subset(res.eb1, p.mod<=0.05)

write.table(res.eb1,
            file = paste0(output_data_dir, script_name, "-kammers.tsv"),
            sep = "\t",
            row.names = FALSE)

print('Done')
###############################################################################
log_stop()
