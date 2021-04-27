# Oculomics analysis, visualizations for proteomics experiment 2
# jinquma/ cgates 
# 3/23/2021

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Oculomics_tomwgard_CU3-power_analysis_cgates_jingquma/cgates_analysis')
source('scripts/common.R')
script_name = "exp2_figures_particle" 
log_start(paste0(output_logs_dir, script_name, ".log")) 

library(ggbeeswarm)
library(tidyverse)


particle_df <-as.tibble(read.table("inputs/particle.tsv", stringsAsFactors = F, header = T)) #, row.names =1))
head(particle_df)

particle_df$subphenotype_f = factor(particle_df$Phenotype,
                                 levels=c('Control', 'PDR-L', 'PDR-M', 'PDR-H')
                                 )
particle_df = particle_df %>% mutate(phenotype =  
                       case_when(
                      subphenotype_f == 'Control' ~ 'Control',
                      TRUE ~ 'PDR'))
particle_df$phenotype_f = factor(particle_df$phenotype)

colors_subphenotype_f = colors_subphenotype[2:5]

# particle_df$percentage_f = factor(format(floor(particle_df$Percentage / 0.05) * 0.05, nsmall = 2))
# sizes_percentage_f = as.numeric(levels(particle_df$percentage_f))*20 * 3
# 
# particle_df$diameter_f = factor(format(floor(particle_df$Diameter / 50) * 50))
# sizes_diameter_f = 3 * as.numeric(levels(particle_df$diameter_f)) / 50
# 
# 
# p=ggplot(particle_df, 
#          aes(x=subphenotype_f, y=Abundance, color=subphenotype_f, size=Diameter)) + 
#   geom_point(alpha=0.8) +
#   scale_color_manual(values = colors_subphenotype_f) +
#   #scale_size_manual(values = sizes_percentage_f) +
#   theme_classic() +
#   theme(text = element_text(size = 30)) +
#   # labs(color = 'phenotype', 
#   #      size = 'percentage',
#   #      x='diameter (nm)',
#   #      y='abundance') + 
#   guides(color = guide_legend(override.aes = list(size=10)))
# 
# plot(p)
# 
# ggsave(filename = paste0(output_figures_dir, script_name, '-abundance_x_diameter_scatter.pdf'), 
#        plot = p, 
#        width = 10, height = 6.6)
# 

big_particle_df = particle_df[particle_df$Percentage>0.10,] # & particle_df$Abundance<(5 * 10^9),]

p=ggplot(big_particle_df, 
         aes(y=Abundance, x=phenotype_f)) + 
  geom_boxplot() +
  geom_beeswarm(aes(color=subphenotype_f), alpha=0.7, size=6, dodge.width = 0.3) +
  scale_color_manual(values = colors_subphenotype_f) +
  #scale_size_manual(values = sizes_diameter_f) +
  theme_classic() +
  theme(text = element_text(size = 30)) +
  labs(color = 'subphenotype',
       x='phenotype',
       y='abundance') #+
#  guides(color = guide_legend(override.aes = list(size=10)))

plot(p)

ggsave(filename = paste0(output_figures_dir, script_name, '-phenotype_x_abundance.pdf'), 
       plot = p, 
       width = 10, height = 6.6)


p=ggplot(big_particle_df, 
         aes(y=Abundance, x=subphenotype_f)) + 
  geom_boxplot() +
  geom_beeswarm(alpha=0.7, size=5) +
  scale_color_manual(values = colors_subphenotype_f) +
  #scale_size_manual(values = sizes_diameter_f) +
  theme_classic() +
  theme(text = element_text(size = 30)) +
  labs(color = 'subphenotype',
       x='phenotype',
       y='abundance') #+
#  guides(color = guide_legend(override.aes = list(size=10)))

plot(p)

ggsave(filename = paste0(output_figures_dir, script_name, '-subphenotype_x_abundance.pdf'), 
       plot = p, 
       width = 10, height = 6.6)

p=ggplot(big_particle_df, 
         aes(y=Diameter, x=phenotype_f, size=Abundance)) + 
  geom_boxplot() +
  geom_beeswarm(aes(color=subphenotype_f), alpha=0.7, dodge.width = 0.3) +
  scale_color_manual(values = colors_subphenotype_f) +
  #scale_size_manual(values = sizes_diameter_f) +
  theme_classic() +
  theme(text = element_text(size = 30))
  # labs(color = 'subphenotype',
  #      x='phenotype',
  #      y='abundance') #+
#  guides(color = guide_legend(override.aes = list(size=10)))

plot(p)


print('Done')
###############################################################################

log_stop()
sessionInfo()
