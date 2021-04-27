Oculomics_tomwgard_CU3-power_analysis
=====================================

Diffex analysis of MS proteomic vitreous data from PDR and MH
clinical samples. Analyses initiated by jinquma and extended
by cgates.

Collaborators:
- Co-I Tom Gardner (tomwgard@med.umich.edu)
- Co-I Jeff Sundstrom (jsundstrom@pennstatehealth.psu.edu)
- Felipe da Veiga Leprevost (felipevl@umich.edu)
- Venky Basrur (vbasrur@med.umich.edu) i
- Sarah Weber (srw5267@gmail.com)
- Chris Gates
- Jingqun Ma

This project encompasses two related experiments:
- experiment 1: demonstrating POC on a MS sample prep for vitreous and 
  examining biological variance among controls: 5 pooled controls + 5 individuals controls
- experiment 2: comparing PDR diseased patients with control patients

The data from experiment 2 were run in four TMT plexes bridged together with a 
reference channel. Felipe normalized the 4 plexes into a unified matrix of protein 
abundance.

Diffex analysis is based on Kammers extension of limma. See:
- Kammers K, Cole RN, Tiengwe C, Ruczinski I. Detecting Significant Changes in Protein Abundance. 
  EuPA Open Proteom. 2015 Jun;7:11-19. doi: 10.1016/j.euprot.2015.02.002. 
  PMID: 25821719; PMCID: PMC4373093.
- http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html


Files
-----
- inputs/ : MS results and sample metadata
- scripts/ : R scripts to process inputs
- outputs/ : figures and result sets

4/27/2021 cgates

