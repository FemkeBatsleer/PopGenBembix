# PopGenBembix
Data and code for population genetics Bembix rostrata
By: Femke Batsleer

## Data and scripts explanation
- **0 - Data assembly**: main script for data assembly and formatting ('Data formatting.R') and a helper script ('data4genind_script.R'). Transforms the main data ('data>microsats_data.csv') into formats and files usable in all following analysis, including metadata from the populations.
- **1 - HW, LD and NA assumptions**: script to calculate Hardy-Weinberg (HW), Linkage Disequilibrium (LD) and Null Alleles (NA) for all loci ('Calculations HW LD NA.R'). This scripts creates output ('Outputs/') that is read in 'Visualisations HW LD NA**.Rmd' to create tables and figures to summarize the assumption-tests. There are two parallel scripts, one where the assumptions were tested only with populations that had more than 10 samples ('*-morethan10samplespops') and one for all. Outputs are only given for the first, but can be created for all pops and used in the script 'Visualisations HW LD NA all populations.Rmd' as check. The file 'Probabilities of number of significant tests for X tests performed.xlsx' gives the calculations of thresholds for multiple testing according to the example in Waples (2015) Journal of Heredity.
- **2 - AMOVA**: script ('AMOVA.Rmd') to perform hierarchical AMOVA
- **3 - Pop stats and gen dist**: script ('Pop stats and genetic distances.Rmd') to calcualte population level statistics and pairwise genetic distances and differentiations. The folder 'Outputs/' holds output created in the script with genetic differentiation values. The output 'Neisgeneticdistance_hierfstat.csv' is used in IBD-analysis below as input.
- **4 - IBD**: script to calculate Isolation-By-Distance. Input 'Neisgeneticdistance_hierfstat.csv' from folder above needed as input. Permutation tests, RDA and Mantel tests are performed to assess IBD statistically.
- **5 - DAPC**: script ('DAPC.R') to perform multivariate DAPC analysis to assess population structure.
- **6 - Assignment tests**: detailed explanation how Geneclass was run for the analysis ('Details for running Geneclass.txt'). Output from geneclass is 'Results_indassign_newpopsoldpops_RM_probRM10000MC.csv'. Script to summarize the export file from Geneclass ('tables_from_output_Geneclass.R'), which also creates output ('Outputs for QGIS visualisations/') to make a flow chart in QGIS.
- **data**: all data used in previous scripts. 'microsats_data.csv' is the main raw datafile. 'genind/' holds genind-objects (*.RDS) used in most R-scripts as the main data format, the R-package adegenet class for individual genotypes. 'genepop/' holds data to be used in the assignment tests (genepop-format). 'GIS/' holds the shapefiles with metadata of samples and populations.


## Packages used
- **genetic analysis packages**: adegenet, diveRsity, graph4lg, hierfstat, pegas, PopGenReport, poppr
- **data handling and manipulation**: dplyr, foreign, reshape2, stringr, tibble, tidyr, tidyverse
- **multivariate analysis**: vegan, adespatial
- **visualisations**: Ggally, ggplot2, ggpubr, ggsci, RColorBrewer, scales
- **other**: car, devtools, MASS, miscTools, pkgload, purr, spdep
