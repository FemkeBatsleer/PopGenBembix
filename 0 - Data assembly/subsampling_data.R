##Script to select subsampling data to check performance of analyses with (un)even sample sizes
#By: Femke Batsleer

library("dplyr")
library("tidyr")
library("tidyverse")
library("tibble")
library("foreign")
library("adegenet")
library("adespatial")
library("spdep")
library("poppr")
library("hierfstat")
library("graph4lg")
set.seed(23)

#load data
data_genind <- readRDS("../data/genind_selection_pop.RDS")
metadatapop <- read.dbf("../../GIS/Selection population genetics/Population_centroids_selection_collapsedyear.dbf")

#data with sample name, population
data_samples <- as.data.frame(data_genind@tab) %>% tibble::rownames_to_column("Sample") %>% add_column(Population = data_genind$pop)
samplepop_df <- data_samples %>% dplyr::select(Sample, Population)

#number of samples per pop currenlty
samples_summary <- data_samples %>% group_by(Population) %>% summarize(n=n())

#decided to subsample 10 samples, then there are a few pops with 7, 8 or 9 samples. Leave out the two pops with 5 samples

#assemble subsample dataset in few steps
#1. get populations with less than 11 samples, those data can go as a whole in the dataset; but delete two pops with 5 samples
populations_to10 <- samples_summary %>% dplyr::filter(n<=10 & n>5)
samplepop_all <- samplepop_df %>% dplyr::filter(Population %in% as.vector(populations_to10$Population)) #first add these to the dataset, which will be completed in next step

#2. loop over all populations with more than 10 samples and randomly select 10 samples; then add them to dataframe
populations_over10 <- samples_summary %>% dplyr::filter(n>10)
for(pop in as.vector(populations_over10$Population)){
  samplepop_sel <- samplepop_df %>% dplyr::filter(Population == pop) %>% slice_sample(n=10, replace=F)
  samplepop_all <- samplepop_all %>% bind_rows(samplepop_sel)
}
samplepop_all %>%  summarize(n=n())

subsamp_df <- samplepop_df %>% add_column(ID_n = 1:nrow(samplepop_df)) %>% dplyr::filter(Sample %in% as.vector(samplepop_all$Sample))

data_genind_subsamp <- data_genind[as.vector(subsamp_df$ID_n)]
#saveRDS(data_genind_subsamp, file = "../data/genind_genind_subsamp.RDS")

View(as.data.frame(data_genind_subsamp@tab) %>% tibble::rownames_to_column("Sample") %>% add_column(Population = data_genind_subsamp$pop) %>%
  group_by(Population) %>% summarize(n=n()))


pop_genind

#transform to data for geneclass2 assignment tests
pop_genind <- readRDS("../data/genind_selection_subsamp.RDS")

genind_to_genepop(pop_genind, "../GeneClass2/subsampling_selectionsamples_all.txt")
#select new populations
pop_genind_newpops <- popsub(pop_genind, sublist = c("Sint-Laureins1", "Sint-Laureins2", "WarandeduinenMiddelkerke", "Veurne",
                                                     "Vloethemveld-Zuid", "Vloethemveld-Noord", "Fort Napoleon", "Keiheuvel",
                                                     "Averbode","Arendschot"))
genind_to_genepop(pop_genind_newpops, "../GeneClass2/subsampling_selectionsamples_newpop.txt")
#select old populations
pop_genind_oldpops <- popsub(pop_genind, exclude = c("Sint-Laureins1", "Sint-Laureins2", "WarandeduinenMiddelkerke", "Veurne",
                                                     "Vloethemveld-Zuid", "Vloethemveld-Noord", "Fort Napoleon", "Keiheuvel",
                                                     "Averbode","Arendschot"))
genind_to_genepop(pop_genind_oldpops, "../GeneClass2/subsampling_selectionsamples_oldpop.txt")
