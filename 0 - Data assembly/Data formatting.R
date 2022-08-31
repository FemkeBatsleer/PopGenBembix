####Aim script####
##This script is to format population genetics dataset to usable files for analyses (structure, genepop, geneclass,...)
##
##By: Femke Batsleer

library("dplyr")
library("tidyr")
library("tidyverse")
library("adegenet")
library("adespatial")
library("spdep")
library("poppr")
library("hierfstat")
library("diveRsity")
library("graph4lg")
library("miscTools")
library("stringr")
library("RColorBrewer")
library("MASS")
library("reshape2")
library("foreign")

###Loading data and making genind objects####
## Loading data
all_data <- read.csv("../data/microsats_data.csv", sep=",")

#metadata of populations
metadatapop <- read.dbf("../data/GIS/Population_centroids_selection_collapsedyear.dbf")
metadatapopyear <- read.dbf("../data/GIS/Populations_centroids_selection_all_Lambert.dbf")

#sort(metadatapop$Population) == sort(unique(all_data$Population)) #check if population names are same in two datasets
length(unique(all_data$Sample)) == nrow(all_data) #check if all samples have a unique code
all_data <- all_data %>% unite(PopYear, c(Population, Year), sep = "", remove = FALSE) #make PopYear variable (combination of year and population)

#add metadata of pop-number, regions
all_data <- all_data %>% left_join(dplyr::select(metadatapop, Population, Number_pop, Region, Region2), by="Population")


#leave out the microsats that did not fullfill assumptions of HW, LD, Null alleles (see script 'Visualisations HW LD NA*')
all_data_assump <- all_data %>% dplyr::select(-AGBro57_a, -AGBro57_b, -AGBro35_a, -AGBro35_b,
                                       -AGBro20_a, -AGBro20_b, -AGBro16_a, -AGBro16_b,
                                       -AGBro419_a, -AGBro419_b, -AGBro111_a, -AGBro111_b,
                                       -AGBro138_a, -AGBro138_b)


####GENIND objects####
#first get own helper function from script to easily convert data to format for genind
source("data4genind_script.R")

#make dataset with loci suitable for genind
loci_data <- data4genind(dplyr::select(all_data, Sample:AGBro218_b)) #with all microsats in (to load for the HW, LD, NA testing) --> assumption testing
loci_data_assump <- data4genind(dplyr::select(all_data_assump, Sample:AGBro218_b)) #with tested microsats left out (after HW, LD and NA assumptions were tested)

#make a genind-object with adegenet package, for population level
pop_genind_noassump <- df2genind(loci_data, sep="/", ind.names=NULL, loc.names=NULL, pop=all_data$Population,
                                 strata = dplyr::select(all_data, Region, Region2, Population, Year, Number_pop), NA.char="NA", ploidy=2, type="codom")
pop_genind <- df2genind(loci_data_assump, sep="/", ind.names=NULL, loc.names=NULL, pop=all_data_assump$Population,
                        strata = dplyr::select(all_data_assump, Region, Region2, Population, Year, Number_pop), NA.char = "NA", ploidy = 2, type= "codom")
#add coordinates to the genind object
pop_genind_noassump@other$xy = as.vector(dplyr::select(all_data, "X_Lambert", "Y_Lambert"))
pop_genind@other$xy = as.vector(dplyr::select(all_data_assump, "X_Lambert", "Y_Lambert"))


#on level of population-year separately
popyear_genind_noassump <- df2genind(loci_data, sep="/", ind.names=NULL, loc.names=NULL, pop=all_data$PopYear,
                                     strata = dplyr::select(all_data, Region, Region2, PopYear, Number_pop), NA.char="NA", ploidy=2, type="codom")
popyear_genind <- df2genind(loci_data_assump, sep="/", ind.names=NULL, loc.names=NULL, pop=all_data_assump$PopYear,
                            strata = dplyr::select(all_data, Region, Region2, PopYear), NA.char = "NA", ploidy = 2, type= "codom")
popyear_genind_noassump@other$xy = as.vector(dplyr::select(all_data, "X_Lambert", "Y_Lambert"))
popyear_genind@other$xy = as.vector(dplyr::select(all_data_assump, "X_Lambert", "Y_Lambert"))


#genind object for populations where year 2018 and 2020 overlap --> AMOVA
pops_2yearsampled <- metadatapopyear %>% group_by(Population) %>% summarise(n=n()) %>% filter(n==2) #make vector which has the population in where in both years was sampled
data_2yearsampled <- all_data_assump %>% filter(Population %in% as.vector(pops_2yearsampled$Population))
nrow(all_data)
nrow(data_2yearsampled)

loci_data_2yearsampled <- data4genind(dplyr::select(data_2yearsampled, Sample:AGBro218_b))

pop_genind_2years <- df2genind(loci_data_2yearsampled, sep="/", ind.names=NULL, loc.names=NULL, pop=data_2yearsampled$Population,
                               strata = dplyr::select(data_2yearsampled, Year, Region, Region2, Population, Number_pop), NA.char="NA", ploidy=2, type="codom")
pop_genind_2years@other$xy = as.vector(dplyr::select(data_2yearsampled, "X_Lambert", "Y_Lambert"))

###Checking duplicated genotypes####
#recapture wetteren
#other recaptures from schipgatduinen were also put in as double blinds, but were removed from the data up front.

#Duplicated genotypes
mlg(pop_genind) #no duplicates present
dups <- mlg.id(pop_genind)
for (i in dups){
  if (length(dups[i]) > 1){print(i)}
}
# Create a vector of individual names without the duplicates
Nodups = indNames(pop_genind)[! indNames(pop_genind) %in% c("122", "519", "520")] #recapture from Wetteren 2018
# Create a new genind object without the duplicates
pop_genind = pop_genind[Nodups, ]
mlg(pop_genind)

pop_genind_noassump = pop_genind_noassump[indNames(pop_genind_noassump)[! indNames(pop_genind_noassump) %in% c("122", "519", "520")],]
pop_genind_2years = pop_genind_2years[indNames(pop_genind_2years)[! indNames(pop_genind_2years) %in% c("122", "519", "520")],]
popyear_genind = popyear_genind[indNames(popyear_genind)[! indNames(popyear_genind) %in% c("122", "519", "520")],]
popyear_genind_noassump = popyear_genind_noassump[indNames(popyear_genind_noassump)[! indNames(popyear_genind_noassump) %in% c("122", "519", "520")],]

mlg(pop_genind_noassump)
mlg(pop_genind_2years)
mlg(popyear_genind)
mlg(popyear_genind_noassump)

# ###saving genind-objects####
saveRDS(pop_genind, file = "../data/genind/genind_pop.RDS")
saveRDS(popyear_genind, file = "../data/genind/genind_popyear.RDS")
saveRDS(popyear_genind_noassump, file = "../data/genind/genind_popyear_noassumptested.RDS")
saveRDS(pop_genind_noassump, file = "../data/genind/genind_pop_noassumtested.RDS")
saveRDS(pop_genind_2years, file = "../data/genind/genind_pops2yearsampled.RDS")

#how to load it again, for in other and following scripts
#a <- readRDS("../data/genind_selection_pop.RDS")

###transform to GENEPOP type####
#suitable for Geneclass

# pop_info_large <- pop_info %>% group_by(LocatieNivYear) %>% summarise(n=n()) %>%
#   filter(n>=18)
# source_genind <- popsub(all_genind, sublist = as.vector(pop_info_large$LocatieNivYear))
#genind_to_genepop(pop_genind, "../data/selectionsamples_genepop.txt")
genind_to_genepop(pop_genind, "../data//genepop/selectionsamples_all.txt")
#new pop and old pop separately
#select new populations
newpops <- c("Sint-Laureins1", "Sint-Laureins2", "WarandeduinenMiddelkerke", "Veurne",
             "Vloethemveld-Zuid", "Vloethemveld-Noord", "Fort Napoleon", "Keiheuvel",
             "Averbode","Arendschot")
pop_genind_newpops <- popsub(pop_genind, sublist = newpops)
genind_to_genepop(pop_genind_newpops, "../data/genepop/selectionsamples_newpop.txt")
#select old populations
pop_genind_oldpops <- popsub(pop_genind, exclude = newpops)
genind_to_genepop(pop_genind_oldpops, "../data/genepop/selectionsamples_oldpop.txt")

# ###transform for gGENALEX####
# genind2genalex(pop_genind, filename="../data/2018westcoast_data_genalex.csv", allstrata=FALSE, overwrite=TRUE, geo=TRUE)