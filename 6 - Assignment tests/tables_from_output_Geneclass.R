###Aim: making summary tables of assignment tests and create output to load into QGIS to make flow charts
##By: Femke Batsleer

library(dplyr)
library(tidyr)
library(stringr)
library(foreign)

###Loading and formatting data####
#load metadata
metadatapop <- read.dbf("../data/GIS/Population_centroids_selection_collapsedyear.dbf") %>%
  mutate(Population = gsub("\\ ", "", Population))

#load results from assignment tests in Geneclass2
RM_P <- read.csv2(file = "Results_indassign_newpopsoldpops_RM_probRM10000MC_forR.csv", header=T) %>%
  rename("AssignedSample"="Assigned.sample") %>%
  dplyr::select(-Used.loci, -Missing.loci) %>%
  separate(col=AssignedSample, into=c("pop", "Sample"), sep="_") %>%
  setNames(gsub("_.*","",names(.))) %>% #delete trailing part of name
  setNames(gsub("\\.", "-", names(.))) %>% #change . into -
  mutate(across(.cols=3:(ncol(.)-1), .fns=as.numeric))
   

#add the Pmax and which pop this is (pop_max)
RM_P_f <- RM_P %>% dplyr::select(pop, Sample) %>% mutate(id=1:n()) %>%
  left_join(RM_P %>% select(3:(ncol(.)-1)) %>%
              mutate(id=1:n()) %>%
              pivot_longer(-id) %>%
              #pivot_wider(names_from=name,values_from=value) %>%
              group_by(id) %>%
              filter(value==max(value)) %>%
              slice(1) %>% #there is one with duplicated value, but both from within coast and low prob
              rename(pop_max=name, Pmax=value), by="id") %>%
  #add Psampled
  left_join(RM_P %>% dplyr::select(pop) %>% mutate(id=1:n()) %>%
              right_join(RM_P %>% select(3:(ncol(.)-1)) %>%
                           mutate(id=1:n()) %>%
                           pivot_longer(-id), by=("id")) %>%
              group_by(id) %>% dplyr::filter(pop==name) %>% dplyr::select(-pop, -name) %>%
              rename(Psampled=value), by="id") %>%
  mutate(pop_max = gsub("\\.", "-", pop_max)) %>%
  left_join(metadatapop %>% dplyr::select(Population, Number_pop, Region) %>%
              rename("Number_pop_sampled"="Number_pop", "Region_sampled"="Region")
            , by=(c("pop"="Population"))) %>%
  left_join(metadatapop %>% dplyr::select(Population, Number_pop, Region) %>%
              rename("Number_pop_max"="Number_pop", "Region_max"="Region")
            , by=(c("pop_max"="Population"))) %>%
  mutate(Lambda=Pmax/Psampled, lnLambda = log(Lambda)) #see Rannala & Mountain for Lambda explanation

#Write a file with all assigned individuals, including the probabilities
# write.csv(RM_P_f, file="new and old pops/list_assigned_all.csv")

#table with connections
connections_pops_df_RM <- RM_P_f %>% group_by(pop, pop_max) %>%
  summarise(n = n())

connections_regions_df_RM <- RM_P_f %>% group_by(Region_sampled, Region_max) %>%
  summarise(n = n())

#make long format out of connections/edges and add info from metadata
edges_pops_RM <- connections_pops_df_RM %>% rename("to"="pop", "from"="pop_max", "intensity"="n") %>%
  left_join(metadatapop %>% dplyr::select(Population, Number_pop, Region) %>% rename("to"="Population", "idto"="Number_pop", "regionto"="Region"),#add id's of to pops
            by="to") %>%
  left_join(metadatapop %>% dplyr::select(Population, Number_pop, Region) %>% rename("from"="Population", "idfrom"="Number_pop", "regionfrom"="Region"), #add id's of from population
            by="from")

#define new pops vector
newpops <- c("Sint-Laureins1", "Sint-Laureins2", "WarandeduinenMiddelkerke", "Veurne",
             "Vloethemveld-Zuid", "Vloethemveld-Noord", "FortNapoleon", "Keiheuvel",
             "Averbode","Arendschot")

###Summary tables####

#new populations are filtered out
#summarise per region, how many within and between
region_summarised_RM <- RM_P_f %>% dplyr::filter(!pop %in% newpops) %>%
  group_by(Region_sampled, Region_max) %>%
  summarise(n=n(), meanPmax = mean(Pmax), sdPmax = sd(Pmax))
region_summarised_RM

#summarise how many in which population/region were assigned to sampled region
sampled_summarised_RM <- RM_P_f %>%  dplyr::filter(!pop %in% newpops) %>%
  mutate(home = ifelse(pop==pop_max, 1, 0)) %>%
  mutate(home_region = ifelse(Region_sampled==Region_max, 1, 0)) %>%
  group_by(Region_sampled) %>% summarise(n_home = sum(home), n_all=n(), n_homeregion=sum(home_region)) %>%
  mutate(ratio_pop = n_home/n_all, ratio_region = n_homeregion/n_all)
sampled_summarised_RM

#new populations
sampled_summarised_RM_new <- RM_P_f %>%  dplyr::filter(pop %in% newpops) %>%
  mutate(home = ifelse(pop==pop_max, 1, 0)) %>%
  mutate(home_region = ifelse(Region_sampled==Region_max, 1, 0)) %>%
  group_by(Region_sampled) %>% summarise(n_home = sum(home), n_all=n(), n_homeregion=sum(home_region)) %>%
  mutate(ratio_pop = n_home/n_all, ratio_region = n_homeregion/n_all)
sampled_summarised_RM_new

#summarised per pop
pops_summarised_RM <- RM_P_f %>% dplyr::filter(!pop %in% newpops) %>%
  mutate(home = ifelse(pop==pop_max, 1, 0)) %>%
  mutate(home_region = ifelse(Region_sampled==Region_max, 1, 0)) %>%
  group_by(pop) %>% summarise(n_home = sum(home), n_all=n(), n_homeregion=sum(home_region)) %>%
  mutate(ratio = n_home/n_all)
pops_summarised_RM

###Nodes and edges for in QGIS####
#write separate csv's for 'new' and 'old' populations
#write csv with only nodes of new pops for use in QGIS
write.csv(metadatapop %>% dplyr::filter(Population %in% newpops), file="Outputs for QGIS visualisations/nodes_pops_new.csv")
write.csv(metadatapop %>% dplyr::filter(!Population %in% newpops), file="Outputs for QGIS visualisations/nodes_pops_old.csv")
#write csv files with edges for use in QGIS
#old pops
#within coast
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-coast"|regionto=="France-Picardie") & (regionfrom=="Flanders-coast"|regionfrom=="France-Picardie") ) %>%
            dplyr::filter(!to %in% newpops),
          file="Outputs for QGIS visualisations/edges_pops_RM_withincoast_old.csv")
#within inland
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-inland"|regionto=="Wallonia") & (regionfrom=="Flanders-inland"|regionfrom=="Wallonia") ) %>%
            dplyr::filter(!to %in% newpops),
          file="Outputs for QGIS visualisations/edges_pops_RM_withininland_old.csv")
#from coast to inland
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-inland"|regionto=="Wallonia") & (regionfrom=="Flanders-coast"|regionfrom=="France-Picardie") ) %>%
            dplyr::filter(!to %in% newpops),
          file="Outputs for QGIS visualisations/edges_pops_RM_fromcoasttoinland_old.csv")
#from inland to coast
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-coast"|regionto=="France-Picardie") & (regionfrom=="Flanders-inland"|regionfrom=="Wallonia") ) %>%
            dplyr::filter(!to %in% newpops),
          file="Outputs for QGIS visualisations/edges_pops_RM_frominlandtocoast_old.csv")

#new pops
#within coast
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-coast"|regionto=="France-Picardie") & (regionfrom=="Flanders-coast"|regionfrom=="France-Picardie") ) %>%
            dplyr::filter(to %in% newpops),
          file="Outputs for QGIS visualisations/edges_pops_RM_withincoast_new.csv")
#within inland
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-inland"|regionto=="Wallonia") & (regionfrom=="Flanders-inland"|regionfrom=="Wallonia") ) %>%
            dplyr::filter(to %in% newpops),
          file="Outputs for QGIS visualisations/edges_pops_RM_withininland_new.csv")
#from coast to inland
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-inland"|regionto=="Wallonia") & (regionfrom=="Flanders-coast"|regionfrom=="France-Picardie") ) %>%
            dplyr::filter(to %in% newpops),
          file="Outputs for QGIS visualisations/edges_pops_RM_fromcoasttoinland_new.csv")
#from inland to coast
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-coast"|regionto=="France-Picardie") & (regionfrom=="Flanders-inland"|regionfrom=="Wallonia") ) %>%
            dplyr::filter(to %in% newpops),
          file="Outputs for QGIS visualisations/edges_pops_RM_frominlandtocoast_new.csv")



