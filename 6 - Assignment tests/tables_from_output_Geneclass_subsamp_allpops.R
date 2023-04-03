library(dplyr)
library(tidyr)
library(stringr)
library(foreign)
library(ggplot2)

#load metadata
metadatapop <- read.dbf("../data/GIS/Population_centroids_selection_collapsedyear.dbf") %>%
  mutate(Population = gsub("\\ ", "", Population))

####RM method####

RM_P <- read.csv2(file = "Results_indassign_subsamp_allpops_RM_probRM10000MCalpha001_forR.csv", header=T) %>%
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

#Write a file with all assigned individuals, including the probabilities and lambda calculations
write.csv(RM_P_f, file="Outputs for QGIS visualisations/subsamp all pops/list_assigned_subsamp.csv")

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
write.csv(edges_pops_RM, file="Outputs for QGIS visualisations/subsamp all pops/edges_pops_RM_subsamp_all.csv") #write to csv file

# #write a file for the nodes, from metadatapop
write.csv(metadatapop, file="Outputs for QGIS visualisations/subsamp all pops/nodes_pops_subsamp_all.csv")


#write separate csv's for within and among coastal or inland regions
write.csv(edges_pops_RM %>%
            dplyr::filter((regionfrom=="Flanders-coast"|regionfrom=="France-Picardie") ),
          file="Outputs for QGIS visualisations/subsamp all pops/edges_pops_RM_fromcoast_subsamp_all.csv")
write.csv(edges_pops_RM %>%
            dplyr::filter((regionfrom=="Flanders-inland"|regionfrom=="Wallonia") ),
          file="Outputs for QGIS visualisations/subsamp all pops/edges_pops_RM_frominland_subsamp_all.csv")

#write csv files with edges for use in QGIS
#within coast
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-coast"|regionto=="France-Picardie") &
                                            (regionfrom=="Flanders-coast"|regionfrom=="France-Picardie") ),
          file="Outputs for QGIS visualisations/subsamp all pops/edges_pops_RM_withincoast_subsamp.csv")
#within inland
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-inland"|regionto=="Wallonia") &
                                            (regionfrom=="Flanders-inland"|regionfrom=="Wallonia") ),
          file="Outputs for QGIS visualisations/subsamp all pops/edges_pops_RM_withininland_subsamp.csv")
#from coast to inland
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-inland"|regionto=="Wallonia") &
                                            (regionfrom=="Flanders-coast"|regionfrom=="France-Picardie") ),
          file="Outputs for QGIS visualisations/subsamp all pops/edges_pops_RM_fromcoasttoinland_subsamp.csv")
#from inland to coast
write.csv(edges_pops_RM %>% dplyr::filter((regionto=="Flanders-coast"|regionto=="France-Picardie") &
                                            (regionfrom=="Flanders-inland"|regionfrom=="Wallonia") ),
          file="Outputs for QGIS visualisations/subsamp all pops/edges_pops_RM_frominlandtocoast_subsamp.csv")


#summarise per region, how many within and between
region_summarised_RM <- RM_P_f %>%
  group_by(Region_sampled, Region_max) %>%
  summarise(n=n(), meanPmax = mean(Pmax), sdPmax = sd(Pmax))
region_summarised_RM

#summarise how many in which population/region were assigned to sampled region
sampled_summarised_RM <- RM_P_f %>%
  mutate(home = ifelse(pop==pop_max, 1, 0)) %>%
  mutate(home_region = ifelse(Region_sampled==Region_max, 1, 0)) %>%
  group_by(Region_sampled) %>% summarise(n_home = sum(home), n_all=n(), n_homeregion=sum(home_region)) %>%
  mutate(ratio_pop = n_home/n_all, ratio_region = n_homeregion/n_all)
sampled_summarised_RM

# #check how many that are assigned to other region, are from either coast or inland
# sampled_other_detail_RM <- RM_P_f %>%
#   mutate(home = ifelse(pop==pop_max, 1, 0)) %>%
#   mutate(home_region = ifelse(Region_sampled==Region_max, 1, 0)) %>%
#   mutate(other_region = ifelse(Region_sampled!=Region_max, 1, 0)) %>%
#   mutate(other_coast = ifelse(Region_sampled!=Region_max &
#                                 Region_max %in% c("Flanders-coast", "France-Picardie"), 1, 0)) %>%
#   mutate(other_inland = ifelse(Region_sampled!=Region_max &
#                                  Region_max %in% c("Flanders-inland", "Wallonie"), 1, 0)) %>%
#   group_by(Region_sampled) %>% summarise(n_other_coast = sum(other_coast), n_other_inland = sum(other_inland),
#                                          n_other = sum(other_region), n_all=n()) %>%
#   mutate(ratio_other_coast = n_other_coast/n_all, ratio_other_inland = n_other_inland/n_all,
#          ratio_other = n_other/n_all)
# sampled_other_detail_RM


#summarised per pop
pops_summarised_RM <- RM_P_f %>%
  mutate(home = ifelse(pop==pop_max, 1, 0)) %>%
  mutate(home_region = ifelse(Region_sampled==Region_max, 1, 0)) %>%
  group_by(pop) %>% summarise(n_home = sum(home), n_all=n(), n_homeregion=sum(home_region)) %>%
  mutate(ratio = n_home/n_all)
pops_summarised_RM



#plot with barplots instead of table
#make dataset
sampled_barplot_df <- RM_P_f %>%
  mutate(home = ifelse(pop==pop_max, 1, 0)) %>%
  mutate(Region_assigned = ifelse(home==1, "Original population", as.character(Region_max))) %>%
  group_by(Region_sampled, Region_assigned) %>%
  summarise(n_links = n()) %>%
  mutate(Region_assigned = recode_factor(Region_assigned, "France-Picardie"="Coastal France", "Flanders-coast"="Coastal Flanders",
                                         "Flanders-inland"="Inland Flanders", "Original population"="Original population"),
         Region_sampled = recode_factor(Region_sampled, "France-Picardie"="Coastal France", "Flanders-coast"="Coastal Flanders",
                                        "Flanders-inland"="Inland Flanders", "Wallonia"="Inland Wallonia"))%>%
  mutate(Region_assigned = factor(Region_assigned, levels=c("Coastal France", "Coastal Flanders","Inland Flanders",  "Original population")),
         Region_sampled = factor(Region_sampled, levels=c("Coastal France", "Coastal Flanders", "Inland Flanders", "Inland Wallonia")))

#labels of total numbers of individuals
labels_ntot <- sampled_barplot_df %>% group_by(Region_sampled) %>% summarise(n_tot = sum(n_links)) %>%
  mutate(n_text = "n = ", n_tot = as.character(n_tot)) %>% mutate(n_tot_text = paste0(n_text, n_tot))

#make plot
assign_plot <- ggplot(sampled_barplot_df, aes(y= n_links, x=Region_sampled, fill=Region_assigned)) +
  geom_bar(position="fill", stat="identity")+
  #geom_text(stat = "count", aes(label = ..count.., y=1),position=position_stack(0.5))+
  geom_text(data = labels_ntot, 
            aes(x = Region_sampled, label =  n_tot_text, y = 1.025, fill=NULL), size = 4)+
  guides(fill=guide_legend(title="Sample assigned to")) +
  scale_fill_manual(values=c("#a6611a", "#dfc27d", "#80cdc1", "#bababa")) +
  labs(x = "Region of sampled population", y = "Ratio of assigned individuals to regions (%)")+  # custom labels
  theme_bw()+  # custom theme
  theme(text = element_text(size = 14), axis.text.x = element_text(face="bold"),
        panel.grid.major.x = element_blank()) +
  coord_cartesian(clip='off')
assign_plot
