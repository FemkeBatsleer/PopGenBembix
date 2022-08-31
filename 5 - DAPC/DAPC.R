# ## DAPC analysis
# Aim: do multivariate DAPC analysis to assess population structure.
# A DAPC analysis can also be used to define number of clusters, we don't go that far (as this is always a matter of scale and interpretation of biological reality)
# Our main goal is to compare population genetic structure between inland and coastal population
# By: Femke Batsleer

# Good tutorials to get going:
# https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
# https://tomjenkins.netlify.app/2020/09/21/r-popgen-getting-started/#4 https://tomjenkins.netlify.app/2020/09/21/r-popgen-getting-started/#5

# Main papers explaining methods: (Jombart 2008; Jombart et al. 2010)
# 
# Package Adegenet  
# In this multivariate statistical approach, variance in the sample is partitioned into a between-group and within- group component, in an effort to maximize discrimination between groups.  
# 
# Multivariate method to define number of clusters (cfr Structure, geneland: Bayesian methods). WITHOUT assumption of panmixia.  
# data is first transformed using a principal components analysis (PCA) and subsequently clusters are identified using discriminant analysis (DA)


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
library("ggplot2")
library("ggpubr")
library("scales")

#loading data
data_genind <- readRDS("../data/genind/genind_pop.RDS")
metadatapop <- read.dbf("../data/GIS/Population_centroids_selection_collapsedyear.dbf")

data_genind <- setPop(data_genind, ~Number_pop)#set population names as numbers

#make coastal and inland genind separate
data_genind_reg <- setPop(data_genind, ~Region2) #make a genind object with region as pop-level, to be able to filter
data_genind_coast <- popsub(data_genind_reg, sublist = c("Coast")) #filter out the coastal population
data_genind_coast <- setPop(data_genind_coast, ~Number_pop) #set the population-level back to population names
data_genind_inland <- popsub(data_genind_reg, sublist = c("Inland"))
data_genind_inland <- setPop(data_genind_inland, ~Number_pop)

#DAPC analysis
#all data together
dapc.all <- dapc(data_genind, var.contrib = TRUE, scale = FALSE, n.pca = 70, n.da = nPop(data_genind) - 1)
scatter(dapc.all, cex = 2, legend = F, clabel = TRUE, posi.leg = "bottomleft", scree.pca = T, posi.pca = "topright", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

#cross-validation to test if enough PC-axisses were retained
set.seed(999)
all.x <- xvalDapc(tab(data_genind, NA.method="mean"),
                  pop(data_genind), n.rep=1000)
all.x[-1]
coastal.x <- xvalDapc(tab(data_genind_coast, NA.method="mean"),
                      pop(data_genind_coast), n.rep=1000) #first larger stretch of possibilities, now narrowed down
coastal.x[-1]
inland.x <- xvalDapc(tab(data_genind_inland, NA.method = "mean"),
                     pop(data_genind_inland))#, n.pca.max = 200, n.pca=20:70, n.rep=50)#, training.set = 0.9, result= "groupMean", scale=F, n.pca=200, n.rep=30, xval.plot=T)
inland.x[-1]


#actual dapc (without cross validation, but takes number of PC's from cross-validation)
coastal.dapc <- dapc(data_genind_coast, var.contrib = TRUE, scale = FALSE, n.pca = 50, n.da = nPop(data_genind_coast) - 1)
scatter(coastal.dapc, cex = 2, legend = F, clabel = TRUE, posi.leg = "bottomleft", scree.pca = T, posi.pca = "topright", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

inland.dapc <- dapc(data_genind_inland, var.contrib = TRUE, scale = FALSE, n.pca = 50, n.da = nPop(data_genind_coast) - 1)
scatter(inland.dapc, cex = 2, legend = F, clabel = TRUE, posi.leg = "bottomleft", scree.pca = T, posi.pca = "topright", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)


#transform to ggplot visualisation, with Tom Jenkins' tutorial https://tomjenkins.netlify.app/2020/09/21/r-popgen-getting-started/#5
ggplot_scatter_dapc <- function(genind.obj, dapc.obj, ax1=1, ax2=2, region=""){
  
  ind_coords <- as.data.frame(dapc.obj$ind.coord) %>%# Create a data.frame containing individual coordinates
    dplyr::select(1:4) #select only first 3 columns 
  colnames(ind_coords) <- c("Axis1","Axis2","Axis3", "Axis4")# select only first 4 columns and Rename columns of dataframe
  ind_coords$ax1 <- ind_coords[,ax1]#select the wanted axes
  ind_coords$ax2 <- ind_coords[,ax2]
  ind_coords$Ind <- indNames(genind.obj)# Add a column containing individuals
  ind_coords$Site <- genind.obj$pop# Add a column with the site IDs
  centroid <- aggregate(cbind(ax1, ax2) ~ Site, data = ind_coords, FUN = mean)# Calculate centroid (average) position for each population
  ind_coords <- ind_coords %>% left_join(centroid, by= "Site", suffix = c("", ".cen"))# Add centroid coordinates to ind_coords dataframe
  percent <- dapc.obj$eig/sum(dapc.obj$eig)*100
  if(region=="coastal"){cols <- c(rep("#a6611a", 3), rep("#dfc27d", 37))}
  else if(region=="inland"){cols <- c(rep("#80cdc1", 11), rep("#018571", 2))}
  else if(region=="all"){cols <- c(rep("#a6611a", 3), rep("#dfc27d", 37), rep("#80cdc1", 11), rep("#018571", 2))}
  else{cols <- colorRampPalette(brewer.pal(6, "Blues"))(nPop(genind.obj)+4)[1:nPop(genind.obj)]}
  xlab = paste("PC " ,ax1,  " (", format(round(percent[ax1], 1), nsmall=1)," %)", sep="")# Custom x and y labels
  ylab = paste("PC " ,ax2,  " (", format(round(percent[ax2], 1), nsmall=1)," %)", sep="")
  ind_coords$Site <- ordered(ind_coords$Site, levels = sort(unique(as.numeric(as.character(ind_coords$Site)))))#reorder the factor of Site, so colours are sequential with numbering of pops
  
  # Scatter plot axis 1 vs. 2
  scat.plot <- ggplot(data = ind_coords, aes(x = ax1, y = ax2))+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    geom_segment(aes(xend = ax1.cen, yend = ax2.cen, colour = Site), show.legend = FALSE)+  # spider segments
    geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+  # points
    geom_label(data = centroid, aes(label = Site, fill = Site), size = 6, show.legend = FALSE)+  # centroids
    scale_fill_manual(values = cols)+# colouring
    scale_colour_manual(values = cols)+
    scale_y_continuous(breaks= pretty_breaks())+
    labs(x = xlab, y = ylab)+  # custom labels
    theme_bw()+  # custom theme
    theme(text = element_text(size = 24))
  return(scat.plot)
}

#coast
ggplot_scatter_dapc(data_genind_coast, coastal.dapc, region="coastal") -> coastal.scatter
#ggsave("coastal_scatter_dapc.png", width=12, height=8, dpi=600) #save scatterplot automatically
#inland
ggplot_scatter_dapc(data_genind_inland, inland.dapc, region="inland") -> inland.scatter
#ggsave("inland_scatter_dapc.png", width=12, height=8, dpi=600)
#all
ggplot_scatter_dapc(data_genind, dapc.all, region="all") -> PC1PC2
#ggsave("all_scatter_dapc_PC1PC2.png", width=12, height=8, dpi=600)
#all axis 3 vs 4
ggplot_scatter_dapc(data_genind, dapc.all, ax1=3, ax2=4, region="all") -> PC3PC4
#ggsave("all_scatter_dapc_PC3PC4.png", width=12, height=8, dpi=600)

#get inland and coastal together in one frame
ggarrange(coastal.scatter, inland.scatter, labels=c("A", "B"), ncol=1, nrow=2, font.label = list(size = 20))
#ggsave("inlandcoast_scatter_dapc.png", width=12, height=16, dpi=600)
#get PC1 vs PC2 and PC3 vs PC4 together in one frame
ggarrange(PC1PC2, PC3PC4, labels=c("A", "B"), ncol=1, nrow=2, font.label = list(size=20))
#ggsave("all_scatter_PC1PC2_PC3_PC4.png", width=12, height=16, dpi=600)
