#Aim: Calculations Null alleles, LD, HW
#Author: Femke Batsleer
#Date: February 2022
##WARNING: some commands take a very long time once you start running them!

library("dplyr")
library("tidyr")
library("tidyverse")
library("purrr")
library("adegenet")#v2.0.0 install.packages("adegenet", version = "2.0.0")
library("poppr")
library("ggplot2")
#library("genepop")
library("graph4lg")
library("pegas")
library("hierfstat")
library("PopGenReport")#v2.2.1 has to be installed
library("stringr")
library("devtools")
library("pkgload")


###Loading data and converting to genind-object for caclulations####
#get own function to easily convert data to format for genind
source("data4genind_script.R")


all_genind_pop <- readRDS("../data/genind/genind_pop_noassumptested.RDS")
all_genind_popyear <- readRDS("../data/genind/genind_popyear_noassumptested.RDS")

#subset of population-year with > 10 samples
pop_subsetmt10samp <- all_data %>% group_by(Population) %>% summarise(n=n()) %>% filter(n >= 10) %>% dplyr::select(Population)
all_genind_pop_mt10samp <- popsub(all_genind_pop, sublist = as.vector(pop_subsetmt10samp$Population))

popyear_subsetmt10samp <- all_data %>% group_by(PopYear) %>% summarise(n=n()) %>% filter(n >= 10) %>% dplyr::select(PopYear)
all_genind_popyear_mt10samp <- popsub(all_genind_popyear, sublist = as.vector(popyear_subsetmt10samp$PopYear))


### Null alleles ####
##WARNING: this function takes quite long to be run!
#With PopGenReport

#If the 95% confidence interval includes zero, it indicates that the frequency of null alleles at a locus
#does not significantly differ from zero.

##Function to calculate null alleles for all pops in genind-object
null.all_test <- function(genind_input){
  #initialize dataframe
  null.all_df <- data.frame(locus=character(), observed=numeric(), median=numeric(),
                            percentile2.5th= numeric(), percentile97.5th=numeric(),
                            pop = character())
  # popname = "Geel-Bel2018"
  # genind_input = all_genind_reg
  for (popname in popNames(genind_input)){
    print(popname)
    pop_sub <- popsub(genind_input, sublist=c(popname)) #take subpopulation per subpopulation
    #run HW test for one population
    print(paste0("Calculating null allele test for ", as.character(match(popname, popNames(genind_input))),
                 " out of ", length(popNames(genind_input)), " populations. Since ", Sys.time()))
    #reload package again?
    #reload(pkg = pkgload::inst("PopGenReport"), quiet=T)
    set.seed(5)
    Nulltest <- null.all(pop_sub) #do the HW-test for this population, gives output per locus
    #bind in dataframe
    print("beginning of dataframe")
    Nulltest_df <- as.data.frame(Nulltest$null.allele.freq$summary2) %>%
      tibble::rownames_to_column(var="locus") %>%
      pivot_longer(-locus) %>%
      pivot_wider(names_from=locus, values_from=value) %>%
      set_names(c("locus", "observed", "median", "percentile2.5th", "percentile97.5th")) %>%
      select(c("locus","observed", "median", "percentile2.5th", "percentile97.5th")) %>%
      mutate(pop = popname)
    null.all_df <- null.all_df %>% bind_rows(Nulltest_df)
  }
  return(null.all_df)
}


# #leave out the ones that give errors...
##Populations level
#WARNING: runs for quite some time...
Nulltest_pop <- null.all_test(popsub(all_genind_pop, exclude = c("Raversijde3", "Averbode2", "Sint-Laureins1", "Raversijde4")))
write.csv(Nulltest_pop, "Outputs/Output Null Alleles population level.csv")
##PopYear level
Nulltest_popyear <- null.all_test(popsub(all_genind_popyear, exclude = c("Raversijde32020", "HB-Oost2020", "Astridlaan2020", "Oosthoekduinen2020",
                                                                         "Westhoek vissersdorp2020", "Schipgatduinen2020", "Raversijde42020",
                                                                         "Ter Yde West2020", "Sint-Laureins12020", "Westhoek NO2020", "Averbode22020")))
write.csv(Nulltest_popyear, "Outputs/Output Null Alleles popyear level.csv")
##subset mt10samp level
Nulltest_pop_mt10samp <- null.all_test(all_genind_pop_mt10samp) #works, no errors
write.csv(Nulltest_pop_mt10samp, "Outputs/Output Null Alleles population morethan10samp level.csv")

Nulltest_popyear_mt10samp <- null.all_test(popsub(all_genind_popyear_mt10samp, exclude= c("Oosthoekduinen2020")))
write.csv(Nulltest_popyear_mt10samp, "Outputs/Output Null Alleles popyear morethan10samp level.csv")

### Linkage Disequilibrium ####

#With poppr

##Function to calculate LD for all pops in genind-object
LD.all_test <- function(genind_input){
  #initialize dataframe
  LD.all_df <- data.frame(pairloci=character(), Ia=numeric(),
                          p.Ia=numeric(), pop = character())
  
  for (popname in popNames(genind_input)){
    # genind_input <- all_genind_reg
    # popname <- "Geel-Bel2018"
    pop_sub <- popsub(genind_input, sublist=c(popname)) #take subpopulation per subpopulation
    #run HW test for one population
    print(paste0("Calculating LD-test for ", as.character(match(popname, popNames(genind_input))),
                 " out of ", length(popNames(genind_input)), " populations. Since ", Sys.time()))
    LDtest <- pair.ia(pop_sub, sample=99, quiet=T, plot=F) #do the LD-test for this population, gives output per locus
    #bind in dataframe
    LDtest_df <- as.data.frame(LDtest) %>% tibble::rownames_to_column("pairloci") %>%
      select(pairloci, Ia, p.Ia) %>%
      set_names(c("pairloci", "Ia", "p.Ia")) %>%
      mutate(pop = popname)
    LD.all_df <- LD.all_df %>% bind_rows(LDtest_df)
  }
  return(LD.all_df)
}

##Population level
LD.all_pop <- LD.all_test(genind_input = all_genind_pop)
write.csv(LD.all_pop, "Outputs/Output LD population level.csv")
##PopYear level
LD.all_popyear <- LD.all_test(genind_input = all_genind_popyear)
write.csv(LD.all_popyear, "Outputs/Output LD popyear level.csv")
##subset mt10samp level
LD.all_pop_mt10samp <- LD.all_test(genind_input = all_genind_pop_mt10samp)
write.csv(LD.all_pop_mt10samp, "Outputs/Output LD population morethan10samp level.csv")
LD.all_popyear_mt10samp <- LD.all_test(genind_input = all_genind_popyear_mt10samp)
write.csv(LD.all_popyear_mt10samp, "Outputs/Output LD popyear morethan10samp level.csv")


### Hardy-Weinberg ####

#With pegas
HW.all_test <- function(genind_input){
  #initialize dataframe
  HW.all_df <- data.frame(locus=character(), chi2=numeric(),
                          df=numeric(), Pr.chi2=numeric(),
                          Pr.exact=numeric(), pop = character())
  
  for (popname in popNames(genind_input)){
    pop_sub <- popsub(genind_input, sublist=c(popname)) #take subpopulation per subpopulation
    #run HW test for one population
    print(paste0("Calculating HW-test for ", as.character(match(popname, popNames(genind_input))),
                 " out of ", length(popNames(genind_input)), " populations. Since ", Sys.time()))
    HWtest <- hw.test(pop_sub, B=1000) #do the HW-test for this population, gives output per locus
    #bind in dataframe
    HWtest_df <- as.data.frame(HWtest) %>% tibble::rownames_to_column("locus") %>%
      set_names(c("locus", "chi2", "df", "Pr.chi2", "Pr.exact")) %>%
      mutate(pop = popname)
    HW.all_df <- HW.all_df %>% bind_rows(HWtest_df)
  }
  return(HW.all_df)
}
##Function to calculate HW for all pops in genind-object


##Population level
HW.all_pop <- HW.all_test(genind_input = all_genind_pop)
write.csv(HW.all_pop, "Outputs/Output HW population level.csv")
##PopYear level
HW.all_popyear <- HW.all_test(genind_input = all_genind_popyear)
write.csv(HW.all_popyear, "Outputs/Output HW popyear level.csv")
##subset mt10samp level
HW.all_pop_mt10samp <- HW.all_test(genind_input = all_genind_pop_mt10samp)
write.csv(HW.all_pop_mt10samp, "Outputs/Output HW population morethan10samp level.csv")
HW.all_popyear_mt10samp <- HW.all_test(genind_input = all_genind_popyear_mt10samp)
write.csv(HW.all_popyear_mt10samp, "Outputs/Output HW popyear morethan10samp level.csv")