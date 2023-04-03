#FIS plot combining full dataset and subsampled dataset
library(dplyr)
library(ggplot2)

subsamp_df <- read.csv(file="Fis_subsampling.csv") %>% mutate(dataset = "subsampled")
full_df <- read.csv(file="Fis_all.csv") %>% mutate(dataset = "full")

FIS_df <- subsamp_df %>% bind_rows(full_df)

plot_Fis <- ggplot(FIS_df, aes(x=reorder(as.character(Number_pop), Number_pop), y=Fis_mean, col=dataset)) +
  geom_point(position= position_dodge(width=0.8), na.rm=TRUE, size=2) +
  geom_errorbar(aes(ymin=Fis_mean-Fis_se, ymax=Fis_mean+Fis_se), position=position_dodge(width=0.8), na.rm=T, linewidth=1) +
  xlab("Population") + ylab(expression(F["IS"])) +
  theme_bw() +
  #theme(axis.text.x  = element_text(angle = 90, vjust=0.5, hjust=1)) +
  geom_hline(yintercept=0) + geom_hline(yintercept=mean(Fis_pop$Fis_mean), linetype="dashed")
plot_Fis
