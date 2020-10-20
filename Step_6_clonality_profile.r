##Oct 2020 : Programme permettant de générer les profils de clonalité (pie charts) + indique le nombre de clones et cellules infectées ##
library(reldist)
library(ggplot2)
library(viridis)

T <- read.csv('20_SFV_co_rmart_no-double_clustered.csv', header=TRUE, flush=TRUE, sep='\t')

#Calcul du nombre total de cellules infectées et de l'index d'oligoclonalité
pielabels <- paste(c('Sample_20\nSFV co-infected baboon','\n','Total clones:', nrow(T),'\n','Total infected cells:', round(sum(T$Sisters),digits = 0),'\n', 'OCI:',round(gini(T$Sisters),digits = 2)),collapse=" ")

#Pie_chart = profil de clonalité
ggplot(T, aes(x="", y=Relabun, fill=Relabun)) + ggtitle(pielabels) +
  geom_bar(stat="identity", width=1, color='black') +
  coord_polar("y", start=0)+theme_bw()+theme(plot.title = element_text(face="bold",hjust = 0.5),
                                                axis.ticks=element_blank(),legend.position="none",
                                                     axis.text.y=element_blank(),
                                                      axis.text.x=element_text(color="white"),
                                                      axis.title=element_blank(),
                                                      panel.grid.major=element_blank(),
                                                      panel.grid.minor=element_blank(),
                                                      panel.border=element_blank())+
                                      scale_fill_continuous(type = "viridis", option="D")           


