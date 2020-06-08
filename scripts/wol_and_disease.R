########################################################################################################################

### WOL lethality analysis to inform on Mendelian phenotypes
### Script: wol.disease.R
### Purpose: explore WoL in the context of Mendelian disease, MoI, type of disease
### Resources: HPO, panelApp, gnomad 
### Latest update : 13.05.2020

## load packages #######################################################################################################

library(plyr);library(dplyr);library(tidyr);library(stringr);library(readr);
library(ggplot2);library(cowplot);library(tidyverse);
library(waffle);library(magrittr);library(hrbrthemes);
library(ggridges)


#########################################################################################################################

## summary viability data DR111 (barplot)  ##############################################################################

dr11.via <- data.frame(Category =c("Lethal","Subviable","Viable"),
                       Number.genes = c(1405,525,3905)) %>%
  mutate(Category = factor(Category,levels =c("Lethal","Subviable","Viable")))



ggplot(dr11.via , aes(x = Category, y = Number.genes,fill=Category)) + 
  geom_bar(stat='identity') + 
  scale_fill_manual(values=c("#D55E00", "#E69F00","#009E73")) +
  theme_minimal() +
  xlab("")+
  ylab("number of genes") +
  theme(axis.text.x = element_text(size=12,face="bold")) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(legend.position="bottom") 



## import windows of lethality (waffle plot)  ###########################################################################

e <- read.delim("./data/embryo.windows.to.share.30042019.txt",stringsAsFactors = F) %>%
  dplyr::select(MGI_ID,HGNC_ID,Window_Category) %>%
  filter(Window_Category!="-") %>%
  distinct()


e %>% group_by(Window_Category) %>% tally()


waffle(c('Early-gestation' = 197, 'Mid-gestation' = 50,
         'Late-gestation' = 153), rows = 17, size=1,pad=1,colors=c("#DCAE1D",  "#7A9D96","#00303F"),
       title="",legend_pos="bottom",flip=TRUE)



#########################################################################################################################

### HPO annotations #####################################################################################################


e.hpo.annotations <- read.delim("./data/aux_data/wol.hpo.annotations.txt",sep="\t",
                                stringsAsFactors = F,header=T) %>%
  group_by(Window_Category,hpo_genes) %>%
  tally()

e.hpo.moi.genes <- read.delim("./data/aux_data/wol.hpo.annotations.txt",sep="\t",
                              stringsAsFactors = F,header=T)  %>%
  group_by(moi,Window_Category) %>%
  tally() %>%
  filter(!is.na(moi)) %>%
  mutate(total = c(77,76,23)) %>%
  mutate(percentage = n/total*100) %>%
  mutate(Window_category = factor(Window_Category,
                                  levels=c("Early-gestation","Mid-gestation","Late-gestation"))) 


ggplot(e.hpo.moi.genes, aes(x = Window_category, y = percentage,fill=Window_category)) + 
  geom_bar(stat='identity') + 
  scale_fill_manual(values=c("#DCAE1D",  "#7A9D96","#00303F")) +
  facet_grid(~moi) +
  theme_minimal() +
  xlab("")+
  ylab("% of genes") +
  theme(strip.text.x = element_text(size = 12,face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size=11,face="bold")) +
  theme(axis.title.y = element_text(size=12)) +
  theme(legend.position="bottom") +
  labs(fill ="Window of Lethality") 


#########################################################################################################################

## gnomad scores ########################################################################################################


e.gnomad <- read.delim("./data/aux_data/wol.gnomad.scores.txt",sep="\t",
                       stringsAsFactors = F,header=T) %>%
  mutate(Window_Category = factor(Window_Category,
                                  levels=c("Early-gestation",
                                           "Mid-gestation","Late-gestation")))


ggplot(e.gnomad, aes(x = pLI, y = Window_Category, fill = Window_Category)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  xlim(0, 1) +
  scale_fill_manual(values =c("#DCAE1D", "#7A9D96","#00303F")) +
  ylab("Window of Lethality") +
  theme(axis.text.y = element_text(size=12,face="bold")) 



ggplot(e.gnomad, aes(x = pRec, y = Window_Category, fill = Window_Category)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  xlim(0, 1) +
  scale_fill_manual(values =c("#DCAE1D", "#7A9D96","#00303F")) +
  ylab("Window of Lethality") +
  theme(axis.text.y = element_text(size=12,face="bold")) 


#########################################################################################################################

