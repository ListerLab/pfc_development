library(readxl)
library(ggplot2)

life_expectancy <- read_xlsx("annotation/Life_expectancy.xlsx")
life_expectancy$lives <- life_expectancy$`live men`+ life_expectancy$`live women`
life_expectancy1 <- life_expectancy[life_expectancy$ages>=90,]
life_expectancy1$probability <- life_expectancy1$lives/sum(life_expectancy1$lives)
life_expectancy2 <- life_expectancy[life_expectancy$ages>=18 & life_expectancy$ages<=79,]
life_expectancy2$probability <- life_expectancy2$lives/sum(life_expectancy2$lives)

atac <- read_xlsx("annotation/Single cell studies database-Specifics.xlsx", sheet=2)
atac$Age <- as.numeric(atac$Age)
atac_fetal <- atac[atac$Fetal=="Yes",]
atac_postnatal <- atac[atac$Fetal=="No",]

gg_fetal_atac <- ggplot(atac_fetal, aes(x=`Age`)) + 
    geom_histogram(bins=34, fill="#512568") + 
    xlim(3,38) + ylim(0,16) + theme_classic() + xlab("Age (pwc)") 
ggsave(gg_fetal_atac, file="supp_figures/SuppFig1_fetal_studies_atac.svg", height=1.68, width=2.96)

gg_adult_atac <- ggplot(atac_postnatal, aes(x=`Age`)) + geom_histogram(bins=83, fill="#F9E51B") + 
    xlim(20,102) + theme_classic() + ylim(0,15) +  xlab("Age (years)") 
ggsave(gg_adult_atac, file="supp_figures/SuppFig1_adult_studies_atac.svg", height=1.68, width=2.96)

postnatal <- read_xlsx("annotation/Single cell studies database-Specifics.xlsx", sheet=3)
ind  <- which(postnatal$Years=="90+")
ind1 <- which(postnatal$Years=="18-79")
set.seed(14)
postnatal$Years[ind] <- sample(life_expectancy1$ages, prob=life_expectancy1$probability, 
                               size=length(ind), replace=TRUE)
postnatal$Years[ind1] <- sample(life_expectancy2$ages, prob=life_expectancy2$probability, 
                               size=length(ind1), replace=TRUE)
fetal <- read_xlsx("annotation/Single cell studies database-Specifics.xlsx", sheet=1)
postnatal$Years <- as.numeric(as.character(postnatal$Years))
fetal$`Age (pcw)` <- as.numeric(as.character(fetal$`Age (pcw)`))
                              

gg_fetal <- ggplot(fetal, aes(x=`Age (pcw)`)) + 
    geom_histogram(bins=34, fill="#512568") + 
    xlim(3,38) + ylim(0,16) + theme_classic() 
ggsave(gg_fetal, file="supp_figures/SuppFig1_fetal_studies.svg", height=1.68, width=2.96)

child <- postnatal[postnatal$Years>1 & postnatal$Years<10,]
gg_childhood <- ggplot(child, aes(x=Years)) + 
    geom_histogram(bins=10, fill="#20988C") + 
    xlim(1,10) + ylim(0,15) + theme_classic() 
ggsave(gg_childhood, file="supp_figures/SuppFig1_child_studies.svg", height=1.68, width=2.96)

adol <- postnatal[postnatal$Years>=10 & postnatal$Years<20,]
gg_adol <- ggplot(adol, aes(x=Years)) + geom_histogram(bins=11, fill="#98CA43") + 
    xlim(10,20) + theme_classic() + ylim(0,15)
ggsave(gg_adol, file="supp_figures/SuppFig1_adol_studies.svg", height=1.68, width=2.96)

adult <- postnatal[postnatal$Years>=20,]
gg_adult <- ggplot(adult, aes(x=Years)) + geom_histogram(bins=83, fill="#F9E51B") + 
    xlim(20,102) + theme_classic() + ylim(0,15)
ggsave(gg_adult, file="supp_figures/SuppFig1_adult_studies.svg", height=1.68, width=2.96)

