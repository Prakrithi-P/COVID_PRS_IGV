library(ggplot2)
library(forcats)
library(AICcmodavg)
library(ggpubr)
library(dplyr)

##PRS- PLINK
system(paste(“plink --bfile ../../../IGV_7l --score info.txt 1 2 3 header --extract snps --out PRS”, sep=” “)   #Fields in info.txt - SNP Risk_Allele OR
## PRS.profile is updated with Population labels (POP)
d<-read.csv("PRS.profile", sep="\t", header=T) ##reading PRS scores 
##median scores for individuals in each population
dd<-d %>%
    group_by(POP) %>%
    summarise(median(SCORE))
    
##boxplot - Distribution of PRS scores across the populations
d%>% ggplot(aes(x=fct_reorder(POP,PRS,.desc = T), y=PRS,fill="red")) + geom_boxplot()+theme(axis.text.x = element_text(angle=60,vjust = 0.7),plot.title = element_text(hjust = 0.5))+geom_jitter(width=0.2,alpha=0.4)+ggtitle("Distribution of PRS scores across IGV Populations")+xlab("Population")

##test for significance- ANOVA & Tukey’s HSD Test
AN1<-aov(PRS~POP, d)
summary(AN1)
thsd<-TukeyHSD(AN1)
tkd<-as.data.frame(thsd[1:1])
write.table(tkd,"ANOVA_1way_tukey_hsd", sep="\t",quote = F, row.names = T)
par(las=1)
par(mar=c(5,15,5,1))
plot(thsd)

## District level COVID19 information for corresponding IGV pops collected from https://www.covid19india.org/ and https://covid19.assam.gov.in/district/
#PRS-Mortality - Spearman correlation
c<-read.csv("prs_deaths", sep="\t",header=T)  ##Fields in prs_deaths - POP,PRS,District,No.of Deaths
ggscatter(c,x="PRS",y="Deaths",add="reg.line",conf.int=T, cor.coef=T, cor.method = "spearman",ylab = "No.of deaths due to COVID",title = "Spearman correlation")  ##correlation plot
