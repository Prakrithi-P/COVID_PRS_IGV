library(ggplot2)
library(forcats)
library(AICcmodavg)
library(ggpubr)
library(dplyr)

##PRS- PLINK
system(paste(“plink --bfile ../../../IGV_7l --score top100_inIGV_summary 1 2 3 header --extract snps --out PRS”, sep=” “)   #Fields in info.txt - SNP Risk_Allele OR
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
#PRS-Mortality - Pearson's correlation
c<-read.csv("prs_deaths", sep="\t",header=T)  ##Fields in prs_deaths - POP,PRS,District,No.of Deaths
ggscatter(c,x="PRS",y="Deaths",add="reg.line",conf.int=T, cor.coef=T, cor.method = "pearson",ylab = "No.of deaths due to COVID")  ##correlation plot

       
########################################################
## VarLD - Comparison of LD around SNPs
##Genotype file 12 format
system(paste(“/opt/apps/plink-1.07-x86_64/plink --bfile CEU_1000G_all  --extract range top100_snps_10MB_list  --recode12 --transpose --tab --noweb --out CEU_100snps_10mb”, sep=” “)
system(paste(“/opt/apps/plink-1.07-x86_64/plink --bfile CEU_1000G_all  --extract range top100_snps_10MB_list  --recode12 --transpose --tab --noweb --out CEU_100snps_10mb”, sep=” “)

##Input SNP pos genotype1 genotype2 …
##Convert Genotype file to 1234 format -- for each chromosome separately
system(paste(“sed  -r 's/(\s+)?\S+//3' CEU_rs10735079_10mb.tped | sed  -r 's/(\s+)?\S+//1' | sed 's/2 2/1/g'  | sed 's/1 1/3/g' | sed 's/2 1/2/g' | sed 's/1 2/2/g' | sed 's/2 0/4/g' | sed 's/0 2/4/g' | sed 's/0 1/4/g' | sed 's/0 0/4/g' | sed 's/\trs/rs/g' > CEU_rs10735079_infile”, sep=” “)
system(paste(“sed  -r 's/(\s+)?\S+//3' ITU_rs10735079_10mb.tped | sed  -r 's/(\s+)?\S+//1' | sed 's/2 2/1/g'  | sed 's/1 1/3/g' | sed 's/2 1/2/g' | sed 's/1 2/2/g' | sed 's/2 0/4/g' | sed 's/0 2/4/g' | sed 's/0 1/4/g' | sed 's/0 0/4/g' | sed 's/\trs/rs/g' > ITU_rs10735079_infile”, sep=” “)
system(paste(“sed  -r 's/(\s+)?\S+//3' PJL_rs10735079_10mb.tped | sed  -r 's/(\s+)?\S+//1' | sed 's/2 2/1/g'  | sed 's/1 1/3/g' | sed 's/2 1/2/g' | sed 's/1 2/2/g' | sed 's/2 0/4/g' | sed 's/0 2/4/g' | sed 's/0 1/4/g' | sed 's/0 0/4/g' | sed 's/\trs/rs/g' > PJL_rs10735079_infile”, sep=” “)

##VarLD: 
system(paste(“for i in {1..22} ; do java -jar rgenetics-1.0.jar -p VarLD CEU_chr$i\_infile ITU_chr$i\_infile -o CEU_ITU_chr$i ; done",sep=" ")
#system(paste(“for i in {1..22} ; do java -jar rgenetics-1.0.jar -p VarLD CHB_chr$i\_infile ITU_chr$i\_infile -o CEU_ITU_chr$i ; done",sep=" ")

##standardize the varLD output with the R script provided by the VarLD package as below     
Standardization
varLD.out <- {}
chr.store <- {}
for (chr in 1:22){
    varLD.temp <- read.table(paste("CEU_ITU_chr", chr, sep=""), sep="\t", header = T)
varLD.out <- rbind(varLD.out, varLD.temp)
chr.store <- c(chr.store, rep(chr, dim(varLD.temp)[1]))
print(paste("completed reading in unstandardized varLD output file for chromosome ", chr, sep=""))
}
varLD.mean <- mean(varLD.out[,"raw_score"])
varLD.sd <- sd(varLD.out[,"raw_score"])
standardized_score <- (varLD.out[,"raw_score"] - varLD.mean)/varLD.sd
varLD.out <- cbind(varLD.out, standardized_score)
varLD.threshold <- quantile(standardized_score, probs = percentile.out)
for (chr in 1:22){
    chr.flag <- which(chr.store == chr)
    write.table(varLD.out[chr.flag,], paste("CEU_ITU_chr", chr, "_standardized.out", sep=""), sep="\t", quote=F, row.names=F)
    print(paste("completed writing out standardized varLD output file for chromosome ", chr, sep=""))
}
n.length.percentile <- length(percentile.out)
for (i in 1:n.length.percentile){
    print(paste("varLD threshold for ", percentile.out[i], " = ", varLD.threshold[i], sep=""))
} 

##varLD thresholds
[1] "varLD threshold for 0.95 = 1.78136769083323"
[1] "varLD threshold for 0.99 = 3.37859820232148"
[1] "varLD threshold for 0.999 = 6.87171244835101"
[1] "varLD threshold for 0.9999 = 12.4041619384797"

## To show variants and their LD patterns across multiple chromosomes
library(qqman)
##add chromosome labels and concatenate the chr-wise standardized output files and use as follows:       
d<-read.csv("varld_standardized.out", header=TRUE, sep="\t")
h<-read.csv("highlight_100", sep="\t", header=T). ## list of positions to highlight
hs<-as.character(h$pos)
manhattan(d, chr = "Chr", bp = "pos", p = "Standardized_score", snp = "pos",
col = c("gray10", "gray60"), chrlabs = NULL,suggestiveline = 3.37, genomewideline = 6.87,highlight = hs, logp = FALSE,ylim=c(-2,20), ylab=”Standardized VarLD score”)
##suggestiveline, genomewideline values are based on the varld standardization script
       
#########################################################
##Spatial Plot-- Modified verion of the already avalable code based on IDW algorithm
library(rgdal)
library(tmap)
library(maptools)
library(tmap)
library(spatstat)
library(gstat) # Use gstat's idw routine
library(sp)    # Used for the spsample function
library(raster)
library(rgeos)
setwd("D:/India Shape/")
s<-read.csv("PRS_lat_long",sep="\t", header=TRUE)  ##PRS scores with corresponding coordinates
UTM32n <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# World Geographic System 1984 (lat/long) - mapping
WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
p <-  SpatialPointsDataFrame(coords = s[,c("Longitude", "Latitude")], 
                                  data = s, 
                                  proj4string = UTM32n)
# Load India boudary map
shp<-readOGR("D:/ancestry/India Shape/IGISMAP/Indian_States.shp")

# Replace point boundary extent with that of India
p@bbox <- shp@bbox

tm_shape(shp) + tm_polygons() +
  tm_shape(p) +
  tm_dots(col="PRS", palette = "YlOrRd",
          title="PRS", size=0.7) +
  tm_text("PRS", just="left", xmod=.5, size = 0.7) +
  tm_legend(legend.outside=TRUE)

th  <-  as(dirichlet(as.ppp(p)), "SpatialPolygons")

# The dirichlet function does not carry over projection information
# requiring that this information be added manually
proj4string(th) <- proj4string(p)

# The tessellated surface does not store attribute information
# from the point data layer. We'll use the over() function (from the sp
# package) to join the point attributes to the tesselated surface via
# a spatial join. The over() function creates a dataframe that will need to
# be added to the `th` object thus creating a SpatialPolygonsDataFrame object
th.z     <- over(th, p, fn=mean)
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)

# Finally, we'll clip the tessellated  surface to the India boundaries
set_RGEOS_CheckValidity(2L)
th.clp   <- intersect(shp,th.spdf)

# Map the data
tm_shape(th.clp) + 
  tm_polygons(col="PRS", palette="Reds", auto.palette.mapping=FALSE,
              title="PRS") +
  tm_legend(legend.outside=TRUE)

grd              <- as.data.frame(spsample(p, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(p) <- proj4string(p) # Temp fix until new proj env is adopted
proj4string(grd) <- proj4string(p)

# Interpolate the grid cells using a power value of 2 (idp=2.0)
P.idw <- gstat::idw(PRS ~ 1, p, newdata=grd, idp=2.0)

# Convert to raster object then clip to Texas
r       <- raster(P.idw)
r.m     <- mask(r, shp)

# Plot
tm_shape(r.m) + 
  tm_raster(n=10,palette = "Reds", auto.palette.mapping = FALSE,
            title="PRS") + 
  tm_shape(p) + tm_dots(size=0.02) +
  tm_legend(legend.outside=TRUE)+tm_legend(legend.outside=TRUE)+tm_text("POP", just="top", xmod=0.7, size = 0.6)
