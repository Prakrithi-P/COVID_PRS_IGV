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

######
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
####op warning: Warning messages:
#1: The argument auto.palette.mapping is deprecated. Please use midpoint for numeric data and stretch.palette for categorical data to control the palette mapping. 
#2: Currect projection of shape shp unknown. Long-lat (WGS84) is assumed. 

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

# Finally, we'll clip the tessellated  surface to the Texas boundaries
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
  tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)+tm_legend(legend.outside=TRUE)+tm_text("POP", just="top", xmod=0.7, size = 0.6)
