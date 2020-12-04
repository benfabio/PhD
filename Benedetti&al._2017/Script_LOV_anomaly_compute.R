##### 27/02/2015 - LOV - Fabio Benedetti
##### Script for : 
#					- Loading 1965-1994 decadal climatologies from WOA13 (Levin's grid, avgSST, stdSST, avgSSS, avgMLD)
#					- Loading 1965-1994 decadal climatologies from NEMOMED8 (Levin's grid, avgSST, stdSST, avgSSS, avgMLD)
#					- Loading scenario decadal climatologies from NEMOMED8
#					- Compute averages for baseline periods for observed and modelled predictors' values
#					- Compute anomalies (= difference) between scenario runs [(20-50) - (65-94)]
#					- Add this anomaly to observed climatologies 


### Last update = 09/03/2015

# ------------------------------------------------------------------------------------------------------------------------------------

library("rgeos")
library("rgdal")
library("raster")
library("sp")
library("spdep")
library("ncdf4")
library("stringr")
library("spdep")
library("parallel")
library("FactoMineR")
library("automap")
library("ggplot2")
library("reshape2")
library("plyr")
library("dplyr")
library("fields")
library("scales")
library("RColorBrewer")
library("matrixStats")

# ------------------------------------------------------------------------------------------------------------------------------------

cl <- read.csv("gshhg_medit.csv")
names(cl) <- c("lon", "lat")
coast <- list(
  # the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  geom_polygon(aes(x=lon, y=lat), data= cl, fill = "white"),
  geom_path(aes(x=lon, y=lat), data= cl, colour = "black", linetype = 1),
  # appropriate projection
  coord_quickmap(),
  # remove extra space around the coast
  scale_x_continuous(name = "Longitude", breaks = c(0,10,20,30), labels= c("0°","10°E","20°E","30°E"), 
                     expand=c(0,0)), 
  scale_y_continuous(name = "Latitude", breaks = c(30,35,40,45), labels= c("30°N","35°N","40°N","45°N")
    , expand=c(0,0)),
  # dark gray background for the panel and legend
  theme(
    panel.background=element_rect(fill="white"),  # background
    legend.key=element_rect(fill="black"),
    panel.grid.major=element_line(colour="white")
  )
)

# MATLAB style maps:
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



##### 1°) Loading NEMOMED8 climatologies for the 3 decades of the chosen baseline period (65-74, 75-84, 85-94), and compute pluridecadal climatology
dir()


### Load each decadal climatology:
ddf_6574 <- read.table("climato_decadal_EH1.7_Levin_65-74.txt", h=T, sep="\t")
dim(ddf_6574)
ddf_7584 <- read.table("climato_decadal_EH1.7_Levin_75-84.txt", h=T, sep="\t")
dim(ddf_7584)
ddf_8594 <- read.table("climato_decadal_EH1.7_Levin_85-94.txt", h=T, sep="\t")
dim(ddf_8594)

# Check colnames:
colnames(ddf_6574)

mat_6594 <- data.frame(x = ddf_6574$x, y = ddf_6574$y,
			SST_6978 = ddf_6574$avgSST, SST_7988 = ddf_7584$avgSST, SST_8998 = ddf_8594$avgSST,
			stdSST_6978 = ddf_6574$stdSST, stdSST_7988 = ddf_7584$stdSST, stdSST_8998 = ddf_8594$stdSST,
			SSS_6978 = ddf_6574$avgSSS, SSS_7988 = ddf_7584$avgSSS, SSS_8998 = ddf_8594$avgSSS,
			MLD_6978 = ddf_6574$avgMLD, MLD_7988 = ddf_7584$avgMLD, MLD_8998 = ddf_8594$avgMLD )
   
gc()			
			
# ggplot(mat_6594) + geom_point(aes(x=x,y=y,colour= MLD_8998)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
			
mat_6594 <- as.matrix(mat_6594)		
dim(mat_6594)
colnames(mat_6594)	
rm(ddf_6574, ddf_7584, ddf_8594)

### Compute 1965-1994 averages: 
ddf_6594_nemomed8 <- data.frame(x = mat_6594[,"x"], y = mat_6594[,"y"], 
								avgSST = rowMeans(mat_6594[,3:5]), 
								stdSST = rowMeans(mat_6594[,6:8]), 
								avgSSS = rowMeans(mat_6594[,9:11]), 
								avgMLD = rowMeans(mat_6594[,12:14]) )
	
gc()	
			
# Check:			
summary(ddf_6594_nemomed8)
dim(ddf_6594_nemomed8)
rm(mat_6594)

# Maps:
sst <- ggplot(ddf_6594_nemomed8) + geom_point(aes(x=x,y=y,colour= avgSST)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggsave(plot= sst, filename = paste("map_SST_HIS_6594_Levin.png", sep=""), width=13, height=7, dpi= 300)

stdsst <- ggplot(ddf_6594_nemomed8) + geom_point(aes(x=x,y=y,colour= stdSST)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggsave(plot= stdsst, filename = paste("map_stdSST_HIS_6594_Levin.png", sep=""), width=13, height=7, dpi= 300)

sss <- ggplot(ddf_6594_nemomed8) + geom_point(aes(x=x,y=y,colour= avgSSS)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggsave(plot= sss, filename = paste("map_SSS_HIS_6594_Levin.png", sep=""), width=13, height=7, dpi= 300)

maxmld <- ggplot(ddf_6594_nemomed8) + geom_point(aes(x=x,y=y,colour= avgMLD)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggsave(plot= maxmld, filename = paste("map_maxMLD_HIS_6594_Levin.png", sep=""), width=13, height=7, dpi= 300)


### Save climatologies for niche modelling projections:
write.table(ddf_6594_nemomed8, "climatos_HISF_6594_Levin.txt", row.names = FALSE, sep="\t", dec=".")

##### 2°) Loading NEMOMED8 climatologies for the 3 decades of the chosen scenario and compute pluridecadal climatology

##### A) 20-50 period
ddf_2029 <- read.table("climato_decadal_EH4.1_Levin_20-29.txt", h=T, sep="\t")
dim(ddf_2029)
ddf_3039 <- read.table("climato_decadal_EH4.1_Levin_30-39.txt", h=T, sep="\t")
dim(ddf_3039)
ddf_4049 <- read.table("climato_decadal_EH4.1_Levin_40-49.txt", h=T, sep="\t")
dim(ddf_4049)

# Check colnames:
colnames(ddf_4049)

mat_2050 <- data.frame(x = ddf_2029$x, y = ddf_2029$y,
					   SST_2029 = ddf_2029$avgSST, SST_3039 = ddf_3039$avgSST, SST_4049 = ddf_4049$avgSST,
					   stdSST_2029 = ddf_2029$stdSST, stdSST_3039 = ddf_3039$stdSST, stdSST_4049 = ddf_4049$stdSST,
					   SSS_2029 = ddf_2029$avgSSS, SSS_3039 = ddf_3039$avgSSS, SSS_4049 = ddf_4049$avgSSS,
					   MLD_2029 = ddf_2029$avgMLD, MLD_3039 = ddf_3039$avgMLD, MLD_4049 = ddf_4049$avgMLD )
   
gc()			
			
ggplot(mat_2050) + geom_point(aes(x=x,y=y,colour= MLD_4049)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
			
mat_2050 <- as.matrix(mat_2050)		
dim(mat_2050)
colnames(mat_2050)	
rm(ddf_2029, ddf_3039, ddf_4049)

### Compute averages: 
ddf_2050_nemomed8 <- data.frame(x = mat_2050[,"x"], y = mat_2050[,"y"], 
								avgSST = rowMeans(mat_2050[,3:5]), 
								stdSST = rowMeans(mat_2050[,6:8]), 
								avgSSS = rowMeans(mat_2050[,9:11]), 
								avgMLD = rowMeans(mat_2050[,12:14]) )
	
gc()	
			
# Check:			
summary(ddf_2050_nemomed8)
dim(ddf_2050_nemomed8)
rm(mat_2050)

# Maps:
ggplot(ddf_2050_nemomed8) + geom_point(aes(x=x,y=y,colour= avgSST)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggplot(ddf_2050_nemomed8) + geom_point(aes(x=x,y=y,colour= stdSST)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggplot(ddf_2050_nemomed8) + geom_point(aes(x=x,y=y,colour= avgSSS)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggplot(ddf_2050_nemomed8) + geom_point(aes(x=x,y=y,colour= avgMLD)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast



##### B) 69-98 period
ddf_6978 <- read.table("climato_decadal_EH4.1_Levin_69-78.txt", h=T, sep="\t")
dim(ddf_6978)
ddf_7988 <- read.table("climato_decadal_EH4.1_Levin_79-88.txt", h=T, sep="\t")
dim(ddf_7988)
ddf_8998 <- read.table("climato_decadal_EH4.1_Levin_89-98.txt", h=T, sep="\t")
dim(ddf_8998)

# Check colnames:
colnames(ddf_8998)

mat_6998 <- data.frame(x = ddf_6978$x, y = ddf_6978$y,
					   SST_6978 = ddf_6978$avgSST, SST_7988 = ddf_7988$avgSST, SST_8998 = ddf_8998$avgSST,
					   stdSST_6978 = ddf_6978$stdSST, stdSST_7988 = ddf_7988$stdSST, stdSST_8998 = ddf_8998$stdSST,
					   SSS_6978 = ddf_6978$avgSSS, SSS_7988 = ddf_7988$avgSSS, SSS_8998 = ddf_8998$avgSSS,
					   MLD_6978 = ddf_6978$avgMLD, MLD_7988 = ddf_7988$avgMLD, MLD_8998= ddf_8998$avgMLD )
   
gc()			
			
#ggplot(mat_6998) + geom_point(aes(x=x,y=y,colour= MLD_8998)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
			
mat_6998 <- as.matrix(mat_6998)		
dim(mat_6998)
colnames(mat_6998)	
rm(ddf_6978, ddf_7988, ddf_8998)

### Compute averages: 
ddf_6998_nemomed8 <- data.frame(x = mat_6998[,"x"], y = mat_6998[,"y"], 
								avgSST = rowMeans(mat_6998[,3:5]), 
								stdSST = rowMeans(mat_6998[,6:8]), 
								avgSSS = rowMeans(mat_6998[,9:11]), 
								avgMLD = rowMeans(mat_6998[,12:14]) )
	
gc()	
			
# Check:			
summary(ddf_6998_nemomed8)
dim(ddf_6998_nemomed8)
rm(mat_6998)

# Maps:
ggplot(ddf_6998_nemomed8) + geom_point(aes(x=x,y=y,colour= avgSST)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggplot(ddf_6998_nemomed8) + geom_point(aes(x=x,y=y,colour= stdSST)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggplot(ddf_6998_nemomed8) + geom_point(aes(x=x,y=y,colour= avgSSS)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
ggplot(ddf_6998_nemomed8) + geom_point(aes(x=x,y=y,colour= avgMLD)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast





##### 3°) Compute anomaly between scanerio and baseline.
##### Watchout for the right baseline : HIS only for A2 runs, HIS-F for all the others.
ls()

# Gather values from both periods: hist and scenario
mat_anom <- data.frame(x = ddf_2050_nemomed8$x, y = ddf_2050_nemomed8$y,
					   SST_6594 = ddf_6594_nemomed8$avgSST,    SST_2050 = ddf_2050_nemomed8$avgSST,
					   stdSST_6594 = ddf_6594_nemomed8$stdSST, stdSST_2050 = ddf_2050_nemomed8$stdSST,
					   SSS_6594 = ddf_6594_nemomed8$avgSSS,    SSS_2050 = ddf_2050_nemomed8$avgSSS,
					   MLD_6594 = ddf_6594_nemomed8$avgMLD,    MLD_2050 = ddf_2050_nemomed8$avgMLD )
					   
# Compute anomaly
ddf_anom <- data.frame(x = mat_anom[,"x"], y = mat_anom[,"y"], 
								anom_SST = mat_anom[,4] - mat_anom[,3] , 
								anom_stdSST = mat_anom[,6] - mat_anom[,5] , 
								anom_SSS = mat_anom[,8] - mat_anom[,7] , 
								anom_MLD = mat_anom[,10] - mat_anom[,9] )
	
gc()	
			
# Check:			
summary(ddf_anom)
dim(ddf_anom)
rm(mat_anom) 

# Maps:
ggplot(ddf_anom) + geom_point(aes(x=x,y=y,colour= anom_SST)) + scale_colour_distiller(palette = "YlOrRd", guide = "colourbar") + coord_quickmap() + theme_bw() + coast
ggplot(ddf_anom) + geom_point(aes(x=x,y=y,colour= anom_stdSST)) + scale_colour_distiller(palette = "Spectral", guide = "colourbar") + coord_quickmap() + theme_bw() + coast
ggplot(ddf_anom) + geom_point(aes(x=x,y=y,colour= anom_SSS)) + scale_colour_distiller(palette = "YlOrRd", guide = "colourbar") + coord_quickmap() + theme_bw() + coast
ggplot(ddf_anom) + geom_point(aes(x=x,y=y,colour= anom_MLD)) + scale_colour_distiller(palette = "Spectral", guide = "colourbar") + coord_quickmap() + theme_bw() + coast

# Save it: 
write.table(ddf_anom, paste("climato_anomalies_B1ARF_2050_Levin.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)










