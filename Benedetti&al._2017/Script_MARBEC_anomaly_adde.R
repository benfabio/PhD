##### 05/03/2015 - MARBEC - Fabio Benedetti
##### Script for : 
#					- Loading 1965-1994 decadal climatologies from WOA13 (Levin's grid, avgSST, stdSST, avgSSS, avgMLD)
#					- Loading anomalies for all scenario runs
#					- Add this anomaly to observed climatologies 
#					- Correct MLD if negative ! 
#					- Map fields


### Last update = 05/03/2015

# ------------------------------------------------------------------------------------------------------------------------------------


library("rgdal")
library("raster")
library("sp")
library("spdep")
library("stringr")
library("parallel")
library("ggplot2")
library("reshape2")
library("plyr")
library("fields")
library("scales")
library("RColorBrewer")
library("matrixStats")


# ------------------------------------------------------------------------------------------------------------------------------------

### Load Mediterranean coastline
cl <- read.csv("gshhg_medit.csv")
names(cl) <- c("lon", "lat")
coast <- list(
 	 	# the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  	  	geom_polygon(aes(x=lon, y=lat), data= cl, fill = "white"),
  	  	geom_path(aes(x=lon, y=lat), data= cl, colour = "black", linetype = 1),
  	  	# appropriate projection
  	  	coord_map(),
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

##### 1°) Loading WOA13 climatologies for the 3 decades of the chosen baseline period (65-74, 75-84, 85-94), and compute pluridecadal climatology
dir()

### Build a ddf containing predictors' values averaged over the each decades
### A) 1965-1974	
# Read avgSST
avgSST <- read.table("avgSST_6574_woa13_Levin.txt", h=T, sep="\t")
# Read stdSST
stdSST <- read.table("stdSST_6574_woa13_Levin.txt", h=T, sep="\t")
# Read avgSSS
avgSSS <- read.table("avgSSS_6574_woa13_Levin.txt", h=T, sep="\t")
# read avgMLD
maxMLD <- read.table("maxMLD_6574_woa13_Levin.txt", h=T, sep="\t")
		
ddf_6574 <- data.frame(x = avgSST$x, y = avgSST$y, avgSST = avgSST$avgSST, stdSST = stdSST$stdSST, avgSSS = avgSSS$avgSSS, maxMLD = maxMLD$maxMLD)
rm('avgSST','avgSSS','maxMLD','stdSST')



### B) 1975-1984	
# Read avgSST
avgSST <- read.table("avgSST_7584_woa13_Levin.txt", h=T, sep="\t")
# Read stdSST
stdSST <- read.table("stdSST_7584_woa13_Levin.txt", h=T, sep="\t")
# Read avgSSS
avgSSS <- read.table("avgSSS_7584_woa13_Levin.txt", h=T, sep="\t")
# read avgMLD
maxMLD <- read.table("maxMLD_7584_woa13_Levin.txt", h=T, sep="\t")
		
ddf_7584 <- data.frame(x = avgSST$x, y = avgSST$y, avgSST = avgSST$avgSST, stdSST = stdSST$stdSST, avgSSS = avgSSS$avgSSS, maxMLD = maxMLD$maxMLD)
rm('avgSST','avgSSS','maxMLD','stdSST')



### C) 1985-1994	
# Read avgSST
avgSST <- read.table("avgSST_8594_woa13_Levin.txt", h=T, sep="\t")
# Read stdSST
stdSST <- read.table("stdSST_8594_woa13_Levin.txt", h=T, sep="\t")
# Read avgSSS
avgSSS <- read.table("avgSSS_8594_woa13_Levin.txt", h=T, sep="\t")
# read avgMLD
maxMLD <- read.table("maxMLD_8594_woa13_Levin.txt", h=T, sep="\t")
		
ddf_8594 <- data.frame(x = avgSST$x, y = avgSST$y, avgSST = avgSST$avgSST, stdSST = stdSST$stdSST, avgSSS = avgSSS$avgSSS, maxMLD = maxMLD$maxMLD)
rm('avgSST','avgSSS','maxMLD','stdSST')


mat_6594 <- data.frame(x = ddf_6574$x, y = ddf_6574$y,
			SST_6574 = ddf_6574$avgSST, SST_7584 = ddf_7584$avgSST, SST_8594 = ddf_8594$avgSST,
			stdSST_6574 = ddf_6574$stdSST, stdSST_7584 = ddf_7584$stdSST, stdSST_8594 = ddf_8594$stdSST,
			SSS_6574 = ddf_6574$avgSSS, SSS_7584 = ddf_7584$avgSSS, SSS_8594 = ddf_8594$avgSSS,
			MLD_6574 = ddf_6574$maxMLD, MLD_7584 = ddf_7584$maxMLD, MLD_8594 = ddf_8594$maxMLD )
   
gc()			
			
mat_6594 <- as.matrix(mat_6594)		
dim(mat_6594)
colnames(mat_6594)	
rm(ddf_6574, ddf_7584, ddf_8594)


### Compute 1965-1994 averages: 
ddf_6594_woa13 <- data.frame(x = mat_6594[,"x"], y = mat_6594[,"y"], 
			avgSST = rowMeans(mat_6594[,3:5]), 
			stdSST = rowMeans(mat_6594[,6:8]), 
			avgSSS = rowMeans(mat_6594[,9:11]), 
			maxMLD = rowMeans(mat_6594[,12:14]) )
	
gc()	
			
# Check:			
summary(ddf_6594_woa13)

mat_woa13 <- as.matrix(ddf_6594_woa13)

##### 2°) Loading NM8 anomalies and add them to WOA13

scenario <- c("A2","A2F","A2RF","A2ARF","A1BARF","B1ARF")
period <- c("2050","6998")

### For testing:
# s <- "A2"
# p <- "6998"


### In a for loop : extract anomalies and add them

for(s in scenario) {
	
		for(p in period) {
		
			### reading anomalies from NM8
			message( paste("Reading",s,p,"anomalies", sep="_") )
			ddf_anom <- read.table(paste("climato_anomalies",s,p,"Levin.txt", sep= "_"), h=T, sep="\t", dec=".")
			# summary(ddf_anom)
			# summary(ddf_6594_woa13)
			mat_anom <- as.matrix(ddf_anom)
			# colnames(mat_anom)
		
			### create new ddf with new predictor variables 
			ddf <- data.frame(x = mat_woa13[,"x"], y = mat_woa13[,"y"], 
							  avgSST = ( mat_anom[,"anom_SST"] + mat_woa13[,"avgSST"] ), 
							  stdSST = ( mat_anom[,"anom_stdSST"] + mat_woa13[,"stdSST"] ), 
							  avgSSS = ( mat_anom[,"anom_SSS"] + mat_woa13[,"avgSSS"] ), 
							  maxMLD = ( mat_anom[,"anom_MLD"] + mat_woa13[,"maxMLD"] ) )
	
			gc()
			
			# Check:			
			# summary(ddf)
			
			# save :
			write.table(ddf, paste("clims",s,p,"final","Levin.txt",sep="_"), sep="\t", row.names=FALSE, quote=FALSE)
			
			### BONUS
			### correct MLD values : negative values in the GoL and in Otranto are due to 
			# - overconvective behaviour of NM8 under the old RCM3 forcings
			# - underestimation of WOA13 climatologies (missing convective events)
			
			# Get absolute values from the negative ones: 
			# vals <- abs(ddf[which(ddf$maxMLD < 0),"maxMLD"])
			# Get row indices of the cells where these negative values occur:
			# ids <- which(ddf$maxMLD < 0)
			# Add absolute value to the negative ones:
			# ddf[ids,"maxMLD"] <- (ddf[ids,"maxMLD"]) + vals
			
			message( paste("Saving_maps",s,p, sep="_") )
			# Maps & save 
			sst <- ggplot(ddf) + geom_point(aes(x=x,y=y, colour = avgSST)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
			ggsave(plot= sst, filename = paste("map_avgSST",s,p,"final_Levin.png", sep="_"), width=13, height=7, dpi= 300)

			stdsst <- ggplot(ddf) + geom_point(aes(x=x,y=y, colour = stdSST)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
			ggsave(plot= stdsst, filename = paste("map_stdSST",s,p,"final_Levin.png", sep="_"), width=13, height=7, dpi= 300)

			sss <- ggplot(ddf) + geom_point(aes(x=x,y=y, colour = avgSSS)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
			ggsave(plot= sss, filename = paste("map_avgSSS",s,p,"final_Levin.png", sep="_"), width=13, height=7, dpi= 300)

			maxmld <- ggplot(ddf) + geom_point(aes(x=x,y=y, colour = maxMLD)) + scale_colour_gradientn(colours = jet.colors(7)) + coord_quickmap() + theme_bw() + coast
			ggsave(plot= maxmld, filename = paste("map_maxMLD",s,p,"final_Levin.png", sep="_"), width=13, height=7, dpi= 300)
			
			
			# Make room:
			rm(ddf,mat_anom,ddf_anom,sst,stdsst,sss,maxmld)
					
		} #eo period
	
} #eo scenario 


### Check for negative MLD values: 
length(ddf[which(ddf$maxMLD <= 0),"maxMLD"])
ggplot(ddf) + geom_point(aes(x=x,y=y,colour= maxMLD)) + scale_colour_gradient2(low="blue", high="red",midpoint= 0, guide = "colourbar") + coord_quickmap() + theme_bw() + coast
ggplot(ddf) + geom_point(aes(x=x,y=y,colour= maxMLD)) + scale_colour_gradientn(colours=jet.colors(7), guide = "colourbar") + coord_quickmap() + theme_bw() + coast



