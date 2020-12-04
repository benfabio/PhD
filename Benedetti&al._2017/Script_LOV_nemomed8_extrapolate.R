##### 22/01/2015 - LOV - Fabio Benedetti
##### Script for : 
#					- Extracting NEMOMED8 21st century projections from netCDF files, while retrieving proper dates
#					- Overing the climatologies on Noam Levin's grid (over() fun) 
#					- Re-interpolate the 1/8° fields on Levin's 0.1°x0.1° grid
#					- Compute decadal seasonal climatologies and mean climatologies (avgSST, avgSSS, avgMLD, stdSST)

### Last update = 27/01/2015

# ------------------------------------------------------------------------------------------------------------------------------------


library(rgeos)
library(rgdal)
library(raster)
library(sp)
library(spdep)
library(ncdf4)
library(stringr)
library(spdep)
library(parallel)
library(FactoMineR)
library(automap)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(scales)
library(fields)
library(RColorBrewer)


# ------------------------------------------------------------------------------------------------------------------------------------


###### 1°) Loading Levin's grid (shapefile)
Layer <- ogrListLayers("grille/grille_med.shp")
grilleMed <- readOGR("grille/grille_med.shp", Layer)
#summary(grilleMed)
# check dimensions : min long, min lat etc.
#P4S <- CRS("+proj=longlat +datum=WGS84 +no_defs")
#grilleMed <- spTransform(grilleMed, P4S) # changing projection system to longlat 

#coordsLevin <- coordinates(grilleMed)
#class(coordsLevin)
#dim(coordsLevin)
#colnames(coordsLevin)[1:2] <- c("x","y")
#head(coordsLevin)
#ggplot(as.data.frame(coordsLevin)) + geom_point(aes(x=x,y=y)) + theme_bw() + coord_quickmap()


###### 2°) Extract NEMOMED8 fields and print raw climatology (containing all month and years) - not re-interpolated on Levin's grid yet

### for looping:

scenario <- c("EH3","EH3.1","EH3.2","EH4.1")
var <- c("SST","SSS","MLD")


# For testin':
#s <- "EH2.1"
#v <- "SST"


for(s in scenario) {
	
	
				for(v in var) {
	
					### Open netCDF proper file and retrieve coordinates:
					message(paste("Reading netCDF:",s, sep="_"))
					library(ncdf)
					
					if(s == "DE9") {
					nc1 <- open.ncdf("NM8-DE9_med_y2000-2099_2D_monmean.nc")
					lat <- get.var.ncdf(nc1,"nav_lat")
					lon <- get.var.ncdf(nc1,"nav_lon")
					dataClimato <- data.frame("x"=c(lon),"y"=c(lat))
					#head(dataClimato)
					rm('lon','lat')
		
					} else if (s == "EH2.1") {
					nc1 <- open.ncdf("NM8-EH2.1_med_y2001-2099_2D_monmean.nc")
					lat <- get.var.ncdf(nc1,"nav_lat")
					lon <- get.var.ncdf(nc1,"nav_lon")
					dataClimato <- data.frame("x"=c(lon),"y"=c(lat))
					#head(dataClimato)
					rm('lon','lat')
					
					} else if (s == "EH3") {
					nc1 <- open.ncdf("NM8-EH3_med_y2001-2099_2D_monmean.nc")
					lat <- get.var.ncdf(nc1,"nav_lat")
					lon <- get.var.ncdf(nc1,"nav_lon")
					dataClimato <- data.frame("x"=c(lon),"y"=c(lat))
					#head(dataClimato)
					rm('lon','lat')
					
					} else if (s == "EH3.1") {
					nc1 <- open.ncdf("NM8-EH3.1_med_y2001-2099_2D_monmean.nc")
					lat <- get.var.ncdf(nc1,"nav_lat")
					lon <- get.var.ncdf(nc1,"nav_lon")
					dataClimato <- data.frame("x"=c(lon),"y"=c(lat))
					#head(dataClimato)
					rm('lon','lat')
					
					} else if (s == "EH3.2") {
					nc1 <- open.ncdf("NM8-EH3.2_med_y2001-2099_2D_monmean.nc")
					lat <- get.var.ncdf(nc1,"nav_lat")
					lon <- get.var.ncdf(nc1,"nav_lon")
					dataClimato <- data.frame("x"=c(lon),"y"=c(lat))
					#head(dataClimato)
					rm('lon','lat')
					
					} else if (s == "EH4.1") {
					nc1 <- open.ncdf("NM8-EH4.1_med_y2001-2099_2D_monmean.nc")
					lat <- get.var.ncdf(nc1,"nav_lat")
					lon <- get.var.ncdf(nc1,"nav_lon")
					dataClimato <- data.frame("x"=c(lon),"y"=c(lat))
					#head(dataClimato)
					rm('lon','lat')
		
		
					} #eo if loop
	

					message(paste("Extracting variable from netCDF:",s,v, sep="_"))
					
					###### Extract proper variable from nc
					if(v == "SST") { vari <- get.var.ncdf(nc1, "sosstsst") } else if(v == "SSS") { vari <- get.var.ncdf(nc1, "sosaline") } else if (v == "MLD") { vari <- get.var.ncdf(nc1, "somixhgt") }

					### Create data frame of monthly climatologies
					l <- mclapply(seq_len(1188), function(i) cbind(c(vari[,,i])))  ######### !!! if s == "DE9" then : seq_len(1200) !!!
					data <- data.frame(l)
					dataClimato <- cbind(dataClimato, data)
					gc()
					rm('vari','l','data')


					### Change column names to month-year
					message(paste("Changing colnames",s,v, sep="_"))
					
					library(ncdf4)
					nc2 <- nc_open("NM8-EH3_med_y2001-2099_2D_monmean.nc")
					sl1 <- nc_slice(nc2, "sosstsst")
					#dim(sl1) # retrieves : lon, lat and time (= month & year)
					# Melt
					mddf <- melt(sl1, .("x","y","time"), na.rm= F)
					# cast to obtain time as columns, as a function of coordinates
					ddf <- dcast(mddf, x+y~time, na.rm=F) 
					ddf <- ddf[,-c(1:2)]
					gc()

					#length(unique(colnames(ddf[,c(1:480)])))
					# https://stat.ethz.ch/R-manual/R-devel/library/base/html/as.POSIXlt.html
					nemoTime <- unique(as.numeric(colnames(ddf[,c(1:1188)])))
					#nc2
					nemoDates <- as.POSIXct(nemoTime, origin="2001-01-01 00:00:0")
					#class(nemoDates)
					nemoDates <- format(nemoDates, format="%b %Y") 
					nemoDates <- str_replace_all(nemoDates, " ","_")  # replace spaces by underscores
					### Change dataClimato's colnames:
					colnames(dataClimato)[3:1190] <- nemoDates # yeah
					#colnames(dataClimato)
					#summary(dataClimato)
					rm('sl1','mddf','nemoTime','nemoDates','ddf','nc2')
					gc()


					###### 4°) Over NEMOMED8 fields on Levin's grid
					#ggplot(dataClimato) + geom_point(aes(x=x,y=y)) + theme_bw() ### NEMOMED8's grid
					#ggplot(dataClimato) + geom_point(aes(x=x,y=y, colour= Dec_2000)) + theme_bw() # good !
					#class(dataClimato)

					ddf <- as.matrix(dataClimato)
					#dim(ddf)

					### Get rid of NAs - but 'ddf' is very large (trouble subsetting) : divide is into 4 periods (120 columns each), get rid pf NAs, and then cbind
					ddf1 <- na.omit(ddf[,c(1:300)])
					ddf2 <- na.omit(ddf[,c(301:601)])
					ddf3 <- na.omit(ddf[,c(602:902)])
					ddf4 <- na.omit(ddf[,c(903:1190)])

					ddf <- cbind(ddf1, ddf2, ddf3, ddf4)
					ddf <- data.frame("ID"=seq(1:nrow(ddf)), ddf)
					rm('dataClimato','ddf1','ddf2','ddf3','ddf4')
					gc()

					### Create the shapefile with spTransform()
					pointsClimato <- SpatialPointsDataFrame(coords=ddf[,2:3], data=ddf[,c(2:1191)], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
					pointsClimato <- spTransform(pointsClimato, CRS(proj4string(grilleMed)))
					gc()
					#writeOGR(pointsClimato, "./2_processed/shp","Points_climato_proj", driver="ESRI Shapefile", overwrite_layer=TRUE)

					### Retrieve climatological values per Levin's cell
					message(paste("Overing",s,v, sep="_"))
					
					ptsInPol <- over(grilleMed, pointsClimato, fn= mean) # On prend la moyenne car certaines cellules contiennent plus d'un point
					#class(ptsInPol)
					#dim(ptsInPol)
					#summary(ptsInPol)
					gc()

					# ggplot(ptsInPol) + geom_point(aes(x=x, y=y, colour= Nov_2030)) + theme_bw() + coord_quickmap() + scale_colour_distiller(palette = "RdYlBu", guide = "colorbar")

					# Print it:
					message(paste("Printing txt",s,v, sep="_"))
					write.table(ptsInPol, paste("climato_",s,"_Levin_",v,"_2001-2099_raw.txt", sep=""), sep=" ", row.names=FALSE)
					gc()
					
					rm('ddf')


			} #eo v loop

	 
} #eo s loop






###### 3°) Re-load the outputs you have just extracted and re-interpolate them on Levin's grid through getMeans() function + while loop

#setwd("/Users/ben_fabio/Desktop/data_montpell/NEMOMED8_outputs/")
ddf <- read.table("climato_EH2.1_Levin_MLD_2001-2099_raw.txt", sep = " ", h=T)
dim(ddf)
summary(ddf)
colnames(ddf) ### Check the columns that meet your target periods (20-50 & 68-98): 
#	- DAE9 = 243:602  &  819:1190
#   - other scenarios = 231:590  &  808:1179
ggplot(ddf) + geom_point(aes(x=x, y=y, colour= Dec_1999)) + theme_bw() + coord_quickmap()


###### 4°) Re-interpolate NEMOMED8 climatologies on Levin's grid:
scenario <- c("EH2.1","EH3","EH3.1","EH3.2","EH4.1")
var <- c("SSS", "MLD", "SST")

# Load Levin's grid...
Layer <- ogrListLayers("grille/grille_med.shp")
grilleMed <- readOGR("grille/grille_med.shp", Layer)



for(s in scenario) {
	
		for(v in var) {

				# Reload NEMOMED8 climatology :
				message(paste("Reading",s,v, sep="_") )
				ptsInPol <- read.table(paste("climato_",s,"_Levin_",v,"_2001-2099_raw.txt", sep=""), header= TRUE, sep=" ")
				#class(ptsInPol)
				#dim(ptsInPol)
				#head(ptsInPol)

				# On crée un fichier texte avec les coordonnées et le sernum de la grille de Noam
				coords_Levin <- cbind("sernum"=grilleMed@data[,1], "X"=getSpPPolygonsLabptSlots(grilleMed)[,1], "Y"=getSpPPolygonsLabptSlots(grilleMed)[,2])
				#class(coords_Levin)
				gc()
				#write.table(coords_Levin, "2_processed/coords_Levin.txt", sep=" ", row.names=FALSE)


				#### MEAN SUR LES 8 CELLULES PLUS PROCHES, ITERATIVEMENT ####
				# On extrait les centroides des polygones de grilleMed pour faire un plus proche voisin sur 3 voisins####
				centroids <- getSpPPolygonsLabptSlots(grilleMed)
				neighbours <- knearneigh(centroids, k=8)
				#class(neighbours)
				#str(neighbours)

				centroids <- data.frame(centroids)
				names(centroids) <- c("x","y")

				neighbours <- data.frame(neighbours$nn)
				#head(neighbours)
				names(neighbours) <- seq(1:ncol(neighbours))


				### Separate periods before neighbour averaging
				varT3 <- ptsInPol[,c(808:1179)] 
				#colnames(varT3)
				#dim(varT3)

				getMeans <- function(i, k=8) {
  			  			vect <- rep(NA, ncol(varT3))
 					   	for (j in 1:length(vect)) {
   						 	ind <- as.numeric(neighbours[i,1:k])
   					 		vect[j] <- mean(c(varT3[ind,j]), na.rm=TRUE)
  						  }#eo for
  						vect
				} #eo fun
				

				# Apply fun in while loop
				message(paste("Re-interpolating",s,v, sep="_") )
				
				while(length(which(is.na(varT3[,1]) )) > 0) {
  				  		# On moyenne les valeurs sur les plus proches voisins
 			   			indices <- which(is.na(varT3[,1]))
  					  	temp <- mclapply(as.list(indices), getMeans, mc.cores = 4)
  					  	temp <- do.call(rbind.data.frame, temp)
 					   	names(temp) <- names(varT3)
  					  	varT3[indices,] <- temp
  					  	rm(temp)
  					  	gc()
				}# eo while loop


				### 5°) Save extracted climatologies:
				#class(varT3)
				varT3$x <- ptsInPol$x
				varT3$y <- ptsInPol$y
				#head(varT3)
				#summary(varT3) # No NAs, excellent
				rm('ptsInPol')
				
				message(paste("Printing",s,v, sep="_") )
				write.table(varT3, paste("climato_",s,"_Levin_",v,"_68-98.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
				gc()

		} # eo for v
		
} #eo for s




# Test:
ggplot(varT3) + geom_point(aes(x=x, y=y, colour= Nov_2049)) + theme_bw() + coord_quickmap() + scale_colour_distiller(palette = "YlOrRd", guide = "colorbar")




###### 4°) Aggregate the ''new'' NEMOMED8 fields within decadal time periods (65-74, 75-84, 85-94) and compute seasonal climatologies 

#setwd("/Users/ben_fabio/Desktop/data_montpell/NEMOMED8_outputs/_process2_Jan2015/HIS")

scenario <- c("EH1.7")
var <- c("SSS", "MLD", "SST")
decade <- c("65-74","75-84","85-94")

### For testing:
#s <- "EH1.7"
#v <- "SST"
#d <- "65-74"


for(s in scenario) {


		for(v in var) {

			 message(paste("Reading",s,v, sep="_") )
			 ptsInPol <- read.table(paste("climato_",s,"_Levin_",v,"_1965-1994.txt", sep=""), header= TRUE, sep="\t")

				
			 for(d in decade) {
				
			      ### Separate each decade :
			 	  message(paste("Choosing decades",s,v,d, sep="_") )
				  
				  if( d == "65-74" ) { t <- ptsInPol[,c(1:120, 361:362)] }  else if (d == "75-84") { t <- ptsInPol[,c(121:240, 361:362)] } else if (d == "85-94") { t <- ptsInPol[,c(241:362)] }
				  gc()
				  

				  ### Compute seasonal climatologies from these using grep() function with multiple pattern
				  message(paste("Computing seasonal climatologies",s,v,d, sep="_") )
				  
				  dates <- colnames(t)
				  t_spring <- as.matrix(t[, grep("Mar|Apr|May", x= dates)])
				  t_summer <- as.matrix(t[, grep("Jun|Jul|Aug", x= dates)])
				  t_fall <- as.matrix(t[, grep("Sep|Oct|Nov", x= dates)])
				  t_winter <- as.matrix(t[, grep("Dec|Jan|Feb", x= dates)])
				  
				  ### rowMeans() to compute seasonal climatologies:
				  t_sp <- rowMeans(t_spring[,c(1:30)], na.rm = F)
				  t_su <- rowMeans(t_summer[,c(1:30)], na.rm = F)
				  t_f <- rowMeans(t_fall[,c(1:30)], na.rm = F)
				  t_w <- rowMeans(t_winter[,c(1:30)], na.rm = F)
				  gc()

				  ddf <- data.frame(Spring = t_sp, Summer = t_su, Fall = t_f, Winter = t_w, x = ptsInPol$x, y = ptsInPol$y)
				  #ggplot(ddf) + geom_point(aes(x=x, y=y, colour= Summer)) + theme_bw() + coord_quickmap() + scale_colour_distiller(palette = "YlOrRd", guide = "colorbar")
				  ### OK
				  message(paste("Printing seasonal climatologies",s,v,d, sep="_") )
				  write.table(ddf, paste("climato_seasons_",s,"_Levin_",v,"_",d,".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
				  	  
				  
			}# eo d
						  
		}# eo v
						  
}# eo s
				
				  
### Check:
ddf <- read.table("climato_seasons_EH1.7_Levin_MLD_65-74.txt", h=T,sep="\t")				  
dim(ddf)				  
summary(ddf)		  
ggplot(ddf) + geom_point(aes(x=x, y=y, colour= Spring)) + theme_bw() + coord_quickmap() + scale_colour_distiller(palette = "YlOrRd", guide = "colorbar")			
# Good


			
			
###### For future (2020-2050 & 2069-2098) projections now: 			
			
scenario <- c("DE9")
#"EH2.1","EH3","EH3.1","EH3.2","EH4.1")
var <- c("SSS", "MLD", "SST")
decade <- c("69-78","79-88","89-98")

### For testing:
# s <- "DE9"
# v <- "SST"
# d <- "69-78"


for(s in scenario) {


		for(v in var) {

			 message(paste("Reading",s,v, sep="_") )
			 ptsInPol <- read.table(paste("climato_",s,"_Levin_",v,"_68-98.txt", sep=""), header= TRUE, sep="\t")
			 # dim(ptsInPol)
			 # colnames(ptsInPol)
				
			 for(d in decade) {
				
			      ### Separate each decade :
			 	  message(paste("Choosing decades",s,v,d, sep="_") )
				  
				  if( d == "69-78" ) { t <- ptsInPol[,c(13:132, 373:374)] }  else if (d == "79-88") { t <- ptsInPol[,c(133:252, 373:374)] } else if (d == "89-98") { t <- ptsInPol[,c(253:374)] }
				  gc()
				  

				  ### Compute seasonal climatologies from these using grep() function with multiple pattern
				  message(paste("Computing seasonal climatologies",s,v,d, sep="_") )
				  
				  dates <- colnames(t)
				  t_spring <- as.matrix(t[, grep("Mar|Apr|May", x= dates)])
				  t_summer <- as.matrix(t[, grep("Jun|Jul|Aug", x= dates)])
				  t_fall <- as.matrix(t[, grep("Sep|Oct|Nov", x= dates)])
				  t_winter <- as.matrix(t[, grep("Dec|Jan|Feb", x= dates)])
				  
				  ### rowMeans() to compute seasonal climatologies:
				  t_sp <- rowMeans(t_spring[,c(1:30)], na.rm = F)
				  t_su <- rowMeans(t_summer[,c(1:30)], na.rm = F)
				  t_f <- rowMeans(t_fall[,c(1:30)], na.rm = F)
				  t_w <- rowMeans(t_winter[,c(1:30)], na.rm = F)
				  gc()

				  ddf <- data.frame(Spring = t_sp, Summer = t_su, Fall = t_f, Winter = t_w, x = ptsInPol$x, y = ptsInPol$y)
				  # ggplot(ddf) + geom_point(aes(x=x, y=y, colour= Summer)) + theme_bw() + coord_quickmap() + scale_colour_distiller(palette = "YlOrRd", guide = "colorbar")
				  ### OK
				  message(paste("Printing seasonal climatologies",s,v,d, sep="_") )
				  write.table(ddf, paste("climato_seasons_",s,"_Levin_",v,"_",d,".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
				  	  
				  
			}# eo d
						  
		}# eo v
						  
}# eo s	
					
			
### Check:
ddf <- read.table("climato_seasons_EH3_Levin_MLD_69-78.txt", h=T,sep="\t")				  
dim(ddf)				  
summary(ddf)		  
ggplot(ddf) + geom_point(aes(x=x, y=y, colour= Winter)) + theme_bw() + coord_quickmap() + scale_colour_distiller(palette = "YlOrRd", guide = "colorbar")			
# Good		
			
			
			
###### 5°) From the seasonal climatologies of each decade, compute decadal average of : SST, SSS, MLD and stdSST.
scenario <- c("EH2.1","EH3","EH3.1","EH3.2","EH4.1","DE9")
decade <- c("69-78","79-88","89-98")

### For testing:
# s <- "DA9"
# d <- "65-74"
 
 
for(s in scenario) { 

		for(d in decade) {

					message(paste("Reading & averaging variables",s,v,d, sep="_") )
					library(matrixStats)
					sst <- read.table(paste("climato_seasons_",s,"_Levin_SST_",d,".txt", sep=""), header= TRUE, sep="\t")
					sst <- as.matrix(sst)
					avgSST <- rowMeans(sst[,1:4])
					stdSST <- rowSds(sst[,1:4])
					rm('sst')

					sss <- read.table(paste("climato_seasons_",s,"_Levin_SSS_",d,".txt", sep=""), header= TRUE, sep="\t")
					sss <- as.matrix(sss)
					avgSSS <- rowMeans(sss[,1:4])
					rm('sss')

					mld <- read.table(paste("climato_seasons_",s,"_Levin_MLD_",d,".txt", sep=""), header= TRUE, sep="\t")
					mld <- as.matrix(mld)
					avgMLD <- rowMeans(mld[,1:4])


					clim <- data.frame(avgSST = avgSST, stdSST = stdSST, avgSSS = avgSSS, avgMLD = avgMLD, x= mld[,"x"], y= mld[,"y"])
					rm('mld','avgSST','stdSST','avgSSS','avgMLD')
					gc()

					# Check:
					# ggplot(clim) + geom_point(aes(x=x,y=y,colour=avgMLD)) + theme_bw() + coord_quickmap() + scale_colour_distiller(palette = "YlOrRd", guide = "colorbar")	
					# Good !
					
  					message(paste("Printing seasonal climatologies",s,v,d, sep="_") )
  				  	write.table(clim, paste("climato_decadal_",s,"_Levin_",d,".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
					
					gc()
					
		}# eo d

}# eo s


### Check :
clim <- read.table("climato_decadal_EH4.1_Levin_40-49.txt", header= TRUE, sep="\t")
dim(clim)
summary(clim)
ggplot(clim) + geom_point(aes(x=x,y=y,colour=stdSST)) + theme_bw() + coord_quickmap() + scale_colour_distiller(palette = "YlOrRd", guide = "colorbar")	







			
			
			