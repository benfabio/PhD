
##### 05/03/2015 - MARBEC - Fabio Benedetti, Tarek Hattab, François Guilhaumon
##### Script for : 
#					- pseudo-absences simumations for all species based on an environmentally- and geographically-weighted method

### Last update : 08/06/2015

### Based on : Hattab & al. (2014) ; Hengl & al. (2009) ; Engler & al. (2004)

# ---------------------------------------------------------------------------------------------------------------------------------

library("rgdal")
#library("raster")
library("maptools")
library("ade4")
library("Hmisc")
library("adehabitat")
library("PresenceAbsence")
library("dismo")
library("SDMTools")
library("spatstat")
library("ggplot2")
library("stringr")
library("parallel")
library("raster")

# ---------------------------------------------------------------------------------------------------------------------------------

### List all species' presences files
setwd("/home/fabio/presence_raster_files/")
species_files <- dir()[grep(".grd", dir())]
species_files <- species_files[-30] # Get rid of clims_woa13_rasters
length(species_files) ; gc()

WD <- getwd() # base WD

### For testing function:
s <- "Temora_stylifera.grd"

### Define function that generates pseudo-absence with environmental and geographical weighting (Script16.1):
# List of runs : 
t <- c(1)

pseudoAbsGenerator <- function(s) {

		for(t in times) {
					
					ss <- str_replace(s, ".grd", "")
					message(paste("Starting simulations of pseudoAbs for ",ss,"_run_",t, sep=""))
					
					# Turn vectors to several objects of class asc:
					ras_sst <- raster("clims_woa13_baseline_04.grd", band = 3)
					sst_asc <- asc.from.raster(ras_sst)

					ras_stdsst <- raster("clims_woa13_baseline_04.grd", band = 4)
					stdsst_asc <- asc.from.raster(ras_stdsst)

					ras_sss <- raster("clims_woa13_baseline_04.grd", band = 5)
					sss_asc <- asc.from.raster(ras_sss)

					# ras_mld <- raster("clims_woa13_baseline_04.grd", band = 6)
					# mld_asc <- asc.from.raster(ras_mld)
					
					#dim(na.omit(as.data.frame(ras_sst, xy = TRUE)))
					#head(na.omit(as.data.frame(ras_sst, xy = TRUE)))

					rm(ras_sst, ras_sss, ras_stdsst)

					# Coerce into list: 
					asc.list <- list(sst_asc, stdsst_asc, sss_asc)
					# Lacking variables' names
					names(asc.list)[1] <- "avgSST"
					names(asc.list)[2] <- "stdSST"
					names(asc.list)[3] <- "avgSSS"

					# Transform asc objects into objects of class kasc
					predictors_enfa <- as.kasc(asc.list) # you-hou
					rm(sst_asc, stdsst_asc, sss_asc)

					### Load presence coordinates (X,Y) of the species (dataframe)
					setwd("/home/fabio/presence_raster_files/")
					ras_sp <- raster(s)
					sp1 <- as.data.frame(ras_sp, xy = TRUE)
					sp1[is.na(sp1$layer), "layer"] <- 0
					### Check nb of presences : 
					# nrow(sp1[which(sp1$layer == 1),])
					# locs <- sp1[which(sp1$layer == 1), c("x","y")]
					# head(locs) ; summary(locs) 
					
					### Filter possible land occurrences with mask : 
					coords_glob <- as.data.frame(ras <- brick("clims_woa13_baseline_04.grd"), xy = TRUE)
					# head(coords_glob)
					coords_glob <- coords_glob[,c(1:2,5)] # !! Loaded coordinates twice !!
					coords_glob$land <- ifelse(is.na(coords_glob$avgSST), 1, 0)
					coords_glob$land2 <- ifelse(coords_glob$land == 1, TRUE, FALSE) # OK
					# Filter presences : 
					sp1[coords_glob$land2, c("layer")] <- 0
					### Re-check nb of presences : 
					# nrow(sp1[which(sp1$layer == 1),])
					locs <- sp1[which(sp1$layer == 1), c("x","y")]
					# dim(locs)
					gc()


					### ENFA : environmental weighting
					# Prepare data for enfa : data2enfa()
					dataenfa1 <- data2enfa(predictors_enfa, locs) 
					#gc()
					# str(dataenfa1)

					# Perform ENFA
					pc <- dudi.pca(dataenfa1$tab, scannf = FALSE, center = TRUE, scale = TRUE)
					enfa <- enfa(pc, dataenfa1$pr, scannf = FALSE)

					pred <- predict(enfa, dataenfa1$index, dataenfa1$attr)
					pred <- as.SpatialGridDataFrame.im(asc2im(pred))

					# Ranking precitions from ENFA in 'rankv' of pred
					slot <- slot(pred, "data")
					# str(slot) # data.frame with 2 columns : v & rankv : 

					pred$rankv <- rank(slot$v, ties.method= "first")
					sum.dist <- summary(pred$v) # retrieves median, mean, max etc.
					# sum.dist[["Median"]]

					### Compute habitat suitability index (HSI) according to the ENFA weighting
					pred$hsit <- ifelse(pred$v < sum.dist[["Median"]], (1 - pred$rankv/max(pred$rankv))*100, (1-pred$rankv/max(pred$rankv))*100)
					pred$hsi <- 100*round((pred$hsit - min(pred$hsit, na.rm=T))/(max(pred$hsit, na.rm=T) - min(pred$hsit, na.rm=T)), 3)
					#  summary(pred)
					rm(dataenfa1, enfa, sum.dist) 
					gc()

					### Buffer map : geographical weighting additionally to the environmental weighting
					sppp <- as.ppp(data.frame(locs$x, locs$y), c(-179.875, 179.875,-89.875, 89.875)) # species'presences coordinates and extent of modeling
					# str(sppp)
					# summary(sppp)

					buffer <- distmap(sppp) ### distance between presences'locations and geographical cells
					buffer <- as.SpatialGridDataFrame.im(buffer)
					buffer <- raster(buffer)
					# plot(buffer)

					# Create a mask via one predictor variable 's raster
					ras_sst <- raster("clims_woa13_baseline_04.grd", band = 3)
					buffer <- raster::resample(buffer, ras_sst, method= "bilinear")
					# Resample transfers values between non matching Raster* objects (in terms of origin and resolution).
					# Use projectRaster if the target has a different coordinate reference system (projection).
					# class(buffer)
					# str(buffer)
					# plot(buffer)
					### Check out buffer's coordinates values:
					# coords_buffer <- as.data.frame(buffer, xy = TRUE) ; head(coords_buffer)
					# rm(coords_buffer)
					## OK  !

					mask <- is.na(ras_sst)  # all cells with NA values 
					mask <- as(mask, "SpatialPixelsDataFrame")
					liste_na <- (mask@data$layer == 1) # TRUE/FALSE vector 
					# (length = global cells), with TRUE for cells which had NA values in sst raster (i.e. non shelf)
					# unique(liste_na)
					# length(liste_na)
					rm(ras_sst)
					gc()

					### Normalize the distance from the observed occurrences by the maximum distance
					buffer <- (buffer/buffer@data@max) * 100
					# str(buffer)
					# str(buffer@data@values)
					# summary(buffer@data@values)
					if( sum(is.na(buffer@data@values)) > 0 ) {
						buffer@data@values[is.na(buffer@data@values)] <- min(na.omit(buffer@data@values))
					}
					
					buffer <- as(buffer, "SpatialPixelsDataFrame")
					# str(buffer)
					# buffer@data[liste_na,][1:1000]
					
					buffer@data$v <- replace(x= c(buffer@data$v), list= liste_na, value = NA)
					# ‘replace’ replaces the values in ‘x’ with indices given in ‘list’
     			   	# by those given in ‘values’. If necessary, the values in ‘values’ are recycled.
					pred$dist <- buffer@data
					pred$weight <- ( (pred$dist + (100-pred$hsi)) / 2)^2  
					### Eq.(11) of Hengl & al. (2009) : 
					### Where square term is made to ensure that pseudo-absences are progressively selected at the edge of low HSI values.
					### This will make the pseudo-absences follow a Poisson density distribution.
					### Hengl & al. (2009) extent the method by considering the locations of occurrence points in geographical space too !
					### Hence the eq. above, where 'pred$dist' is the normalized distance in the range (from 0 to 100%) i.e. the distance
					### from the observation points divided by the maximum distance (definition of the 'buffer' above).
					# str(pred)
					# summary(pred@data)

					dens.weight <- as.im(as.image.SpatialGridDataFrame(pred[6]))  # pred[6] because contains weights !
					# converts to pixel image

					rm(buffer, sppp, liste_na, pc)
					enfa <- pred["hsi"]
					rm(pred, enfa) 
					gc()

					### Pseudo-absence simulations : 
					# Pseudo_absence_simulation
					# ?spatstat::rpoint
					# Generates a random point pattern containing n independent, identically distributed random points with any specified distribution.
					# This function generates ‘n’ independent, identically distributed random points with common probability density proportional to ‘f’ argument.

					l <- nrow(locs)*50
					pseudo <- rpoint(l, f = dens.weight) ### Generates as many pseudo-absences as presences, 
								     			 		 ### according to an envrionemntally- and geographically-weighted probability function !
					# rm(dens.weight)

					# Create the dataframe 
					pseudo.absences <- data.frame(x = pseudo$x, y = pseudo$y, p = rep(0,l) )
					# coordinates(pseudo.absences) <- ~x+y
					# proj4string(pseudo.absences) <- proj4string(ras_sp)
					
					### !¡! Issue with coordinates, not projected correctly on WOA13's grid for some reason !¡!
					### Need to re-project pseudo-absences on correct 
					# ras_sst <- raster("clims_woa13_baseline_04.grd", band = 3)
					
					true.coords <- coords_glob[,c("x","y")]
					your.coords <- pseudo.absences[,c("x","y")]
					your.coords$ids <- paste(pseudo.absences$x, pseudo.absences$y, sep="_") # Add some IDS for lapplying

					### Retreive all true 'x' and 'y' values that correspond to the closest values of absences's coordinates (combining which() and abs())
					new.coords <- lapply(your.coords$ids, function(id) {
									# First, find the proper value of 'x' among 'true.coords$x'
									new.x <- which(abs(true.coords$x - your.coords[which(your.coords$ids == id),"x"]) == min(abs(true.coords$x - your.coords[which(your.coords$ids == id),"x"])))
									new.x <- unique(true.coords[new.x,"x"])
				
									# Now, find the proper value of 'y' among 'true.coords$y'
									new.y <- which(abs(true.coords$y - your.coords[which(your.coords$ids == id),"y"]) == min(abs(true.coords$y - your.coords[which(your.coords$ids == id),"y"])))
									new.y <- unique(true.coords[new.y,"y"])
									# Return
									new.coords <- data.frame(x = new.x, y = new.y)
									return(new.coords)
					}) #eo lapply
					new.coords <- do.call(rbind, new.coords)
					
					# dim(new.coords)
					# head(new.coords)
					# dim(your.coords)
					# head(your.coords)

					# Replace DataSpecies' absences coordinates by new.coords  
					pseudo.absences[,c("x","y")] <- new.coords
					
					# Ok, check with mask if there are any land points...
					# pseudo.absences[coords_glob$land2, c("p")] <- 9999
					### Seems ok !
					
					# XXX <- rasterize(pseudo.absences[,1:2], ras_sp, field = 9999, background = NA)
					# coords <- as.data.frame(XXX, xy = TRUE) ; head(coords)
					# dim(coords[which(coords$layer == 9999),])
					# not OK
					
					# pseudo.absences <- coords[which(coords$layer == 1),]
					# coordinates(pseudo.absences) <- ~x+y
					# proj4string(pseudo.absences) <- proj4string(ras_sp)

					# Combine species' actual presence with the newly generated peudo-absences ! 
					locs$p <- 1
					# Species <- SpatialPointsDataFrame(coordinates(locs[,1:2]), data = locs)
					# proj4string(Species) <- proj4string(ras_sp)
					# str(Species)
					# str(pseudo.absences)
					occurence.all <- rbind(locs, pseudo.absences)
					
					# Check for duplicates :
					duplicates <- duplicated(occurence.all[,c("x","y")])
					occurence.all[duplicates,] <- NA
					occurence.final <- na.omit(occurence.all)
					# duplicated(occurence.all[,c("x","y")])
					# duplicated(occurence.final[,c("x","y")])
					### Ok ! 
			
					# Save as table :
			   		# Create directory for each species where you will store each of the 10 pseudoAbs generations
		   			if(t == 1) { dir.create(ss)
						         setwd(paste(WD,"/",ss,"/", sep="")) 
					   } else { 
							setwd(paste(WD,"/",ss,"/", sep="")) 
					} #eo ifelse
					
					# Save
					message(paste("Saving",ss,"run",t, sep="_"))
					write.table(occurence.final,paste("pabs_",ss,"_run",t,".txt", sep = ""), sep="\t")
					# Make room for next
					rm(occurence.all, locs, pseudo.absences, l, pseudo, asc.list, predictors_enfa, ras_sp, occurence.final, true.coords, your.coords, 
					duplicates) 
					gc()
		
					setwd(WD)	
			
			} #eo t
			
} #eo function


### apply function 10 times on species_files :
#species_files <- "Parapontella_brevicornis.grd"
mclapply(X = species_files, FUN = pseudoAbsGenerator, mc.cores = 12)


