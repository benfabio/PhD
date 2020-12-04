


##### 15/04/2016 - R CMD BATCH SCRIPT_species_shifte.R
##### Script for : 
#		- extracting the previously-created communities for T0 and T2 in the MED for the 106 spp.
#		- compute species' range centroids shifts
#		- examine species that experience the most changes, link it with prevalence patterns (Script24.1)
#		- find a nice way to plot centroid distribution shifts


# ------------------------------------------------------------------------------------------------------------------------------------

library("stringr")
library("reshape2")
library("sp")
library("matrixStats")
library("RColorBrewer")
library("scales")
library("fields")
library("dplyr")
library("parallel")
library("ggplot2")
library("RColorBrewer")
library("geosphere")
library("rgdal")
library("raster")
library("geosphere")

# ------------------------------------------------------------------------------------------------------------------------------------

### Basic stuff to start with
WD <- getwd()
setwd("./niche_modelling")
second.wd <- getwd()
setwd(second.wd)
sp.names <- dir() ; length(sp.names)
sp.names <- sp.names[-c(6:10,48,56,65,89,116)]
setwd(WD) ; length(sp.names)

# List of pseudo-absences runs
#runs <- as.character(c(1:3))
# List of SDMs
SDMs <- c("SRE","RF","MARS","GLM","ANN","MAXENT")
# List of scenarios 
scenarios <- c("A2ARF","A1BARF","B1ARF")
# List pf prevalences
preval <- c("n1","n10","n50")
# Cross-evaluation runs
RUNS <- c("RUN1","RUN2","RUN3")


# ------------------------------------------------------------------------------------------------------------------------------------

### Adopt the following approach, for each species create a multidimensional array with dimensions for : 
#	- number of cells : 3 values (T0, T2, ∆)
#	- pseudo-absence run : 10 values
#	- SDM : 10 values
#	- scenario : 6 values
#	- community number : 1:100

### Community tables: species' P/A table in the Med Sea. Just lacking coordinates...let's get 'em back
setwd(paste(WD,"/","ANOVA_tables","/", sep=""))
#  dir()
numbers <- c(1:10)
res <- lapply(numbers, function(n) {
			d <- get(load(paste("avg_all_dSR_T2_subtable_",n,".Rdata", sep= "")))
			ids <- data.frame(unique(d$id))
			colnames(ids) <- "id"
			return(ids)
} ) # eo lapply

ids <- do.call(rbind, res) ; rm(res)
# head(ids)
# unique(ids)
### Retrieve coordinates
coords <- data.frame(str_split_fixed(ids$id, "_", 2))
colnames(coords)[1:2] <- c("Lon","Lat")
coords$Lon <- as.numeric(levels(coords$Lon))[coords$Lon]
coords$Lat <- as.numeric(levels(coords$Lat))[coords$Lat]
rm(ids)

# ------------------------------------------------------------------------------------------------------------------------------------

# For testing:
 p <- "n1"
 sdm <- "SRE"
 scenar <- "A2ARF"
 R <- "RUN1"


centroid_shifter <- function(sdm = SDMs) {
					
						# For each SDMs
						for(scenar in scenarios) {
	
							# For RCM bounday forcing condition 
							for(p in preval) {

								# For each random community 
								for(R in RUNS) {
									
										# Useless message
										message( paste("Doing ",sdm,"_",scenar,"_",p,"_",R,".Rdata", sep= "") )
											
										# Loading community tables
										setwd(paste(WD,"/","Communities_10_05_16","/", sep= ""))
										#Community_table_T2_run_n50_1_GAM_B1ARF_RUN1.Rdata
										commT0 <- get(load(paste("Community_table_T0_run_",p,"_1_",sdm,"_",scenar,"_",R,".Rdata", sep="")))
										commT2 <- get(load(paste("Community_table_T2_run_",p,"_1_",sdm,"_",scenar,"_",R,".Rdata", sep="")))
										# Need to re-supply cells' ids
										rownames(commT0) <- paste(coords$Lon, coords$Lat, sep= "_")
										rownames(commT2) <- paste(coords$Lon, coords$Lat, sep= "_")
										
										# Computing species range centroid at T0 and T2: average lon + average lat
										centroidsT0 <- apply(commT0, 2, function(l) {
													# Get ids of the presences
													names(l) <- rownames(commT0)
													pres_id <- names(l)[l==1]
													pres_coord <- data.frame(str_split_fixed(pres_id, "_", 2))
													colnames(pres_coord)[1:2] <- c("x","y")
													avg_x <- mean(as.numeric(levels(pres_coord$x) )) 
													avg_y <- mean(as.numeric(levels(pres_coord$y) ))
													
													return(data.frame(mean_lon_T1 = avg_x, mean_lat_T1 = avg_y, prevalence = p, 
													ENM = sdm, forcing = scenar, cv_run = R))
										}) # eo apply
										
										# Rbind
										centroidsT0 <- do.call(rbind, centroidsT0)
										
										centroidsT2 <- apply(commT2, 2, function(l) {
													# Get ids of the presences
													names(l) <- rownames(commT2)
													pres_id <- names(l)[l==1]
													pres_coord <- data.frame(str_split_fixed(pres_id, "_", 2))
													colnames(pres_coord)[1:2] <- c("x","y")
													avg_x <- mean(as.numeric(levels(pres_coord$x) )) 
													avg_y <- mean(as.numeric(levels(pres_coord$y) ))
										
													return(data.frame(mean_lon_T2 = avg_x, mean_lat_T2 = avg_y, prevalence = p, 
													ENM = sdm, forcing = scenar, cv_run = R))
										}) # eo apply
										
										# Rbind
										centroidsT2 <- do.call(rbind, centroidsT2)
										
										### Ok, now that you have both centroids (future and present), need to compute distance between coordinates
										centroids <- cbind(centroidsT0, centroidsT2)
										centroids <- centroids[,-c(9:12)] # drop useless columns
										
										# Use mutate() to add columns
										centroids$distm <- NA
										require("geosphere")
										
										for(i in 1:nrow(centroids)) {
											centroids[i,"distm"] <- distm(centroids[i,c("mean_lon_T1", "mean_lat_T1")], centroids[i,c("mean_lon_T2", "mean_lat_T2")], fun = distHaversine)
										} # eo dist for loop
										# Convert to kilometers
										centroids$distm <- (centroids$distm)/1000
											
										### Now, need to tompute the angle so we will be able to plot future range centroid from a (0,0) point
										# Need to compute latitudinal distance 
										centroids$dist_lat <- NA
										for(i in 1:nrow(centroids)) {
											centroids[i,"dist_lat"] <- distm(centroids[i,c("mean_lon_T1", "mean_lat_T1")], centroids[i,c("mean_lon_T1", "mean_lat_T2")], fun = distHaversine)
										} # eo dist for loop
										# Convert to kilometers
										centroids$dist_lat <- (centroids$dist_lat)/1000
										
										# Need to compute longitudinal distance 
										centroids$dist_lon <- NA
										for(i in 1:nrow(centroids)) {
											centroids[i,"dist_lon"] <- distm(centroids[i,c("mean_lon_T1", "mean_lat_T1")], centroids[i,c("mean_lon_T2", "mean_lat_T1")], fun = distHaversine)
										} # eo dist for loop
										# Convert to kilometers
										centroids$dist_lon <- (centroids$dist_lon)/1000
											
										### Useless message (again)
										#message(paste("Saving list for run ",r, sep=""))
										setwd(paste(WD,"/","Species_shift_T2","/", sep=""))
										save(centroids, file = paste("ddf_",p,"_",scenar,"_",sdm,"_",R,".Rdata", sep=""))
											
										# Make room
										rm(commT0, commT2, resChanges, centroids)
											
								} # eo communities		
									
						} # eo scenarios
							
			} # eo sdm in SDMs
		
} # eo FUN


### Apply function within mclapply
mclapply(SDMs, centroid_shifter, mc.cores = 10)
gc()

### 19/07/16 - Check results and make plots
setwd(paste(WD,"/","Species_shift_T2","/", sep= ""))
files <- dir() ; files

res <- lapply(files, function(f) {
			d <- get(load(f))
			return(d)
} ) # eo lapply

ddf <- do.call(rbind, res) ; rm(res)
dim(ddf)
head(ddf)

ddf$sp.name <- gsub('[0-9]+', '', rownames(ddf))

library("dplyr")
table <- data.frame(ddf[which(ddf$ENM != "SRE"),] %>%
  group_by(sp.name) %>%
  summarise(mean_lon_T0 = mean(mean_lon_T1, na.rm = T), mean_lat_T0 = mean(mean_lat_T1, na.rm = T), 
  mean_lon_T2 = mean(mean_lon_T2, na.rm = T), mean_lat_T2 = mean(mean_lat_T2, na.rm = T),
  avg_distm = mean(distm, na.rm = T)) )

#quartz()
plot <- ggplot(table) + 
		geom_point(aes(x = mean_lon_T0, y = mean_lat_T0), shape = 21, colour = "black", fill = "blue", size = 2, alpha = 0.6) +
		geom_point(aes(x = mean_lon_T2, y = mean_lat_T2), shape = 21, colour = "black", fill = "red", size = 2, alpha = 0.6) +
		xlab("Mean longitude (°)") + ylab("Mean latitude (°)") + theme_bw()

ggsave("test_plot.pdf", plot = plot, width = 7, height = 7)


table <- data.frame(ddf[which(ddf$ENM == "SRE"),] %>%
  group_by(sp.name) %>%
  summarise(mean_lon_T0 = mean(mean_lon_T1, na.rm = T), mean_lat_T0 = mean(mean_lat_T1, na.rm = T), 
  mean_lon_T2 = mean(mean_lon_T2, na.rm = T), mean_lat_T2 = mean(mean_lat_T2, na.rm = T),
  avg_distm = mean(distm, na.rm = T)) )

#quartz()
plot <- ggplot(table) + 
		geom_point(aes(x = mean_lon_T0, y = mean_lat_T0), shape = 21, colour = "black", fill = "blue", size = 2, alpha = 0.6) +
		geom_point(aes(x = mean_lon_T2, y = mean_lat_T2), shape = 21, colour = "black", fill = "red", size = 2, alpha = 0.6) +
		xlab("Mean longitude (°)") + ylab("Mean latitude (°)") + theme_bw()

ggsave("test_plot2.pdf", plot = plot, width = 7, height = 7)


save(ddf, file= paste("species_shifts_T2.Rdata", sep= ""))