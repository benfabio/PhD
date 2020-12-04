
##### Script for computing the species' changes in occupancy patterns in the Mediterranean Sea due CC

# - based on the revised communities submitted to Ecography back in July.
# - based on 6 ENMs (chosen according to Fig. 6 of the article), 3 prevalences (may be turned don to 2), 3 forcings (SRES) and the 3 CV runs


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
library("doParallel")
library("ggplot2")
library("RColorBrewer")

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
# p <- "n1"
# sdm <- "SRE"
# scenar <- "A2ARF"
# R <- "RUN1"

species_ranger_2 <- function(sdm = SDMs) {
							
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
									commT1 <- get(load(paste("Community_table_T2_run_",p,"_1_",sdm,"_",scenar,"_",R,".Rdata", sep="")))
									
									# Need to re-supply cells' ids
									rownames(commT0) <- paste(coords$Lon, coords$Lat, sep= "_")
									rownames(commT1) <- paste(coords$Lon, coords$Lat, sep= "_")
									ids <- rownames(commT1)
									
									### Make clever indices to quickly identify the spp. that left/arrived/stayed within the cell after CC
									commT1 <- commT1 * 2
									commsum <- commT0 + commT1
									# Use them in apply to find sum of 1,2,3 per species
									
									resChanges <- apply(commsum, 2, function(l) {
									
											# sum of cells between T0 and T2 == 3 (where not affected by CC in the cell)
											nCommons <- sum(l == 3)
											# sum of cells the species left because of CC == 1
											nLeft <- sum(l == 1)
											# sum of cells the species entered thx to of CC == 1
											nArrived <- sum(l == 2)
											# sum of cells the species was always absent 
											nAbsent <- sum(l == 0)
										
											return(data.frame(commons = nCommons, left = nLeft, arrived = nArrived, absent = nAbsent,
											prevalence = p, ENM = sdm, forcing = scenar, cv_run = R))
											
									}) # eo apply
									
									ddf <- do.call(rbind, resChanges)
									rownames(ddf) <- gsub(".", "_", rownames(ddf), fixed= T)
									
									setwd(paste(WD,"/","Species_changes_T2","/", sep=""))
									save(ddf, file = paste("ddf_",p,"_",scenar,"_",sdm,"_",R,".Rdata", sep= "") )
										
									# Make room
									rm(commT0, commT1, resChanges, ddf)
										
								} # eo R in RUNS
									
							} # eo p in prevalences

						} # eo scenar in scenarios

} # eo FUN		
								
library("parallel")
mclapply(X = SDMs, FUN = species_ranger_2, mc.cores = 10)
gc()




### Ok, let's examine the results	
setwd(paste(WD,"/","Species_changes_T2","/", sep= ""))
length(dir())
files <- dir()

res <- lapply(files, function(f) {
			d <- get(load(f))
			return(d)
} ) # eo lapply

ddf <- do.call(rbind, res) ; rm(res)
#dim(ddf)
ddf$sp.name <- gsub('[0-9]+', '', rownames(ddf))

library("dplyr")
table <- data.frame(ddf %>%
  group_by(sp.name) %>%
  summarise(avg_cells_common= mean(commons), avg_cells_left= mean(left), 
  			avg_cells_arrived= mean(arrived), avg_cells_absent= mean(absent)))
			
table$ratio <- (table$avg_cells_left / table$avg_cells_arrived)			

#table[order(table$ratio, decreasing = T),]
setwd(WD)
save(ddf, file= paste("species_changes_T2.Rdata", sep= ""))




