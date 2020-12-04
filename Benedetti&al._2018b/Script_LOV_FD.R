

##### 25/07/16 - Script for computing FD from the 'FD' package and estimate ∆FD from S. Villéger's indices (FRic, FEve and FDis)


# ------------------------------------------------------------------------------------------------------------------------------------

# libraries
library("stringr")
library("reshape2")
library("sp")
library("matrixStats")
library("scales")
library("fields")
library("dplyr")
library("parallel")
library("doParallel")
library("ggplot2")
library("RColorBrewer")

library("vegan")
library("geometry")
library("FD")

# ------------------------------------------------------------------------------------------------------------------------------------

### Basic stuff to start with
WD <- getwd()
setwd("./niche_modelling")
second.wd <- getwd()
setwd(second.wd)
sp.names <- dir() ; length(sp.names)
sp.names <- sp.names[-c(6:10,48,56,65,89,116)]
setwd(WD) ; length(sp.names)

### Modelling parameters
# List of SDMs
SDMs <- c("SRE","RF","MARS","GLM","ANN","MAXENT")
# List of scenarios 
scenarios <- c("A2ARF","A1BARF","B1ARF")
# List pf prevalences
preval <- c("n1","n10","n50")
# Cross-evaluation runs
RUNS <- c("RUN1","RUN2","RUN3")

### Read functional traits matrix:
setwd(WD)
traits <- read.csv("table_traits_biogeography.csv", h = T, dec = ",", sep = ";")
traits <- traits[,c(1,65:68)]
str(traits)
# Need to properly define the traits: as.numeric vs. as.ordered() vs. as.factor()
rownames(traits) <- traits$sp.name

### Give it a try without the randomization :
p <- "n50"
sdm <- "MAXENT"
scenar <- "A2ARF"
R <- "RUN1"

setwd(paste(WD,"/","Communities_10_05_16","/", sep= ""))
commT0 <- get(load( paste("Community_table_T0_run_",p,"_1_",sdm,"_",scenar,"_",R,".Rdata", sep= "") ))
commT2 <- get(load( paste("Community_table_T2_run_",p,"_1_",sdm,"_",scenar,"_",R,".Rdata", sep= "") ))

### Check consistency between species names:
# colnames(commT2) ; traits$sp.name
colnames(commT0) <- gsub(".", "_", colnames(commT0), fixed = TRUE)
colnames(commT2) <- gsub(".", "_", colnames(commT2), fixed = TRUE)

FD_T0 <- FD::dbFD(x = traits[,c(2:5)], a = commT0[26000:26200,], w.abun = F, calc.FRic = T, corr = "none")
str(FD_T0)
FD_T2 <- FD::dbFD(x = traits[,c(2:5)], a = commT2, w.abun = F, calc.FRic = T, corr = "none")
str(FD_T2)

# ----------------------------------------------------------------------------------------------------------------------------

##### Important notes and FUNCTIONS
### Second, prepare the function that will get you:
# - distribution of the 999 null estimates of ∆FD/∆PD
# - the values of observed ∆FD/∆PD
# - the p-value corresponding to the proportion of the 999 null estimates that are >= the observed ∆FD/∆PD
# - the Standardized Effect Size (SES), which can also be used to test the significance of values to a null model
### And store the p-values and the SES for each cell, and save a list of length == 26490 containing all the cells and the 2 stats


# Within each cell projection (~ ENM + SRES + FORCING + pseudo-abs run + simulated comm), need to identify which species left, arrived etc. (done above)
# From this, one can identify species affected by CC and those that will not be. Leave the latter alone. 
# To make null estimates of ∆FD/∆PD, one will have to derive null estimates of absolute PD and FD for the future depending on the pattern of SR:
#	- pure decrease in richness (∆SR < 0, no species turn-over)
#   - pure increase in richness (∆SR > 0, no species turn-over)
#	- change in richness with species replacement/ turn-over

### 1st case: ∆SR < 0, no species turn-over
# If you are just loosing species, then your measure of future FD/PD will be equivalent to the FD/PD of the species that are common between T2 and T0. 
# As you cannot randomize the species that remain in the future (= common bulk), you have to randomize the positions of the species being lost for the FD/PD T0.
# You therefore need a null measure of the absolute PD/FD @ T0.
# For that: randomize the position (on the dendrograms) of the species being lost with species that never occurred in the assemblage.

### 2nd case: ∆SR > 0, no species turn-over
# If you are just gaining species, then your measure of FD/PD at T0 will be equivalent to the FD/PD of the species that are common between T2 and T0. 
# As you cannot randomize the species that remain in the present, you have to randomize the positions of the species being gained for the FD/PD at T2.
# You therefore need a null measure of the absolute PD/FD @ T2. 
# For that: randomize the position (on the dendrograms) of the species being gained with species that never occurred in the assemblage.

### 3rd case: ∆SR > 0 or ∆SR < 0, with species turn-over
# Here, you are gaining and loosing certain species at the same time, so you need to randomize both future and present absolute FD/PD estimates.
# Identify the ones that are gone, the ones that have arrived because of CC. 
# Randomize the position of the winners while not allowing common species and lost species to be switched on the dendrograms !
# (i.e. if a species has been lost, you cannot make it possible for it to be switched with a species that has been gained).
# Plus, randomize the tips of the loosers, while not allowing common species and gained species to be switched on the dendrograms !


tot_edges <- function(phy, dropped.tips) {
 	 		# get edges
 		   	e <- phy$edge
  		 	Nedge <- nrow(e)
  		  	# remove some tips
 		 	# keep <- !logical(Nedge)
  	   		# keep[match(dropped.tips, e[,2])] <- FALSE
  	 		keep <- !(e[,2] %in% dropped.tips)
  		  	# detect internal nodes (i.e. not tips)
  		  	Ntip <- length(phy$tip.label)
  		  	is_node <- e[,2] > Ntip
  		  	# remove internal nodes which become tips (i.e. have no children left)
  		  	# repeat {
  		  	#     sel <- !(e[,2] %in% e[,1][keep]) & is_node & keep
  		  	#     if (sum(sel) == 0)
  		  	#         break
  		  	#     keep[sel] <- FALSE
  		  	# }
  		  	sel <- 1
  		  	while(sum(sel) != 0) {
   			 		sel <- !(e[,2] %in% e[,1][keep]) & is_node & keep
    				keep[sel] <- FALSE
  			} # eo while loop
  
  		  	# compute FD
			
			
} # eo tot_edges


### Quick check:
# sum(drop.tip(phy, 1:10)$edge.length)
# tot_edges(phy, 1:10)
# i <- sample(1:nsp, 10)
# sum(drop.tip(phy, i)$edge.length)
# tot_edges(phy, i)
# # -> OK

# the FUN takes 3 arguments: 
# - 'x' : the community table (equivalent to commsum)
# - 'phy' : the functional dendrogram/ the phylogenetic tree ; object of class "phylo"
# - 'n' : number of iterations, for obtaining null distribution of ∆FD/∆PD (ideally, should be 999)
fd_test <- function(x, phy, n= 500) {
  
  			# First, compute observed ∆PD: future PD - present PD
  		  	pd0_obs <- tot_edges(phy, dropped.tips= which(!x %in% c(1,3)) )
  			pd2_obs <- tot_edges(phy, dropped.tips= which(!x %in% c(2,3)) )
  		  	delta_pd_obs <- pd2_obs - pd0_obs
			
  		  	# null hypothesis PD
  		  	set.seed(1)
 		   	# 1) if some species disappear, randomise among the absent (0) and the ones that left (1)
  		  	if (1 %in% x) {
    			sp_to_randomise <- x %in% c(1, 0)
    			x0 <- x
    			pd0 <- numeric(n)
				
    			for (i in 1:n) {
      			  	x0[sp_to_randomise] <- sample(x[sp_to_randomise])
      		  		pd0[i] <- tot_edges(phy, dropped.tips=which(!x0 %in% c(1, 3)))
    			} # eo for loop
  		  	} else {
    			pd0 <- rep(pd0_obs, times= n)
  		  	} # eo if else loop for randomizing T0 tree
  
  		  	# 2) if some species appear, randomise which ones
  		  	if (2 %in% x) {
   			 	sp_to_randomise <- x %in% c(2, 0)
    			x2 <- x
    			pd2 <- numeric(n)
				
    			for (i in 1:n) {
      			  	x2[sp_to_randomise] <- sample(x[sp_to_randomise])
      		  		pd2[i] <- tot_edges(phy, dropped.tips=which(!x2 %in% c(2, 3)))
    			} # eo for loop
  		  	} else {
    			pd2 <- rep(pd2_obs, times= n)
  		  	} # eo if else loop
  		  	
			# Compute null distrib of ∆PD/∆FD
  		  	delta_pd <- (pd2 - pd0)
			# Compute p-value and standardized effect size
 		   	p.value <- sum(delta_pd >= delta_pd_obs) / n
  		  	ses <- (delta_pd_obs - mean(delta_pd)) / sd(delta_pd)

  		  	return( c(p.value, ses) )

} # eo pd_test


### Paralleling according to simulated communities' id

# For testing:
p <- "50"
sdm <- "MARS"
scenar <- "A2ARF"
R <- "RUN3"

fd_null_modeler <- function(sdm = SDMs) {
							
					# For each scenarios
					for(scenar in scenarios) {

						# For each chosen prevalence
						for(p in preval) {

								# For each cv run
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
									
									### Make clever indices to quickly identify the spp. that left/arrived/stayed within the cell after CC
									# - spp that have 1 --> spp present @ T0 only, disappeared becoze of CC
									# - spp that have 2 --> spp present @ T2 only, arrived becoze of CC
									# - spp that have 3 --> spp present @ T0 & T2, not affected by CC
									# - spp that have 0 --> spp that are never in the darn cell
									
									commT2 <- commT2 * 2
									commsum <- commT0 + commT2
									# head(commsum)
									
									### Apply pd_test() described above
									system.time( 
									resNullModel <- apply(commsum, 1, function(x, phy= dendro, n = 500) {
												  
											#cat(".")
											# Fo' testin' 
											#  x <- commsum[1000,]
											#  phy <- dendro
						  					# First, compute observed ∆PD: future PD - present PD
						  		  			pd0_obs <- tot_edges(phy, dropped.tips= which(!x %in% c(1,3)) )
						  					pd2_obs <- tot_edges(phy, dropped.tips= which(!x %in% c(2,3)) )
						  		  			delta_pd_obs <- pd2_obs - pd0_obs
			
						  		  			# null hypothesis PD
						  		  			# set.seed(1)
						  		  			if (1 %in% x) {
						    						sp_to_randomise <- x %in% c(1, 0)
						    					  	x0 <- x
						    						pd0 <- numeric(n)
						    						for (i in 1:n) {
						      			  					x0[sp_to_randomise] <- sample(x[sp_to_randomise])
						      		  						pd0[i] <- tot_edges(phy, dropped.tips = which(!x0 %in% c(1,3)))
						    						} # eo for loop
													
						  		  			} else {
						    					  	pd0 <- rep(pd0_obs, times = n)
						  		  			} # eo if else loop for randomizing T0 tree
  
						  		  			# 2) if some species appear, randomise which ones
						  		  			if (2 %in% x) {
						   			 				sp_to_randomise <- x %in% c(2,0)
						    					  	x2 <- x
						    						pd2 <- numeric(n)
						    					  	for (i in 1:n) {
						      			  					x2[sp_to_randomise] <- sample(x[sp_to_randomise])
						      		  						pd2[i] <- tot_edges(phy, dropped.tips=which(!x2 %in% c(2,3)))
						    						} # eo for loop
													
						  		  			} else {
						    					  	pd2 <- rep(pd2_obs, times= n)
						  		  			} # eo if else loop
  		  	
											# Compute null distrib of ∆PD/∆FD
						  		  		 	delta_pd <- (pd2 - pd0)
											# Compute p-value and standardized effect size from the null distribution od estimates ('delta_pd')
						 		   			p.value <- sum(delta_pd >= delta_pd_obs) / (n+1)
						  		  			ses <- (delta_pd_obs - mean(delta_pd)) / sd(delta_pd)

						  		  			return( c(p.value, ses, delta_pd_obs, pd0_obs, pd2_obs) )
									
									}) # eo apply
									)
									
									gc()
																			
									# go save tables in proper directory
									setwd(paste(WD,"/","Null_Delta_FD_T2","/", sep=  ""))
									save(resNullModel, file = paste("delta_FD_", p, "_", sdm, "_", scenar, "_" , R, ".Rdata", sep= "") )
										
									# Make room
									rm(commT0, commT2, resNullModel)
										
							} # eo cv run
									
					} # eo p in prevalences

			} # eo scenar in scenarios

} # eo FUN										





