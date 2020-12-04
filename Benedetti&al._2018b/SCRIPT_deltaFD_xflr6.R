

### Fabio Benedetti - François Guilhaumon - Jean-Olivier Irisson
### UPMC/LOV             IRD/MARBEC/UM2

### R SCRIPT BATCH to compute ∆FD and ∆PD indices on xflr6

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
runs <- as.character(c(1:3))
# List of SDMs
SDMs <- c("SRE","CTA","RF","MARS","FDA","GLM","GAM","ANN","GBM","MAXENT")
# List of scenarios 
scenarios <- c("A2","A2F","A2RF","A2ARF","A1BARF","B1ARF")
# List of potential communities :
communities <- sample(c(51:100), size = 5)

### Community tables: species' P/A table in the Med Sea. Just lacking coordinates...let's get 'em back
setwd(paste(WD,"/","ANOVA_tables","/","All","/", sep=""))
numbers <- c(1:10)
res <- lapply(numbers, function(n) {
			d <- get(load(paste("means_SRT0_subtable_T2_",n,".Rdata", sep= "")))
			ids <- data.frame(unique(d$id))
			colnames(ids) <- "id"
			return(ids)
			}
) # eo lapply

ids <- do.call(rbind, res) ; rm(res)
#head(ids)
#unique(ids)
### Retrieve coordinates
coords <- data.frame(str_split_fixed(ids$id, "_", 2))
colnames(coords)[1:2] <- c("Lon","Lat")
coords$Lon <- as.numeric(levels(coords$Lon))[coords$Lon]
coords$Lat <- as.numeric(levels(coords$Lat))[coords$Lat]
rm(ids)


# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------


### First, load the phylogenetic tree and functional dendrogram containing the 106 spp.
#setwd(second.wd)
#sp.names <- dir() ; length(sp.names)
#sp.names <- sp.names[-c(6:10,48,56,65,89,116)]
#setwd(WD) ; length(sp.names)

### A) Phylogenetic tree
setwd(WD)
library("ape")
phylo <- ape::read.nexus("Tree_cons_189.nex")
phylo$tip.label
# Make sure species names are the same for all datasets
phylo$tip.label <- gsub("_", ".", phylo$tip.label , fixed= T)
phylo$tip.label[1:6] <- c("Acartia.Acanthacartia.bifilosa", "Acartia.Acartiura.clausi", 
						  "Acartia.Acartia.danae", "Acartia.bifilosa", "Acartia.Acartia.negligens", "Acartia.Acanthacartia.tonsa")
						  
phylo$tip.label[65:74] <- c("Corycaeus.brehmi","Corycaeus.clausi","Corycaeus.Agetus.flaccus","Corycaeus.furcifer","Corycaeus.giesbrechti",
							"Corycaeus.latus","Corycaeus.Agetus.limbatus","Corycaeus.ovalis","Corycaeus.speciosus","Corycaeus.Agetus.typicus")
						  
# Drop the tips that do not belong to your Gamma diversity (106 spp only)
dropped <- phylo$tip.label[!(phylo$tip.label %in% sp.names)]
phylo2 <- drop.tip(phy= phylo, tip= dropped )
phylo2$tip.label
### Only 88 tips on the phylogeny, not 106.
# Who is lacking, let's see:
sp.names[!sp.names %in% phylo2$tip.label]
# yep...makes sense, noone of these genera have been sequenced.


### B) Functional dendrogram
fct <- read.csv("table_traits_biogeography.csv", h=T, dec=",", sep=";")
fct <- fct[,c(1,65,66,67,68)]
# Make sure species names are the same for all datasets
fct$sp.name <- sp.names
### Need to make a functional dendrogram for future randomization of the specie stips 
distance <- dist((fct[,c("length_max")]), method = "euclidean") # distance matrix
fit <- hclust(distance, method= "ward") 
groups <- cutree(fit, k=4) # cut tree into 4 clusters
fct$size_class <- factor(groups)

# Convert trophism into new classes
fct$Carnivore <- (fct$trophism=="Carnivore")|(fct$trophism=="Omnivore-Carnivore")
fct$Omnivore <- (fct$trophism=="Omnivore")|(fct$trophism=="Omnivore-Carnivore")|(fct$trophism=="Omnivore-Detritivore")|(fct$trophism=="Omnivore-Herbivore")
fct$Detritivore <- (fct$trophism=="Omnivore-Detritivore")
fct$Herbivore <- (fct$trophism=="Omnivore-Herbivore")
rownames(fct) <- fct$sp.name
library("FactoMineR")
resMCA <- MCA(fct[,c(4:10)], ncp=4, na.method="Average")
resMCA$ind
dist_mca <- dist(data.frame(resMCA$ind$coord), method = "euclidean")
my.traits.mca <- resMCA$ind$coord

# Clustering to create dendrogram
library("fastcluster")
clust_mca <- fastcluster::hclust(dist_mca,method="average")
clust_mca <- as.phylo.hclust(clust_mca)
#clust_mca$tip.label

# Check common tip labels between fct table and phylo2$tip.label
fct$sp.name[!fct$sp.name %in% phylo2$tip.label] # 18, same as phylo2 versus sp.names object == colnames in assemblage tables
# gut

#save(clust_mca, file = "fct_dendrogram.Rdata")

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


### We will re-use the scripts functions implemented by J-O Irisson (test_tree_computation.R) to apply them for our data. 
### The point is to shortcut the drop.tip() FUN from 'ape' which takes the longest --> tot_edges()
# 1) Compute total phylogenetic distance on a partial tree
# 2) Some tips are dropped
# 3) When all childs of an ancestor are dropped, the edge leading to the ancestor must also be removed

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
  
  		  	# compute total phylogenetic distance with sum()
  		  	tot_pd <- sum(phy$edge.length[keep])
  		  	return(tot_pd)  
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
pd_test <- function(x, phy, n= 500) {
  
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


### takes ~ 23min per commsum. 
### 23min x 10 ENMs x 3 SRES x 3 pabs x 10 comms x 2 indices x time periods / 15 cores --> ~ 1840 min --> ~ 31hrs --> ~ 1,3 day
# If you consider 3 pabs runs + 20 core parallelisation --> 2.9 days
# If you consider only end of century (T2 vs T0) + 20 cores --> 1.5 days
# If you consider only end of century (T2 vs T0) + 15 cores --> 1.9 days

### Choosing nrandom == 500 goes faster than == 999 (of course)
# if n == 500 --> ~ 23min per table (see stats above)
# if n == 999 --> ~ 60min per table ! 
# When n == 999 you might want to reduce number of ENM algorithms (7-8 instead of 10) and parallelise on 20 cores.

### Paralleling according to simulated communities' id
null_modeler <- function(sdm = SDMs) {
							
						# For each SDMs
						for(comm in communities) {
							
							# For RCM bounday forcing condition 
							for(scenar in scenarios) {
		
								# For each pseudo-absence iteration 
								for(r in runs) {
										
									# Useless message
									message( paste("Doing_run",r,"_",sdm,"_",scenar,"_",comm,".Rdata", sep= "") )
									
									# Loading community tables
									setwd(paste(WD,"/","Communities_28_07_15","/", sep=""))
									commT0 <- get(load(paste("Community_table_T0_run",r,"_",sdm,"_",scenar,"_",comm,".Rdata", sep="")))
									commT2 <- get(load(paste("Community_table_T2_run",r,"_",sdm,"_",scenar,"_",comm,".Rdata", sep="")))
									
									# Need to re-supply cells' ids
									rownames(commT0) <- paste(coords$Lon, coords$Lat, sep= "_")
									rownames(commT2) <- paste(coords$Lon, coords$Lat, sep= "_")
									ids <- rownames(commT2)
									
									### Make clever indices to quickly identify the spp. that left/arrived/stayed within the cell after CC
									commT2 <- commT2 * 2
									commsum <- commT0 + commT2
									# head(commsum)
									
									### Apply pd_test() described above
									#system.time( 
									resNullModel <- apply(commsum, 1, function(x, phy= clust_mca, n= 500) {
												  
											#cat(".")
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
						      		  						pd2[i] <- tot_edges(phy, dropped.tips=which(!x2 %in% c(2,3)))
						    						} # eo for loop
													
						  		  			} else {
						    					  	pd2 <- rep(pd2_obs, times= n)
						  		  			} # eo if else loop
  		  	
											# Compute null distrib of ∆PD/∆FD
						  		  		 	delta_pd <- (pd2 - pd0)
											# Compute p-value and standardized effect size
						 		   			p.value <- sum(delta_pd >= delta_pd_obs) / (n+1)
						  		  			ses <- (delta_pd_obs - mean(delta_pd)) / sd(delta_pd)

						  		  			return( c(p.value, ses, delta_pd_obs) )
									
									}) # eo apply
									#)
									
									gc()
																			
									# go save tables in proper directory
									setwd(paste(WD,"/","Null_Delta_FD_T2","/", sep=""))
									save(resNullModel, file = paste("mat_deltaFD",comm,sdm,scenar,r,".Rdata", sep="_") )
										
									# Make room
									rm(commT0, commT2, resNullModel)
										
								} # eo r in runs
									
							} # eo scenar in scenarios

						} # eo sdm in SDMs

} # eo FUN										

mclapply(X= SDMs, FUN= null_modeler, mc.cores = 15)
gc()
