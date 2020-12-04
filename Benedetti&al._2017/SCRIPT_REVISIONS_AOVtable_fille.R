
##### 10/05/2016 - ETHZ - Fabio Benedetti 
##### Script for : 
#		- extracting the binomials outputs for T0, T1 & T2 in the MED for the 106 spp.
#		- create communities
#		- compute species richnesses (SRT0 and SRT2), ∆SR (T2-T1) and Jaccard index + components (turn-over + nestedness)
#		- fill multidimensional arrays and save 'em
#		- IMPORTANT : in order to gain time and parallelize the approach, to it for each pabs run separately and parallelize on 10 cores 
#					  (one for each pabs run)
#		- save the 10 marrays as .Rdata ; then read each one, melt and rbind all !!


### Last update : 11/05/2016

# ------------------------------------------------------------------------------------------------------------------------------------

library("stringr")
library("raster")
library("reshape2")
library("sp")
library("matrixStats")
library("ggplot2")
library("RColorBrewer")
library("scales")
library("fields")
library("parallel")
library("biomod2")
library("betapart")
library("plyr")
library("doParallel")
#library("FactoMineR")
#library("vegan")

# ------------------------------------------------------------------------------------------------------------------------------------

### Basic stuff to start with
WD <- getwd()
# For Mediterranean coordinates :
varT0 <- read.table("climatos_woa13_6594_Levin.txt", header=T, sep="\t", dec=".")
# Vector containing cells'id : 
ids <- paste(varT0[,"x"], y = varT0[,"y"], sep="_")

setwd("./niche_modelling")
second.wd <- getwd()
setwd(second.wd)
sp.names <- dir() ; length(sp.names)
sp.names <- sp.names[-c(6:11,49,65,89)]
setwd(WD) ; length(sp.names)

# List of pseudo-absences runs
runs <- as.character(c(1:10)) ; runs
# Cross-validation run
eval_runs <- c("RUN1","RUN2","RUN3")
# List of SDMs
SDMs <- c("SRE","CTA","RF","MARS","FDA","GLM","GAM","ANN","GBM","MAXENT")
# List of binomial experiments 
experiments <- c("exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10")
# List of scenarios 
scenarios <- c("A2","A2F","A2RF","A2ARF","A1BARF","B1ARF")
# List of prevalence levels:
prevalences <- c("n1","n1.5","n2","n10","n50")
# Period list:
#period <- c("T0","T1","T2")
## List of indices that need to be computed : 
indices <- c("SRT0","SRT2","∆SR","Jac","Jtu","Jne","Jratio")
#
experiments <- c("exp1","exp2","exp3","exp4","exp5")

# ------------------------------------------------------------------------------------------------------------------------------------

### For testing : 
 r <- "7"
 sdm <- "GAM"
 scenar <- "A2"
 preval <- "n10"
 R <- "RUN1"
 sp <- "Oncaea.ornata"
 s <- "A2F"

##### !! Missing one species i Binomials: S.bradyi
##### And missing Oncaea ornata projections at T1 and T2....

### Define a function that will fill the marray :

aov_table_filler <- function(r = runs) {

				# Useless message
				message( paste("Extracting values and computing species richness + Jaccard indices for run ",r, sep="") )
				
				marray <- array(data = NA, 
								dim = c(26490, 1, 5, 10, 3, 6, 7), 
					            dimnames= list(ids,
								r,
								prevalences,
					            SDMs,
								eval_runs,
								scenarios,
								indices) 
						  	  	) # eo empty array

				# For each prevalence level
				for(preval in prevalences) {
				
					# For each ENM
					for(sdm in SDMs) {
							
							# For each forcing condition 
							for(scenar in scenarios) {
		
									# For cross-evaluation run
									for(R in eval_runs) {
			
												##### 1°) First step, create community matrices for T0 and T2, 
												##### randomly picking the number of the binomial experiment
												
												### Extract and concatenate species' presence-absences for T0
												p <- "T0"
												speciesT0 <- lapply(sp.names, function(sp) {	
																# Read each species file
																setwd(paste(second.wd,"/","Binomials_proba_T0","/",sp,"/", sep="") ) 
																var <- get(load(paste(sp,"_",preval,"_run",r,"_",p,"_",scenar,".Rdata", sep=""))) 	
																# Sample one og the 5 binomial experiment randomly
																exp <- sample(size= 1, experiments)					
																spT0 <- var[,1,1,sdm,R,exp]
																return(spT0)
												}) #eo lapply T0	
												ddfT0 <- data.frame(do.call(cbind, speciesT0))
												colnames(ddfT0) <- sp.names
												rm(speciesT0)
											
												### Extract and concatenate species' presence-absences for T2
												p <- "T2"
												speciesT2 <- lapply(sp.names, function(sp) {	
																# Read each file
																setwd(paste(second.wd,"/","Binomials_proba_T2","/",sp,"/", sep="") ) 
																var <- get(load(paste(sp,"_",preval,"_run",r,"_",p,"_",scenar,".Rdata", sep="")) )
																# Sample one og the 5 binomial experiment randomly
																exp <- sample(size= 1, experiments)							
																spT2 <- var[,1,1,sdm,R,exp]
																return(spT2)
									
												}) # eo lapply T2	
												ddfT2 <- data.frame(do.call(cbind, speciesT2))
												colnames(ddfT2) <- sp.names
												rm(speciesT2)
												gc()
												
												# Few NAs to be turned into zeros (for pseudo abs runs n° : 1,3 & 7)
												if ( length(ddfT0[is.na(ddfT0)]) > 0 ) {
													ddfT0[is.na(ddfT0)] <- 0
												}
												if ( length(ddfT2[is.na(ddfT2)]) > 0 ) {
													ddfT2[is.na(ddfT2)] <- 0
												} # eo if loops
												
												# Save communities for later...
												setwd(paste(WD,"/","Communities_proba_12_05_16","/", sep=""))
												save(ddfT0, file = paste("Community_table_T0_run_",preval,"_",r,"_",sdm,"_",scenar,"_",R,".Rdata", sep="") )
												save(ddfT2, file = paste("Community_table_T2_run_",preval,"_",r,"_",sdm,"_",scenar,"_",R,".Rdata", sep="") )
												setwd(WD)
												
												##### 2°) Compute vectors of SRT0/T2 (RowSums) and ∆SR (difference) from these communities
												# Compute species richness and add in new dataframe: 
												SRT0 <- rowSums(ddfT0[,c(1:106)], na.rm = TRUE)
												SRT2 <- rowSums(ddfT2[,c(1:106)], na.rm = TRUE)
												# Compute difference in SR between T2 and T0:
												dSR <- (SRT2 - SRT0)
									
												# Supply species richness values to marray
												# marray[,1,"GLM",33,"A2","SRT0"]
												marray[,r,preval,sdm,R,scenar,"SRT0"] <- SRT0
												marray[,r,preval,sdm,R,scenar,"SRT2"] <- SRT2
												marray[,r,preval,sdm,R,scenar,"∆SR"] <- dSR
												
												rm(SRT0,SRT2,dSR) 
												gc()
												
												##### 3°) Compute Jaccard indices and supply them as well : 
												require("dplyr")
												m <- nrow(ddfT0)
												n <- ncol(ddfT0)
												beta.div <- data.frame()
												d <- cbind(ddfT0, ddfT2) 
												d$bit <- cut(1:m, 50, labels=FALSE)
												# Compute ßjaccard + turn-over and nestedness metrics
												registerDoParallel(cores = 2) # Cannot allow too many cores because we are aparalleling on 10 already...
												beta.div <- ddply(d, ~ bit, function(x) {
																beta.temp(x[,1:n], x[,(n+1):(2*n)], "jaccard")
												}, .parallel = TRUE)										
												
												# Supply values to marray : 
												marray[,r,preval,sdm,R,scenar,"Jac"] <- beta.div$beta.jac  # Jaccard index
												marray[,r,preval,sdm,R,scenar,"Jtu"] <- beta.div$beta.jtu  # Turn-over component
												marray[,r,preval,sdm,R,scenar,"Jne"] <- beta.div$beta.jne  # Nestedness component
												marray[,r,preval,sdm,R,scenar,"Jratio"] <- (beta.div$beta.jne) / (beta.div$beta.jac)  # ß ratio
												
												rm(d, beta.div, ddfT0, ddfT2)
												gc()											
														
								} # eo R in eval_runs (usually takes less than 1 min)
									
						} # eo scenar in scenarios (~ a few minutes)
							
				  } # eo sdm in SDMs (~ 25 minutes)
				
			} # eo preval in prevalences	
			
			### Go to scores directory, and print 'eval_marray' as Rdata :
			setwd(paste(WD,"/","ANOVA_tables_proba","/", sep= ""))
			save(marray, file = paste("ANOVA_table_run",r,"_T2.Rdata", sep=""))
			rm(marray)
			setwd(WD)
			
} # eo FUN


### Apply function within mclapply
#runs <- c(7:10)
runs <- as.character(c(7:10)) ; runs
mclapply(runs, aov_table_filler, mc.cores = 10)
gc()



### Considering system time above, it should ~ 2 hours and a half ?

### Check out results ?
#setwd(paste(WD,"/","ANOVA_tables_proba","/", sep= ""))

#d <- get(load(dir()[5]))
#class(d)
#dim(d)
#dimnames(d)
# ...

#d[,"10","n10","SRE","RUN1","A2ARF","Jac"]
#summary( d[,,,"CTA","RUN1","A2","∆SR"] )
#summary(d[,"10","n1","SRE","RUN1","A2ARF","∆SR"])


#count <- lapply(sp.names, function(sp) {
			#setwd(paste(second.wd,"/","Binomials_proba_",p,"/",sp,"/", sep="") )
			#n <- length(dir())
			#return(c(sp,n))
#}) 

#counts <- do.call(rbind, count)
#counts



##### REST IS DONE ON SCRIPT22.4: RESHAPING ANOVA_TABLES WITH dplyr
