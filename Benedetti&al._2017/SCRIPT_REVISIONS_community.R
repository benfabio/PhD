
#
##### 10/05/2015 - ETHZ - Fabio Benedetti
##### Script for : 
#					- Retrieving the TSS bins (P/A maps) for each SDM projections
#					- Save community matrices for each combination of: ENM, forcing, SRES, prevalence etc.
#					- Compure ∆SR and betadiv metrics ! 

 
### Last update : 10/05/2016
# ---------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("sp")
library("ggplot2")
library("stringr")
library("reshape2")
library("dplyr")
library("PresenceAbsence")
library("SDMTools")
library("biomod2")
library("parallel")
library("ggplot2")
library("fields")

WD <- getwd()

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Loading environmental predictors, present and future Mediterranean fields, and species' presences/ pseudoAbsences files
varT0_glob <- brick("clims_woa13_baseline_04.grd")
varT0_glob <- dropLayer(varT0_glob, c(1,2,6))
# varT0_glob

# Pour projection/prédiction : climatos WOA13 Mediterranée et NM8
varT0 <- read.table("climatos_woa13_6594_Levin.txt", header=T, sep="\t", dec=".")

# On crée le dossier dans lequel on va stocker tous les résultats de la modélo et on y dépose gentiment le fichier java de maxente
initial.wd <- getwd()
setwd("./niche_modelling")
second.wd <- getwd()

# Liste des scenarios
scenarios <- c("A2","A2F","A2RF","A2ARF","A1BARF","B1ARF")
# List des runs de pseudoAbs
runs <- c(1:10)
# ENMs list:
SDMs <- c('SRE','CTA','RF','GBM','GLM','GAM','MARS','FDA','ANN','MAXENT')
# List of eval_runs:
eval_runs <- c("RUN1","RUN2","RUN3") ## !! need an underscore because biomod2 adds one by default in pred and calib_lines !! 
# List of prevalence levels:
prevalences <- c("n1","n1.5","n2","n10","n50")
# Period list:
period <- c("T0","T1","T2")

# Species names list:
setwd(paste(second.wd,"/", sep=""))
sp.names <- dir() ; length(sp.names)
sp.names <- sp.names[-c(6,7,8,46,62,86)]
#sp.names

setwd(second.wd)

################################################################# ¡ For testing !

# sp <- 'Calanus.helgolandicus'
# r <- 1
# s <- "A2"
# preval <- "n1.5"
# sdm <- "GAM"
# p <- "T2"

################################################################################


### Defining the function : 
getTSSbins <- function(sp = sp.names) {

			for(p in period) {
			message(paste("Extracting TSS bins for ",sp," @ ",p, sep=""))
			
				# For each prevalence level
				for(r in runs) {
			
					# For each run of pseudo-absences
					for(preval in prevalences) {
					  
					  		# For each scenario
					  		for(s in scenarios) {
								
						  			# Create a multidimensional array for HSI values
									marray <- array(data = NA, 
													dim = c(26490, 1, 1, 10, 3), 
	              									dimnames= list(
													c(1:26490),
													c("x"), 
													c("y"),
	                								SDMs,
													eval_runs))

									### Supply coordinates to the multidimensional array first & once
									marray[,1,,,] <- varT0$x
									marray[,,1,,] <- varT0$y
								
									# Fill it within the SDMs loop:				
									for(sdm in SDMs) {
								
											# Go to species directory depending on considered period :
											# Beware: T0 projectiosn do not have a scenario index ! 
											if( p == "T0" ) {						
													setwd( paste(second.wd,"/",sp,"/","proj_",sp,"_",preval,"_run_",r,"_",p,"/", sep="") )
													# Load proper model ouputs
													var <- get(load(paste("proj_",sp,"_",preval,"_run_",r,"_",p,"_",sp,"_TSSbin.RData", sep="")) ) 														
	
													# On extrait les résultats des différents runs et on sort la valeur moyenne sans tenir compte des NAs												
													res <- var[,sdm,,]
													res <- apply(res, 2, as.integer, na.rm = F) 
													
													### Watchout for C. helgolandicus and its messed up runs...
													if(sp == "Calanus.helgolandicus" && sdm == "GAM" && r == 10 && preval == "n1.5") {
															# Supply P/A to marray accounting for differences in RUN name :
															marray[,,,sdm,"RUN1"] <- res[,"RUN1"]
															marray[,,,sdm,"RUN2"] <- res[,"RUN1"]
															marray[,,,sdm,"RUN3"] <- res[,"RUN1"]		
													} else {
															# Supply P/A to marray normally :
															marray[,,,sdm,"RUN1"] <- res[,"RUN1"]
															marray[,,,sdm,"RUN2"] <- res[,"RUN2"]
															marray[,,,sdm,"RUN3"] <- res[,"RUN3"]
													}
													# Make room, biatch
													rm(res)		
													
											} else {
													setwd( paste(second.wd,"/",sp,"/","proj_",sp,"_",preval,"_run_",r,"_",p,"_",s,"/", sep="") )
													# Load proper model ouputs
													var <- get(load(paste("proj_",sp,"_",preval,"_run_",r,"_",p,"_",s,"_",sp,"_TSSbin.RData", sep="")) ) 														

													# On extrait les résultats des différents runs et on sort la valeur moyenne sans tenir compte des NAs												
													res <- var[,sdm,,]
													res <- apply(res, 2, as.integer, na.rm = F) 
												
													### Watchout for C. helgolandicus and its messed up runs...
													if(sp == "Calanus.helgolandicus" && sdm == "GAM" && r == 10 && preval == "n1.5") {
															# Supply P/A to marray accounting for differences in RUN name :
															marray[,,,sdm,"RUN1"] <- res[,"RUN1"]
															marray[,,,sdm,"RUN2"] <- res[,"RUN1"]
															marray[,,,sdm,"RUN3"] <- res[,"RUN1"]
													} else {
															# Supply P/A to marray normally :
															marray[,,,sdm,"RUN1"] <- res[,"RUN1"]
															marray[,,,sdm,"RUN2"] <- res[,"RUN2"]
															marray[,,,sdm,"RUN3"] <- res[,"RUN3"]
													}
													# Make room, biatch
													rm(res)		
											
											} # eo period if else loop									
												
									} # eo sdm in SDMs 		
									
									# marray[,,,"SRE","RUN1"]
								
									# Go to proper dir and save files there:  					  				  
									setwd(paste(second.wd,"/","Binomials","_",p,"/", sep="") )
									# if first run, create the species dir, if not, got to it
									if( r == 1 ) {
											dir.create(paste(sp)) 
											setwd(paste(second.wd,"/","Binomials","_",p,"/",sp,"/", sep="") )
									} else {
											setwd(paste(second.wd,"/","Binomials","_",p,"/",sp,"/", sep="") )
											# length(dir())
									} # eo if else loop
								
									# Save multi-array as Rdata:
									message(paste( paste(sp,"_",preval,"_run",r,"_",p,"_",s, sep="") ))
									save(marray, file = paste(sp,"_",preval,"_run",r,"_",p,"_",s,".Rdata", sep="") )
									rm(marray)
									gc()
							
						} #eo s in scenarios	
					
				} #eo r in runs
		
		} #eo preval in prevalences
		
	} # eo p in period
			
	# Go back to base WD
	setwd(WD)
			
} #eo FUN


### Apply fun within mclapply:
mclapply(sp.names, getTSSbins, mc.cores = 20)



