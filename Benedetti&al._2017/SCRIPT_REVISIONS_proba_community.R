
#
##### 10/05/2015 - ETHZ - Fabio Benedetti
##### Script for : 
#					- Retrieving the soecies HSI and making binomial experiments on them 
#					--> Will allow to compare thresholded communities to probabilistic communities

 
### Last update : 13/05/2016
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
# List of binomial experiments
experiments <- c("exp1","exp2","exp3","exp4","exp5")
# Period list:
period <- c("T0")

# Species names list:
setwd(paste(second.wd,"/", sep=""))
sp.names <- dir() ; length(sp.names)
sp.names <- sp.names[-c(6:11,49,65,89)]
# sp.names

setwd(second.wd)

################################################################# ¡ For testing !

# sp <- 'Centropages.typicus'
# r <- 1
# s <- "A2F"
# preval <- "n10"
# sdm <- "GAM"
#p <- "T0"
# exp <- 5
# R <- "RUN1"

#issue ith : Oncaea.ornata_n10_run7_T2_A2

################################################################################

setwd(second.wd)


### Defining the function : 
binomial_experiments <- function(sp = sp.names) {
			
			# For each of the 3 time periods:
			for(p in period) {
					message(paste("Running binomials experiments for ",sp," at ",p, sep= ""))
					# For each run of pseudo-absences
					for(r in runs) {
						
						# For each run of pseudo-absences
						for(preval in prevalences) {
				
							# For each scenario
					  		for(s in scenarios) {
								
			  						# Create a multidimensional array for HSI values
									marray <- array(data = NA, 
													dim = c(26490, 1, 1, 10, 3, 5), 
        											dimnames= list(
													c(1:26490),
													c("x"), 
													c("y"),
          											SDMs,
													eval_runs,
													experiments))

									### Supply coordinates to the multidimensional array first & once
									marray[,1,,,,] <- varT0$x
									marray[,,1,,,] <- varT0$y
					
								
									# Fill it within the SDMs loop:				
									for(sdm in SDMs) {
									
											# Go to species directory depending on considered period :
											setwd( paste(second.wd,"/",sp,"/","proj_",sp,"_",preval,"_run_",r,"_",p,"/", sep= "") )
											# Load proper model ouputs
											var <- load(paste("proj_",sp,"_",preval,"_run_",r,"_",p,"_",sp,".RData", sep=""))    															
											# Get results
											resModel <- get(var)
											rm(list=(var))
													
											for(R in eval_runs) {			
													# On extrait les résultats des différents runs et on sort la valeur moyenne sans tenir compte des NAs												
													res <- resModel[,sdm,R,]
													#res <- apply(res, 1, mean, na.rm = TRUE) 
													res <- (res/1000) 
													### Supply HSI values to the multidimensional array !
													marray[,1,1,sdm,R,] <- res
													
													# Supply the binomials experiments to the multidimensional array !
													for(i in 1:26490) {
															marray[i,1,1,sdm,R,experiments] <- t(rbinom(n= 5, size= 1, prob = marray[i,1,1,sdm,R,] ))
													} # for i in cells
																										
											} #eo for R in eval_runs				
											gc()
										# Check 
										# marray[,1,1,sdm,R,] 
										# Make room									
												
									} #eo sdm in SDMs (~ Less than 30 seconds)		
								
									rm(resModel, res)				
								
									# Go to proper dir and save files there:  			
									message(paste("Saving binomials experiments for ",sp," at ",p," for prevalence ",preval, sep=""))
										  				  
									setwd(paste(second.wd,"/","Binomials_proba_",p,"/", sep="") )
									# if first run, create the species dir, if not, got to it
									if(r == 1) {
										dir.create(paste(sp)) 
										setwd(paste(second.wd,"/","Binomials_proba_",p,"/",sp,"/", sep="") )
									} else {
										setwd(paste(second.wd,"/","Binomials_proba_",p,"/",sp,"/", sep="") )
									}
								
									# Save multi-array as Rdata:
									save(marray, file = paste(sp,"_",preval,"_run",r,"_",p,"_",s,".Rdata",sep="") )
									rm(marray)
							
						} #eo s in scenarios (~ 1min30sec)	
					
					} #eo preval in prevalence
					
				} #eo r in runs
		
			} #eo p in period 
			
		# Go back to base WD
		setwd(WD)
			
} #eo FUN


### Apply fun within mclapply:
mclapply(sp.names, binomial_experiments, mc.cores = 25)



