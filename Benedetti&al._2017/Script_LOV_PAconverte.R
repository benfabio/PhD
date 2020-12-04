
##### 23/04/2015 - LOV - Fabio Benedetti 
##### Script for : 
#					- creating a function that executes 10 binomial experiments for each projection and saves them in a ddf
#					- execute the function with a mclapply()


### Last update : 08/06/2015

# ------------------------------------------------------------------------------------------------------------------------------------

library("stringr")
library("raster")
library("plyr")
library("reshape2")
library("sp")
library("matrixStats")
library("ggplot2")
library("scales")
library("fields")
library("parallel")

# ------------------------------------------------------------------------------------------------------------------------------------

WD <- getwd()

# Pour coordonnées en Méd :
varT0 <- read.table("climatos_woa13_6594_Levin.txt", header=T, sep="\t", dec=".")

# ------------------------------------------------------------------------------------------------------------------------------------

##### First, retreive a dataframe of HSI and convert it to raster : 

# register second.wd
setwd("./niche_modelling")
second.wd <- getwd()
setwd(second.wd)
#dir.create("Binomials_T0_glob")
sp.names <- dir() ; length(sp.names)
sp.names <- sp.names[-c(6,7,45,53,62,86)]
setwd(WD) ; length(sp.names)

# List of pseudo-absences runs
runs <- as.character(c(1:10)) ; runs
# List of emission scenarios
scenarios <- c("A2","A2F","A2RF","A2ARF","A1BARF","B1ARF")
# List of periods/ times:
# times <- c("T0","T1","T2")
# List of SDModels
SDMs <- c("SRE","CTA","RF","MARS","FDA","GLM","GAM","ANN","GBM","MAXENT")
# List of binomial experiments
experiments <- c("exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10") ; gc()


# ------------------------------------------------------------------------------------------------------------------------------------

### 24/04/2015
### First version of the script takes way too long, especially the loading and printing of the big data tables (read.table & write.table functions)
### Try to correct that by saving multidimensional arrays as Rdata ; for instance, you can create an empty multi-array of the right dimensions
### (i.e. Lon/Lat, SDM, binomial experiement based on HSI) before the sdm for loop. 
### 	- 1 dimension for longitude & latitude 
###		- 1 dimension for SDM choice (10 categories) with HSI values 
###		- 1 dimension for binomial experiments (10 categories too) based on the HSI values

##### 1st step : how to build a multidimensional array ? Yeah sounds dumb but I've never done it before.
?array
marray <- array(data = NA, 
				dim = c(26490, 2, 10, 10), 
              	dimnames= list(
				c(1:26490),
				c("Lon", "Lat"), 
                SDMs,
                experiments))
				
dim(marray)
str(marray)
length(marray)
summary(marray)
head(marray[,2,3,9])

### Test :
setwd(second.wd)
sp <- "Pseudocalanus.elongatus"
sdm <- "GBM"
r <- 1
s <- "A2"


### Defining the function for T1 & T2 : 

Rose <- function(sp = sp.names) {

			message(paste("Running binomials experiments for ",sp, sep=""))
			
			# For each run of pseudo-absences
			for(r in runs) {
				
					  for(s in scenarios) {
								
						  		# Create a multidimensional array for HSI values
								marray <- array(data = NA, 
												dim = c(26490, 1, 1, 10, 10), 
			              						dimnames= list(
												c(1:26490),
												c("x"), 
												c("y"),
			                					SDMs,
												experiments))

								### Supply coordinates to the multidimensional array first & once
								marray[,1,,,] <- varT0$x
								marray[,,1,,] <- varT0$y
								# dim(marray)						
								
								# Fill it within the SDMs loop:				
								
								for(sdm in SDMs) {
								
											# Go to species directory depending on considered period :
											setwd( paste(second.wd,"/",sp,"/","proj_",sp,"_run_",r,"_","T2","_",s,"/", sep="") )

											# Load proper model ouputs
											var <- load(paste("proj_",sp,"_run_",r,"_","T2","_",s,"_",sp,".RData", sep=""))    
															
											# Get results
											resModel <- get(var)
											rm(list=(var))
												
											# On extrait les résultats des différents runs et on sort la valeur moyenne sans tenir compte des NAs												
											res <- resModel[,sdm,,]
											res <- apply(res, 1, mean, na.rm = TRUE) 
											res <- (res/1000) 
											ddf <- data.frame(x = varT0$x, y = varT0$y, p = res)
											# summary(ddf)
										
											### Supply HSI values to the multidimensional array !
											marray[,1,1,sdm,] <- ddf$p
											
											# Supply the binomials experiments to the multidimensional array !
											for(i in 1:26490) {
												marray[i,1,1,sdm,experiments] <- t(rbinom(n= 10, size= 1, prob = marray[i,1,1,sdm,] ))
											}
											gc()

											# Make room
											rm(ddf, resModel, res)											
												
								} #eo sdm in SDMs 		
								
								# Go to proper dir and save files there:  					  				  
								setwd(paste(second.wd,"/","Binomials_T2","/", sep="") )
								# if first run, create the species dir, if not, got to it
								if(r == 1) {
									dir.create(paste(sp)) 
									setwd(paste(second.wd,"/","Binomials_T2","/",sp,"/", sep="") )
								} else {
									setwd(paste(second.wd,"/","Binomials_T2","/",sp,"/", sep="") )
								}
								
								# Save multi-array as Rdata:
								save(marray, file = paste(sp,"_run",r,"_T2_",s,".Rdata",sep="") )
								rm(marray)
							
					} #eo s in scenarios	
					
		} #eo r in runs
			
		# Go back to base WD
		setwd(WD)
			
} #eo FUN



### Apply fun within mclapply:
# sp.names <- sp.names[-c(35)]
sp.names <- c("Corycaeus.furcifer","Subeucalanus.crassus")		
# sp.names <- c("Triconia.minuta")
mclapply(sp.names, Rose, mc.cores = 2)		

ls()
dim(marray)
marray[,,,"GLM","exp1"]

ddf1 <- data.frame(x = varT0$x, y = varT0$y, p = marray[,1,1,"GLM","exp1"]) 
ddf2 <- data.frame(x = varT0$x, y = varT0$y, p = marray[,1,1,"GAM","exp1"]) 

### BINUS : Create a custom color scale
colScale <- scale_colour_manual(name = "Modelled\nPresence/Absence", values = c("#63B8FF","#FF4040")) # Blue & Red
colScale <- scale_colour_manual(name = "Modelled\nPresence/Absence", values = c("#E0EEE0","#FF4040")) # Whitish & Red
colScale <- scale_colour_manual(name = "Modelled\nPresence/Absence", values = c("#87CEFF","#FF4040")) #  Skyblue & Red

plot1 <- ggplot(ddf1) + geom_point(aes(x=x, y=y, colour= factor(p)), alpha = 0.8) + colScale + coast + coord_quickmap() ; plot1
plot2 <- ggplot(ddf2) + geom_point(aes(x=x, y=y, colour= factor(p))) + coast + colScale + coord_quickmap() ; plot2

ggsave(plot = plot1, filename= paste("map","P/A","T1",sp,r,"GLM.png", sep="_"),  width=13, height=7, dpi= 300)

# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------


### 28/04/2015
### Defining the function for T0 Mediterranean scale : 

### Test :
setwd(second.wd)
#sp <- "Temora.stylifera"
#sdm <- "GBM"
#r <- 1


Sun <- function(sp = sp.names) {

			message(paste("Running binomials experiments for ",sp, sep=""))
			
			# For each run of pseudo-absences
			for(r in runs) {
								
					# Create a multidimensional array for HSI values
					marray <- array(data = NA, 
									dim = c(26490, 1, 1, 10, 10), 
			              			dimnames= list(
									c(1:26490),
									c("x"), 
									c("y"),
			                		SDMs,
									experiments))

					### Supply coordinates to the multidimensional array first & once
					marray[,1,,,] <- varT0$x
					marray[,,1,,] <- varT0$y
					# dim(marray)						
								
					# Fill it within the SDMs loop:				
								
					for(sdm in SDMs) {
								
								# Go to species directory depending on considered period :
								setwd( paste(second.wd,"/",sp,"/","proj_",sp,"_run_",r,"_","T0","/", sep="") )

								# Load proper model ouputs
								var <- load(paste("proj_",sp,"_run_",r,"_","T0","_",sp,".RData", sep=""))    
															
								# Get results
								resModel <- get(var)
								rm(list=(var))
												
								# On extrait les résultats des différents runs et on sort la valeur moyenne sans tenir compte des NAs												
								res <- resModel[,sdm,,]
								res <- apply(res, 1, mean, na.rm = TRUE) 
								res <- (res/1000) 
								ddf <- data.frame(x = varT0$x, y = varT0$y, p = res)
								# summary(ddf)
										
								### Supply HSI values to the multidimensional array !
								marray[,1,1,sdm,] <- ddf$p
											
								# Supply the binomials experiments to the multidimensional array !
								for(i in 1:26490) {
										marray[i,1,1,sdm,experiments] <- t(rbinom(n= 10, size= 1, prob = marray[i,1,1,sdm,] ))
								}
								gc()

								# Make room
								rm(ddf, resModel, res)											
												
					} #eo sdm in SDMs 		
								
					# Go to proper dir and save files there:  					  		  
					setwd(paste(second.wd,"/","Binomials_T0","/", sep="") )
					# if first run, create the species dir, if not, got to it
					if(r == 1) {
						dir.create(paste(sp)) 
						setwd(paste(second.wd,"/","Binomials_T0","/",sp,"/", sep="") )
					} else {
						setwd(paste(second.wd,"/","Binomials_T0","/",sp,"/", sep="") )
					}
								
					# Save multi-array as Rdata:
					save(marray, file = paste(sp,"_run",r,"_T0_",".Rdata",sep="") )
					rm(marray)
					
			} #eo r in runs
			
		# Go back to base WD
		setwd(WD)
			
} #eo FUN


### Apply fun within mclapply:
# sp.names <- sp.names[-c(35,100)]
sp.names <- c("Corycaeus.furcifer","Subeucalanus.crassus")	
mclapply(sp.names, Sun, mc.cores = 2)	

### mapping
ddf1 <- data.frame(x = varT0$x, y = varT0$y, p = marray[,1,1,"SRE","exp2"])
ddf2 <- data.frame(x = varT0$x, y = varT0$y, p = marray[,1,1,"GLM","exp1"])

ggplot(ddf1) + geom_point(aes(x=x, y=y, colour= factor(p))) + coast + coord_quickmap()
ggplot(ddf2) + geom_point(aes(x=x, y=y, colour= factor(p))) + coast + coord_quickmap()

### Checking C.furcifer :
marray[i,1,1,sdm,experiments]
ddf1 <- data.frame(x = varT0$x, y = varT0$y, p = marray[,1,1,"GBM","exp1"])
ggplot(ddf1) + geom_point(aes(x=x, y=y, colour= factor(p))) + coast + coord_quickmap()



# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------

### Defining the function for T0 global scale : 

### Test :
setwd(second.wd)
sp <- "Parapontella.brevicornis"
# sdm <- "GBM"
# runs <- c(6)
r <- 2

# 
Peak <- function(sp = sp.names) {

			message(paste("Running global binomials experiments for ",sp, sep=""))
			
			# For each run of pseudo-absences
			for(r in runs) {	
									
					# Load proper model ouputs				
					setwd( paste(second.wd,"/",sp,"/","proj_",sp,"_run_",r,"_","T0_glob","/", sep="") )
					ras <- brick(paste("proj_",sp,"_run_",r,"_T0_glob_",sp,".grd", sep= ""))
					ddf <- as.data.frame(ras, xy= TRUE)
					# dim(ddf)
								
					# Give simpler column names :
					colnames <- names(ras)
					colnames <- str_replace_all(colnames, paste(sp,"_AllData_", sep=""), "") # getting rid of unused text
					L <- length(ddf)	
					colnames(ddf)[3:L] <- colnames
								
					# summary(ddf)
					# str(ddf)
					# dim(ddf)
					ddf <- na.omit(ddf)
					n <- nrow(ddf)	
					n			
					
					# Create a multidimensional array for HSI values with right dimensions: 
					marray <- array(data = NA, 
									dim = c(n, 1, 1, 10, 10), 
			              			dimnames= list(
									c(1:n),
									c("x"), 
									c("y"),
			                		SDMs,
									experiments))		
									
					### Supply coordinates to the multidimensional array first & once
					coords <- data.frame(x= ddf$x, y= ddf$y)
					marray[,1,,,] <- ddf$x
					marray[,,1,,] <- ddf$y						
								
					# Fill it within the SDMs loop:					
					for(sdm in SDMs) {
						
								# Convert to numerical matrix for quicker computing
								mat <- as.matrix(ddf[,c(1:L)]) 			
								gc()

								# On extrait les résultats des différents runs et on sort la valeur moyenne sans tenir compte des NAs	
								##### !!! Sometimes, mean cannot be computed because we only have 1 SDM RUN available (GBM usually)
								##### Need to add a if loop : do not compute mean if only 1 column, rather keep the values from the
								##### only one available.
								
								# If only one column, then keep value, otherwise compute rowmeans
								l <- length(grep(sdm,colnames(mat)))
								
								if(l == 1) {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = mat[,grep(sdm,colnames(mat))])
								} else {
																			
										if(sdm == "SRE") { 
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm,colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "CTA") { 
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "RF") {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "GBM") {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "GLM") {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "GAM") {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "MARS") {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "FDA") {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "ANN") {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} else if (sdm == "MAXENT") {
											ddf2 <- data.frame(x= mat[,"x"], y= mat[,"y"], p = rowMeans(mat[,grep(sdm, colnames(mat))], na.rm = TRUE) )
										} # eo if else # 2
										
								} # eo if else # 1
								
								gc() # purge !
								
								# summary(ddf2)
								
								### Divide the probabilities by 1000 to scale between 0 and 1:
								ddf2$p <- (ddf2$p)/1000
										
								### Supply HSI values to the multidimensional array !
								marray[,1,1,sdm,] <- ddf2$p
											
								# Supply the binomials experiments to the multidimensional array !
								for(i in 1:n) {
									  marray[i,1,1,sdm,experiments] <- t(rbinom(n= 10, size= 1, prob = marray[i,1,1,sdm,] ))
								}
								gc()								
												
					} #eo sdm in SDMs 		
								
					# Go to proper dir and save files there:  					  
					setwd(paste(second.wd,"/","Binomials_T0_glob","/", sep="") )
					# if first run, create the species dir, if not, got to it
					if(r == 1) {
						dir.create(paste(sp)) 
						setwd(paste(second.wd,"/","Binomials_T0_glob","/",sp,"/", sep="") )
					} else {
						setwd(paste(second.wd,"/","Binomials_T0_glob","/",sp,"/", sep="") )
					} # eo if else
								
					# Save multi-array as Rdata in species dir :
					save(marray, file = paste(sp,"_run",r,"_T0_glob",".Rdata",sep="") )
					rm(ddf2, L, ras, mat, l, n, marray, ddf) # make room
					
			} #eo r in runs
			
		# Go back to base WD
		setwd(WD)
			
} #eo FUN


### Apply fun within mclapply:
sp.names <- "Mormonilla.phasma"
runs <- c(7) ; runs
mclapply(sp.names, Peak, mc.cores = 1)


### Checking :
setwd(paste(second.wd,"/","Binomials_T0","/",sp,"/", sep="") ) ; dir()
# Load proper model ouputs
var <- load(paste(sp,"_run",r,"_T0_",".Rdata",sep=""))    						
# Get results
marray <- get(var)
# Map
ddf3 <- data.frame(x = coords$x, y = coords$y, p = marray[,1,1,"GBM","exp1"])
ggplot(ddf3) + geom_point(aes(x=x, y=y, colour= factor(p))) + coord_quickmap() + theme_bw() + world_coast



