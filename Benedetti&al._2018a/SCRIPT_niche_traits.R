
##### 17/01/2017 - LOV - Fabio Benedetti
##### Script for : 
#					- conducting the sensitivity analyses prior to manucsript submission to L&O
#					- assessing niche traits variations due to smoothing parameter
#					- assessing niche traits distortion due to data omission (80%, 70% and 50% of the occurrence data)

### Last update : 22/02/2017

# -------------------------------------------------------------------------------------------------------------

library("PresenceAbsence")
library("biomod2")
#library("raster")
#library("rgdal")
library("sp")
#library("plyr")
library("dplyr")
library("stringr")
library("reshape2")
library("cclust")
library("fields")
library("scales")
library("gbm")
library("ggplot2")
library("RColorBrewer")
library("boot")

# ---------------------------------------------------------------------------------------------------------------------------------------------

WD <- getwd()
setwd(WD)

# List of chosen predictors:
vars <- c("SST","dSST","SSS","MLD1","MLPAR1","logChla")
# List of species names: 
setwd(paste(WD,"/","niche.modelling","/","pabs_files_matched","/", sep=""))
sp.names <- dir()
sp.names <- str_replace_all(sp.names, "pabs_", "")
sp.names <- str_replace_all(sp.names, "_all.txt", "")
sp.names <- sp.names[13:106]

# List of smoothing parameters, k values
#klevels <- c(3,10,15)

# List of data omission degrees
omis <- c("80p","70p","50p")

setwd(WD)
setwd("./niche.modelling")
second.wd <- getwd()


###############################################################

#sp <- "Calanus_helgolandicus"
#v <- "SST"
#omi <- "70p"


for(sp in sp.names) {
	
		# Useless message #1
		message(paste("Doing ",sp, sep="")) 
		#myRespName <- sp
		# Go to proper directory
		setwd( paste(second.wd,"/","pabs_files_matched","/", sep="") )
		DataSpecies <- read.table(paste("pabs_",sp,"_all.txt", sep=""), header = T)
			
		for(v in vars) {
		
				# Useless message #2 
				message(paste("Doing ",sp," & ",v, sep="")) 
				# Need to re-load DataSpecies with full predictors
				setwd( paste(second.wd,"/","pabs_files_matched","/", sep="") )
				DataSpecies <- read.table(paste("pabs_",sp,"_all.txt", sep=""), header = T)
				# restrict to variable of interest
				DataSpecies <- DataSpecies[,c("x","y","p",v)]
				
				# Omit data accrding to data omission degree ; but first, save the P and psAbs separately (you will sample randomly among each of them)
				Pres <- DataSpecies[which(DataSpecies$p == 1),]
				nPres <- nrow(Pres)
				psAbs <- DataSpecies[which(DataSpecies$p == 0),]
				npsAbs <- nrow(psAbs)
				
				# For each level of data omission, randomly sample among pres and psAbs and then concatenate them into a new dataset ("DataSpecies2")
				for(omi in omis) {
					
						if( omi == "50p") {
							
							newPres <- dplyr::sample_n(tbl = Pres, size = nPres/2 )
							newpsAbs <- dplyr::sample_n(tbl = psAbs, size = npsAbs/2 )
							# rbind
							DataSpecies2 <- rbind(newPres, newpsAbs)
							
						} else if (omi == "70p") {
							
							newPres <- dplyr::sample_n(tbl = Pres, size = nPres*0.7)
							newpsAbs <- dplyr::sample_n(tbl = psAbs, size = npsAbs*0.7)
							# rbind
							DataSpecies2 <- rbind(newPres, newpsAbs)
							
						} else if (omi == "80p") {
							
							newPres <- dplyr::sample_n(tbl = Pres, size = nPres*0.8)
							newpsAbs <- dplyr::sample_n(tbl = psAbs, size = npsAbs*0.8)
							# rbind
							DataSpecies2 <- rbind(newPres, newpsAbs)
							
						} # eo if else loop
									
						# Formating Data
						setwd(second.wd)
						require("biomod2")
						myRespCoord <- DataSpecies2[,c(1,2)]
						myExplVar <- DataSpecies2[,c(v)]
						myResp <- as.numeric(DataSpecies2$p)
						myRespName <- sp
						myRespName <- str_replace_all(myRespName, "_", ".") 
				
					
						# Define FUN to use for bootstrapping
						univariate_niche_model <- function(DataSpecies2, i) {
		
									DataSpecies2 <- DataSpecies2[i,]
		
									myRespCoord <- DataSpecies2[,c(1,2)]
									myExplVar <- DataSpecies2[,c(v)]
									myResp <- as.numeric(DataSpecies2$p)
									myRespName <- sp
									myRespName <- str_replace_all(myRespName, "_", ".") 	
			
									# Formatting Data
									setwd(second.wd)
									require("biomod2")
							
									# Modeling options :
									myBiomodOption2 <- BIOMOD_ModelingOptions( GAM = list( algo = 'GAM_mgcv',
							 									type = 's_smoother',
							 									k = 5,
							 									interaction.level = 0,
							 									myFormula = NULL,
							 									family = binomial("logit"),
							 									method = 'GCV.Cp',
							 									optimizer = c('outer','newton'),
							 									select = FALSE,
							 									knots = NULL,
							 									paraPen = NULL,
							 									control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
							 									, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
							 									, rank.tol = 1.49011611938477e-08
							 									, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
							 									, optim = list(factr=1e+07)
							 									, newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
							 									, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE)), )
			
									myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
				                                     expl.var = myExplVar,
				                                     resp.xy = myRespCoord,
				                                     resp.name = myRespName )
					
									# Niche modeling with GAMs
									myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
				               							models = c("GAM"), 
				               							models.options = myBiomodOption2, 
				               							NbRunEval = 5, 
				               	 						DataSplit = 80, 
				               							Prevalence = 0.5,
				               							VarImport = 0, 
				               							models.eval.meth = c('TSS'),
				               							SaveObj = FALSE,
				               							do.full.models = FALSE )
	   
									gc()

									scores <- data.frame(get_evaluations(myBiomodModelOut))
			
									# Retrieve TSS score, threshold approach
									if( length(which(is.na(scores))) != 0 ) {
											rm(scores)
											setwd(WD)
			
									} else {
			
											tss.scores.table <- data.frame(species = myRespName, 
													TSS_R1 = scores["TSS",c(1)], 
													TSS_R2 = scores["TSS",c(5)], 
													TSS_R3 = scores["TSS",c(9)],
													TSS_R4 = scores["TSS",c(13)],
													TSS_R5 = scores["TSS",c(17)] )

											# Compute average TSS:
											avgTSS <- rowMeans(tss.scores.table[,c(2:6)])  	

											# Load response plot
											response.plot <- BIOMOD_LoadModels(myBiomodModelOut, models= "GAM")

											# Get response plot from biomod2's special function
											myRespPlot2D <- response.plot2(models= response.plot,
														Data = get_formal_data(myBiomodModelOut,'expl.var'), 
														show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
														do.bivariate = FALSE,
														fixed.var.metric = 'median',
														col = c("purple", "blue", "red", "orange", "green"),
														legend = TRUE,
														data_species = get_formal_data(myBiomodModelOut,'resp.var') )

											# Transform into a data.frame for ggploting							   
											resp <- data.frame(myRespPlot2D)
											sdm <- "GAM"
											colnames(resp)[1:6] <- c(v, paste(v,sdm,"RUN1",sep="_"), paste(v,sdm,"RUN2",sep="_"), 
																	paste(v,sdm,"RUN3",sep="_"), paste(v,sdm,"RUN4",sep="_"), 
																	paste(v,sdm,"RUN5",sep="_") )
																	
											resp$avg_HSI <- rowMeans(resp[,c(2:6)], na.rm= F) # average across the 5 CV runs		


											### Now, need to compute univariate niche center and breadth (+ lower and upper limits for information):
											### If variable == Salinity, then only integrate from 30 psu 
											### (open ocean ; avoids spurious response curves which have high HSI at very low SSS levels)
			
											if( v == "SSS") {
				
												resp <- resp[which(resp[,v] >= 30 ), ]
				
												### a) Center = median
												resp$area <- cumsum(resp$avg_HSI)
												medarea <- (resp[length(resp$area),"area"] / 2)
												resp$dist1 <- abs(resp$area - medarea)
			
												center <- resp[which(resp$dist1 == min(resp$dist1)), c(v)]
												if( length(center) > 1) {
														center <- mean(center)
												} # eo if loop

												### b) Breadth = interquantile range between the 10th and the 90th
												q1 <- (resp[length(resp$area),"area"] / 10)
												q3 <- (9*(resp[length(resp$area),"area"]) / 10)

												resp$dist2 <- abs(resp$area - q1)
												resp$dist3 <- abs(resp$area - q3)

												lower <- resp[which(resp$dist2 == min(resp$dist2)), c(v)]
												if( length(lower) > 1) {
														lower <- mean(lower)
												} # eo if loop

												upper <- resp[which(resp$dist3 == min(resp$dist3)), c(v)]
												if( length(upper) > 1) {
														upper <- mean(upper)
												} # eo if loop

												breadth <- upper - lower
			
				
											# If v != SSS, then -------------------------------------------------------------------------------------------------------------
											} else {
			
												### a) Center = median
												resp$area <- cumsum(resp$avg_HSI)
												medarea <- (resp[length(resp$area),"area"] / 2)
												resp$dist1 <- abs(resp$area - medarea)
			
												center <- resp[which(resp$dist1 == min(resp$dist1)), c(v)]
												if( length(center) > 1) {
														center <- mean(center)
												}

												#plot <- ggplot(resp) + geom_path(aes(x= SST, y= avg_HSI), colour = "#3288bd", size = 1, linetype="solid") + theme_linedraw()
												#pdf(paste(myRespName,sdm,"SST_curve_k5.pdf", sep="_"), width = 6, height = 6)
												#plot
												#dev.off()

												### b) Breadth = interquantile range between the 10th and the 90th
												q1 <- (resp[length(resp$area),"area"] / 10)
												q3 <- (9*(resp[length(resp$area),"area"]) / 10)

												resp$dist2 <- abs(resp$area - q1)
												resp$dist3 <- abs(resp$area - q3)

												lower <- resp[which(resp$dist2 == min(resp$dist2)), c(v)]
													if( length(lower) > 1) {
														lower <- mean(lower)
												} # eo if loop

												upper <- resp[which(resp$dist3 == min(resp$dist3)), c(v)]
												if( length(upper) > 1) {
														upper <- mean(upper)
												} # eo if loop

												breadth <- upper - lower
			
											} # eo if else loop
				
											### For boot, need to return a vector containing all statistics
											return( c(avgTSS, center, breadth, lower, upper) )

									} # eo if else	
									
							} # eo FUN		
			
							myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
									            	 expl.var = myExplVar,
									            	 resp.xy = myRespCoord,
									             	resp.name = myRespName )
													
							myBiomodOption2 <- BIOMOD_ModelingOptions( GAM = list( algo = 'GAM_mgcv',
											 						type = 's_smoother',
											 						k = 5,
											 						interaction.level = 0,
											 						myFormula = NULL,
											 						family = binomial("logit"),
											 						method = 'GCV.Cp',
											 						optimizer = c('outer','newton'),
											 						select = FALSE,
											 						knots = NULL,
											 						paraPen = NULL,
											 						control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
											 						, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
											 						, rank.tol = 1.49011611938477e-08
											 						, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
											 						, optim = list(factr=1e+07)
											 						, newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
											 						, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE)), )
					
							# Niche modeling with GAMs
							myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
									            	 models = c("GAM"), 
									            	 models.options = myBiomodOption2, 
									             	NbRunEval = 5, 
									             	  DataSplit = 80, 
									             	 Prevalence = 0.5,
									             	VarImport = 0, 
									             	models.eval.meth = c('TSS'),
									             	SaveObj = FALSE,
									             	do.full.models = FALSE )
	   
	   					
							# Load response plot
							response.plot <- BIOMOD_LoadModels(myBiomodModelOut, models= "GAM")

							# Get response plot from biomod2's special function
							myRespPlot2D <- response.plot2(models = response.plot,
								 						Data = get_formal_data(myBiomodModelOut,'expl.var'), 
								 						show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
								 						do.bivariate = FALSE,
								 						fixed.var.metric = 'median',
								 						col = c("purple", "blue", "red", "orange", "green"),
								 						legend = TRUE,
								 						data_species = get_formal_data(myBiomodModelOut, 'resp.var') )

							# Transform into a data.frame for ggploting							   
							resp <- data.frame(myRespPlot2D)
							sdm <- "GAM"
							colnames(resp)[1:6] <- c(v, paste(v,sdm,"RUN1", sep= "_"), 
														paste(v,sdm,"RUN2", sep= "_"), 
														paste(v,sdm,"RUN3", sep= "_"), 
														paste(v,sdm,"RUN4", sep= "_"), 
														paste(v,sdm,"RUN5", sep= "_") )
														
							resp$avg_HSI <- rowMeans(resp[,c(2:6)], na.rm= F) # average across the 5 CV runs	
				
							rm(resp, sdm, myRespPlot2D, response.plot)
			
							require("boot")
							myboot <- boot(data = DataSpecies2, statistic = univariate_niche_model, R = 30) # parallel= "multicore", ncpus= 12)
		
							# Turn into a dataframe: t1 = average TSS, t2 = niche center, t3 = niche breadth, t4 = lower niche boundary, t5 = upper niche boundary
							# str(myboot)
		
							# myboot$t0
							stats <- data.frame(sp.name = sp, 
												avgTSS = myboot$t0[1], sdTSS = sd(myboot$t[,1]), 
												center = myboot$t0[2], sdcenter = sd(myboot$t[,2]), 
												breadth = myboot$t0[3], sdbreadth = sd(myboot$t[,3]), 
												lower = myboot$t0[4], sdlower = sd(myboot$t[,4]), 
												upper = myboot$t0[5], sdupper = sd(myboot$t[,5]) )
			
		
							# Got tonproper dir and save 'stats' object
							if( omi == "50p") {
								setwd(paste(second.wd, "/", "niche_stats_16_01_17", "/","omission_test","/",omi,"/", sep= ""))
							} else if (omi == "70p") {
								setwd(paste(second.wd, "/", "niche_stats_16_01_17", "/","omission_test","/",omi,"/", sep= ""))
							} else if (omi == "80p") {
				    			setwd(paste(second.wd, "/", "niche_stats_16_01_17", "/","omission_test","/",omi,"/", sep= ""))	
							} # eo if else loop
				
							# Save stats in dir
							save(stats, file= paste(sp,"_",v,".Rdata", sep=""))
							gc()
							rm(stats, myboot)
						
				} # eo omi in omis


		} # eo v in vars


} # eo sp in sp.names



###############################################################################################################################################################################

### 22/02/2017 : For each test (2 ENMs; 3 k levels; 3 omission levels)

setwd(paste(WD,"/","niche.modelling","/","pabs_files_matched","/", sep=""))
sp.names <- dir()
sp.names <- str_replace_all(sp.names, "pabs_", "")
sp.names <- str_replace_all(sp.names, "_all.txt", "")
sp.names


# ENMs : 
setwd(paste(second.wd,"/","niche_stats_16_01_17","/","ENM_test","/","MAXENT","/", sep = ""))
dir()
#sp <- "Parapontella_brevicornis"

# For testing first: 
for(sp in sp.names) {
	
		message(paste("Doing ", sp, sep = ""))
		
		sst <- get(load(paste(sp,"_","SST",".Rdata", sep="")))
		dsst <- get(load(paste(sp,"_","dSST",".Rdata", sep="")))
		sss <- get(load(paste(sp,"_","SSS",".Rdata", sep="")))
		mld <- get(load(paste(sp,"_","MLD1",".Rdata", sep="")))
		mlpar <- get(load(paste(sp,"_","MLPAR1",".Rdata", sep="")))
		logChla <- get(load(paste(sp,"_","logChla",".Rdata", sep="")))
		
		# Cbind
		table <- cbind(sst, dsst, sss, mld, mlpar, logChla)
}

### Try to fox the following BRTs : 
# Candacia varicans # 17 ; logChla model > supply NAs          	# DONE
# Euaugaptilus_hecticus # 37 ; logChla models --> supply NAs   	# DONE
# Oithona nana # 70 ; logChla model --> supply NAs				# DONE
# Parapontella_brevicornis # 85 ; MLPAR1 model --> supply NAs	# DONE & DONE for MAXENT as well
# Phaenna spinifera # 87 ; logChla model --> supply NAs			# DONE
# Sapphirina nigromaculata # 94 ; logChla model --> supply NAs	# 


### 23/02/2017 : supplying NAs
setwd(paste(second.wd,"/","niche_stats_16_01_17","/","ENM_test","/","MAXENT","/", sep = ""))
d <- data.frame(sp.name = "Parapontella_brevicornis", avgTSS = NA, sdTSS = NA, center = NA, sdcenter = NA, breadth = NA, sdbreadth = NA, lower = NA, sdlower = NA, upper = NA, sdupper = NA)
d
save(d, file = paste("Parapontella_brevicornis","_","MLPAR1",".Rdata", sep = ""))



### 23/02/2017: For each sensitivity tests, concatenate all model outputs and save the resulting dataframe as a text file
# MAXENT : done
# BRT : done

# k3 : done
# k10 : done
# k15 : done

# 80p : done
# 70p : 
# 50p : 

# Proper dir
setwd(paste(second.wd,"/","niche_stats_16_01_17","/","ENM_test","/","BRT","/", sep = ""))

res <- lapply(sp.names, function(sp) {

				# Load SST, SSS etc.
				sst <- get(load(paste(sp,"_","SST",".Rdata", sep="")))
				dsst <- get(load(paste(sp,"_","dSST",".Rdata", sep="")))
				sss <- get(load(paste(sp,"_","SSS",".Rdata", sep="")))
				mld <- get(load(paste(sp,"_","MLD1",".Rdata", sep="")))
				mlpar <- get(load(paste(sp,"_","MLPAR1",".Rdata", sep="")))
				logChla <- get(load(paste(sp,"_","logChla",".Rdata", sep="")))
				
				# Cbind
				table <- cbind(sst, dsst, sss, mld, mlpar, logChla)
				# Re-name columns
				table <- table[,-c(12,23,34,45,56)]
				
				colnames(table)[2:11] <- c("avgTSS_SST","sdTSS_SST","center_SST","sdcenter_SST","breadth_SST",
											"sdbreadth_SST","lower_SST","sdlower_SST","upper_SST","sdupper_SST")
				
				colnames(table)[grep(".1", colnames(table))] <- c("avgTSS_dSST","sdTSS_dSST","center_dSST","sdcenter_dSST","breadth_dSST",
																  "sdbreadth_dSST","lower_dSST","sdlower_dSST","upper_dSST","sdupper_dSST")
																  
				colnames(table)[grep(".2", colnames(table))] <- c("avgTSS_SSS","sdTSS_SSS","center_SSS","sdcenter_SSS","breadth_SSS",
												  				"sdbreadth_SSS","lower_SSS","sdlower_SSS","upper_SSS","sdupper_SSS")
																
				colnames(table)[grep(".3", colnames(table))] <- c("avgTSS_MLD1","sdTSS_MLD1","center_MLD1","sdcenter_MLD1","breadth_MLD1",
																"sdbreadth_MLD1","lower_MLD1","sdlower_MLD1","upper_MLD1","sdupper_MLD1")
																
				colnames(table)[grep(".4", colnames(table))] <- c("avgTSS_MLPAR1","sdTSS_MLPAR1","center_MLPAR1","sdcenter_MLPAR1","breadth_MLPAR1",
																"sdbreadth_MLPAR1","lower_MLPAR1","sdlower_MLPAR1","upper_MLPAR1","sdupper_MLPAR1")	
				
				colnames(table)[grep(".5", colnames(table))] <- c("avgTSS_logChla","sdTSS_logChla","center_logChla","sdcenter_logChla","breadth_logChla",
																"sdbreadth_logChla","lower_logChla","sdlower_logChla","upper_logChla","sdupper_logChla")
																
				return(table)																																																								  

} )

data <- do.call(rbind, res)
dim(data) # 106  61
data
summary(data)

setwd(WD)
write.table(file = "niche_traits_BRT_23_02_17.txt", x = data, sep = ";", dec = ".")

# t <- read.table("niche_traits_MAXENT_23_02_17.txt", h = T, sep = ";", dec = ".")
