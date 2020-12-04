
##### 17/01/2017 - LOV - Fabio Benedetti
##### Script for : 
#				- conducting the sensitivity analyses prior to manucsript submission to L&O
#				- assessing niche traits variations due to ENM choice
#				- compute niche traits not from GAMs, but from : GLMs, MAXENT and/ or BRT 
#				- plot the SST response curves from diverse models

### Last update : 28/02/2017

# --------------------------------------------------------------------------------------------------------------------------------------------

library("PresenceAbsence")
library("biomod2")
#library("raster")
#library("rgdal")
library("sp")
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

# List of species names:é
setwd(paste(WD,"/","niche.modelling","/","pabs_files_matched","/", sep=""))
sp.names <- dir()
sp.names <- str_replace_all(sp.names, "pabs_", "")
sp.names <- str_replace_all(sp.names, "_all.txt", "")
#sp.names

### 19/01/2017 : C.plumatus logChla models tend to fail
#sp.names <- sp.names[103:106]

# List of ENM choices
#SDMs <- c("GLM","GBM")

setwd(WD)
setwd("./niche.modelling")
second.wd <- getwd()


###############################################################

# For single modelling:
sp <- "Clausocalanus_mastigophorus"
v <- "SST"
sdm <- "ANN"


for(sp in sp.names) {
	
		# Useless message #1
		message(paste("Doing ",sp, sep="")) 
		#myRespName <- sp
		# Go to proper directory
		setwd( paste(second.wd,"/","pabs_files_matched","/", sep = "") )
		DataSpecies <- read.table(paste("pabs_",sp,"_all.txt", sep = ""), header = T)
			
		for(v in vars) {
		
				# Useless message #2 
				message(paste("Doing ",sp," & ",v, sep="")) 
				# Need to re-load DataSpecies with full predictors
				setwd( paste(second.wd,"/","pabs_files_matched","/", sep="") )
				DataSpecies <- read.table(paste("pabs_",sp,"_all.txt", sep=""), header = T)
				# restrict to variable of interest
				DataSpecies <- DataSpecies[,c("x","y","p",v)]
									
				# Formating Data
				setwd(second.wd)
				require("biomod2")
				myRespCoord <- DataSpecies[,c(1,2)]
				myExplVar <- DataSpecies[,c(v)]
				myResp <- as.numeric(DataSpecies$p)
				myRespName <- sp
				myRespName <- str_replace_all(myRespName, "_", ".") 
				
				# For each smoothing parameter
				for(sdm in SDMs) {
					
						# Define FUN to use for bootstrapping
						univariate_niche_model <- function(DataSpecies, i) {
									
									# Formatting Data
									DataSpecies <- DataSpecies[i,]
									myRespCoord <- DataSpecies[,c(1,2)]
									myExplVar <- DataSpecies[,c(v)]
									myResp <- as.numeric(DataSpecies$p)
									myRespName <- sp
									myRespName <- str_replace_all(myRespName, "_", ".") 	
			
									# Formatting Data
									setwd(second.wd)
									require("biomod2")
							
									# Provide modelling options for GLMs, MAXENT and BRTs	
									myBiomodOption2 <- BIOMOD_ModelingOptions( MAXENT = list( path_to_maxent.jar = second.wd,
																		maximumiterations = 500,
																		visible = FALSE,
																		linear = TRUE, 
																		quadratic = TRUE,
																		product = TRUE,
																		threshold = FALSE,
																		hinge = TRUE,
																		lq2lqptthreshold = 80,
																		l2lqthreshold = 10,
																		hingethreshold = 15,
																		beta_threshold = -1,
																		beta_categorical = -1,
																		beta_lqp = -1,
																		beta_hinge = -1,
																		defaultprevalence = 0.5),
							
																	 GBM = list( distribution = 'bernoulli',
																	 	n.trees = 2000,
																	  	interaction.depth = 1,
																	  	n.minobsinnode = 5,
																	  	shrinkage = 0.001,
																	  	bag.fraction = 0.5,
																	  	train.fraction = 1,
																	  	cv.folds = 0,
																	  	keep.data = FALSE,
																	  	verbose = FALSE,
																	  	perf.method = "OOB"), 
																		
																	 GLM = list( type = 'quadratic',
													           		 	interaction.level = 0,
													           			myFormula = NULL,
													           			test = 'AIC',
													            		family = binomial("logit"),
													            		mustart = 0.5,
													            		control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE))
																			  
														) # eo BIOMOD_ModelingOptions
			
									myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
				                                     expl.var = myExplVar,
				                                     resp.xy = myRespCoord,
				                                     resp.name = myRespName )
					
									# Niche modeling with GAMs
									myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
				               							models = sdm, 
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
											response.plot <- BIOMOD_LoadModels(myBiomodModelOut, models = sdm)

											# Get response plot from biomod2's special function
											myRespPlot2D <- response.plot2(models= response.plot,
														Data= get_formal_data(myBiomodModelOut,'expl.var'), 
														show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
														do.bivariate= FALSE,
														fixed.var.metric= 'median',
														col= c("purple", "blue", "red", "orange", "green"),
														legend= TRUE,
														data_species= get_formal_data(myBiomodModelOut,'resp.var') )

											# Transform into a data.frame for ggploting							   
											resp <- data.frame(myRespPlot2D)
											colnames(resp)[1:6] <- c(v, paste(v,sdm,"RUN1",sep="_"), paste(v,sdm,"RUN2",sep="_"), 
																	paste(v,sdm,"RUN3",sep="_"), paste(v,sdm,"RUN4",sep="_"), 
																	paste(v,sdm,"RUN5",sep="_") )
																	
											resp$avg_HSI <- rowMeans(resp[,c(2:6)], na.rm = F) # average across the 5 CV runs		


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
									
							} # eo FUN univariate_niche_model for boostrapping		
			
							# Data formatting	
							myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
									            	 expl.var = myExplVar,
									            	 resp.xy = myRespCoord,
									             	resp.name = myRespName )
													
							# Provide modelling options for GLMs, MAXENT and BRTs						
							myBiomodOption2 <- BIOMOD_ModelingOptions( MAXENT = list( path_to_maxent.jar = second.wd,
																		maximumiterations = 500,
																		visible = FALSE,
																		linear = TRUE, 
																		quadratic = TRUE,
																		product = TRUE,
																		threshold = FALSE,
																		hinge = TRUE,
																		lq2lqptthreshold = 80,
																		l2lqthreshold = 10,
																		hingethreshold = 15,
																		beta_threshold = -1,
																		beta_categorical = -1,
																		beta_lqp = -1,
																		beta_hinge = -1,
																		defaultprevalence = 0.5),
							
																	 GBM = list( distribution = 'bernoulli',
																	 	n.trees = 2500,
																	  	interaction.depth = 1,
																	  	n.minobsinnode = 5,
																	  	shrinkage = 0.001,
																	  	bag.fraction = 0.5,
																	  	train.fraction = 1,
																	  	cv.folds = 0,
																	  	keep.data = FALSE,
																	  	verbose = FALSE,
																	  	perf.method = "OOB"), 
																		
   																	 GLM = list( type = 'quadratic',
   													           		 	interaction.level = 0,
   													           			myFormula = NULL,
   													           			test = 'AIC',
   													            		family = binomial("logit"),
   													            		mustart = 0.5,
   													            		control = glm.control(epsilon = 1e-08, maxit = 100, trace = FALSE)) 
																		
												) # BIOMOD_ModelingOptions
					
							# Niche modeling with chosen sdm
							myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
									            	 models = sdm, 
									            	 models.options = myBiomodOption2, 
									             	 NbRunEval = 5, 
									             	 DataSplit = 80, 
									             	 Prevalence = 0.5,
									             	 VarImport = 0, 
									             	 models.eval.meth = c('TSS'),
									             	 SaveObj = FALSE,
									             	 do.full.models = FALSE )   
	   					
							# Load response plot
							response.plot <- BIOMOD_LoadModels(myBiomodModelOut, models = sdm)

							# Get response plot from biomod2's special function
							myRespPlot2D <- response.plot2(models = response.plot,
								 						Data = get_formal_data(myBiomodModelOut,'expl.var'), 
								 						show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
								 						do.bivariate = FALSE,
								 						fixed.var.metri = 'median',
								 						col = c("purple", "blue", "red", "orange", "green"),
								 						legend = TRUE,
								 						data_species = get_formal_data(myBiomodModelOut,'resp.var') )

							# Transform into a data.frame for ggploting							   
							resp <- data.frame(myRespPlot2D)
							# dim(resp)
							
							colnames(resp)[1:6] <- c(v, paste(v,sdm,"RUN1", sep = "_"), 
														paste(v,sdm,"RUN2", sep = "_"), 
														paste(v,sdm,"RUN3", sep = "_"), 
														paste(v,sdm,"RUN4", sep = "_"), 
														paste(v,sdm,"RUN5", sep = "_") )
																								
							resp$avg_HSI <- rowMeans(resp[,c(2:6)], na.rm = F) # average across the 5 CV runs	
							
							### 28/02/2017 : plot SST response curves for C. helgolandicus, according to 9 ENMs
							plot_resp_curves <- ggplot(data = resp) + 
												geom_path(aes(x = SST, y = SST_ANN_RUN1), linetype = "solid", colour = "#d53e4f") + 
												geom_path(aes(x = SST, y = SST_ANN_RUN2), linetype = "longdash", colour = "#fdae61") + 
												geom_path(aes(x = SST, y = SST_ANN_RUN3), linetype = "dashed", colour = "#66c2a5") + 
												geom_path(aes(x = SST, y = SST_ANN_RUN4), linetype = "twodash", colour = "#3288bd") + 
												geom_path(aes(x = SST, y = SST_ANN_RUN5), linetype = "solid", colour = "#5e4fa2") + 
												xlab("Sea Surface Temperature (°C)") + ylab("Habitat Suitability Index (HSI)") + theme_bw()
							
							quartz()
							plot_resp_curves
							# Save plot
							setwd(WD)
							ggsave(filename = "resp_curve_SST_Clausocalanus_mastigophorus_ANN.pdf", plot = plot_resp_curves, dpi = 300, height = 7, width = 7)	
				
							rm(resp, myRespPlot2D, response.plot)
			
							require("boot")
							myboot <- boot(data = DataSpecies, statistic = univariate_niche_model, R = 20) #parallel = "multicore", ncpus = 10)
		
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
							if( sdm == "GLM") {
								setwd(paste(second.wd, "/", "niche_stats_16_01_17", "/","ENM_test","/","GLM","/", sep = ""))
							} else if (sdm == "MAXENT") {
								setwd(paste(second.wd, "/", "niche_stats_16_01_17", "/","ENM_test","/","MAXENT","/", sep = ""))
							} else if (sdm == "GBM") {
							    setwd(paste(second.wd, "/", "niche_stats_16_01_17", "/","ENM_test","/","BRT","/", sep = ""))	
							} # eo if elseoop
				
							# Save stats in dir
							save(stats, file = paste(sp, "_", v, ".Rdata", sep = ""))
							rm(stats, myboot)
						
				} # eo sdm in SDMs


		} # eo v in vars


} # eo sp in sp.names




