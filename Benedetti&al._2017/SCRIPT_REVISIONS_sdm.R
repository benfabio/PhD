#
##### 10/03/2015 - LOV - Fabio Benedetti, Guilhem Marre, François Guilhaumon
##### Script for : 
#					- Defining the function ResNicheModelling function that will launch the biomod2's SDMs on all species, and project models for each
#					  pseudo-absence run, each period (T0, T1, T2), each scenario 
#					- Prints models' evaluation scores (TSS & AUC), and all models' outputs in separate folders
#					- Save training & testing datasets (with coordinates)
#					- Contrary to Guilhem's script, probabilities are not converted to presences and absences again ; we are are keeping raw probabilities 
#					  (are SDMs rescaled between 0 and 1000 ; will be trescaled to 0-1 when extracting ouputs)
#					- apply function in mclapply or dplyr for parallel computation

 
### Last update : 09/05/2016
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("rJava")
library("rgdal")
library("rgeos")
library("ggplot2")
library("stringr")
library("foreign")
library("PresenceAbsence")
library("gtools")
library("SDMTools")
library("biomod2")
library("parallel")

#old_gbm <- "http://cran.r-project.org/src/contrib/Archive/gbm/gbm_2.0-8.tar.gz"
#install.packages(old_gbm, repos=NULL, type="source")
sessionInfo()

library("gbm")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Loading environmental predictors, present and future Mediterranean fields, and species' presences/ pseudoAbsences files
# Tu charges ici tes variables environnementales
# Pour calibration : climatos WOA13, globales, 1965-1994
varT0_glob <- brick("clims_woa13_baseline_04.grd")
varT0_glob <- dropLayer(varT0_glob, c(1,2,6))
# varT0_glob

# Pour projection/prédiction : climatos WOA13 Mediterranée et NM8
varT0 <- read.table("climatos_woa13_6594_Levin.txt", header=T, sep="\t", dec=".")
# Check avant
# summary(varT0)
# ggplot(varT0) + geom_point(aes(x=x, y=y, colour = avgSST)) + theme_bw() + coord_map()

# On crée le dossier dans lequel on va stocker tous les résultats de la modélo et on y dépose gentiment le fichier java de maxente
initial.wd <- getwd()
setwd("./niche_modelling")
second.wd <- getwd()

# Liste des scenarios
scenario <- c("A2","A2F","A2RF","A2ARF","A1BARF","B1ARF")
# List des runs de pseudoAbs
runs <- c(1:10)
# ENMs list:
SDMs <- c('SRE','CTA','RF','GBM','GLM','GAM','MARS','FDA','ANN','MAXENT')
# List of eval_runs:
eval_runs <- c("RUN1","RUN2","RUN3") ## !! need an underscore because biomod2 adds one by default in pred and calib_lines !! 
# List of prevalence levels:
prevalences <- c("n1","n1.5","n2","n10","n50")

# Species names list:
setwd(paste(second.wd,"/","pabs_files","/", sep=""))
sp.names <- dir() ; length(sp.names)
setwd(second.wd)


################################################################# ¡ For testing !

# sp.n <- 'Calanus_helgolandicus'
# r <- 1
# s <- "A2"
# preval <- "n2"
# splitValidation <- 80

################################################################################

##### La fonction à appliquer :

SpeciesNicheModelling <- function(sp.n = sp.names, 
							  splitValidation = 80, 
                              sdm.models = SDMs ) {
							  
				# Pour chaque niveau de prévalence
				for(preval in prevalences) {
				
				# Pour chaque run de pseudo-absence :	
				for(r in runs) {

							############################# 1°) Préparer les data
  							myRespName <- sp.n
  
  			  				#cat('\n', myRespName, 'modeling...')
  
  			  				### Définition des data
							setwd( paste(second.wd,"/","pabs_files","/",myRespName,"/", sep="") )
							DataSpecies <- read.csv(paste("pabs_",myRespName,"_",preval,"_run",r,".txt", sep=""), header = T, sep="\t", dec=".")
							setwd(second.wd)
							
							# Charger les présences/pseudoabs ainsi que leurs coordonnées respectives
  			  				myResp <- as.numeric(DataSpecies[,3])
  			  				myRespCoord <- DataSpecies[,c(1,2)]
  
  			  				### Initialisation: formatage des données
							# Changer le nom de l'espèce afin de conserver le format BIOMOD, et ainsi éviter des problèmes d'arborescence
							myRespName <- str_replace_all(myRespName, "_", ".")
							
							# formatage des données
  			  				myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                      	 			 	expl.var = varT0_glob,
                                       			  	 	resp.xy = myRespCoord,
                                       			  	 	resp.name = myRespName )
  
  							### Définition des options (default options)
							### More permissive SRE
							### No cv.folds in GBM, may cause issues with socketConnections when mclapplying
							#   https://r-forge.r-project.org/forum/forum.php?thread_id=28450&forum_id=995&group_id=302
							#   http://cran.r-project.org/web/packages/gbm/gbm.pdf
							
  			  				myBiomodOption <- BIOMOD_ModelingOptions(SRE = list(quant = 0.010), 
							
																	MAXENT = list( path_to_maxent.jar = second.wd,
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
																		
																	GAM = list(algo = 'GAM_mgcv',
																		type = 's_smoother',
																		k = 4,
																		interaction.level = 0,
																		myFormula = NULL,
																		family = binomial("logit"),
																		method = 'GCV.Cp',
																		optimizer = c('outer','newton'),
																		select = FALSE,
																		knots = NULL,
																		paraPen = NULL,
																		control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
																		,maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
																		,rank.tol = 1.49011611938477e-08
																		,nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
																		,optim = list(factr=1e+07)
																		,newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
																		,outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE)),	
							
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
																	  	perf.method = "OOB") )
  
							############################# 2°) Modéliser l'habitat suitability des espèces sur les climatos WOA13 (65-94)
							
 			   				### Modelling
							# library("gbm")
							
  			  				myBiomodModelOut <- BIOMOD_Modeling(
    										myBiomodData,
    										models = SDMs,
    										models.options = myBiomodOption,
    										NbRunEval = length(eval_runs),
    										DataSplit = splitValidation,
    										Prevalence = 0.5, # this is to ensure presence ans pseudoabsences are weighted equally
    										VarImport = 0,
    										models.eval.meth = c('TSS','ROC'),
    										SaveObj = TRUE,
    										rescal.all.models = TRUE,
    										do.full.models = FALSE,
    										modeling.id = paste(myRespName,preval,"run",r, sep = "_") )						
  
  							### On sauve les statistiques d'évaluation du modèle
  			  				capture.output( get_evaluations(myBiomodModelOut), 
											file = file.path(myRespName, paste(myRespName,"_formal_models_evaluation_",preval,"_run",r,".txt", sep="") ) )
							
							gc()
							
							############################# 3°) Comment récupérer les traing et testing sets pour calculer le TSSmax des tirages binomiaux
							### https://r-forge.r-project.org/forum/forum.php?thread_id=31252&forum_id=995&group_id=302
							### © Damien Georges
							pred <- get_predictions(myBiomodModelOut)
							calib_lines <- get_calib_lines(myBiomodModelOut) ## coded as TRUE: calib pt ; FALSE: eval pt; NA: non used pt (not in the PA dataset)
							# Get presences-absences coordinates : 
							resp_with_PA <- myBiomodData@data.species ## coded as 1: presence ; 0: absence ; NA: pseudo-absence
							coordinates_with_PA <- myBiomodData@coord 
							
							### Faire le reste dans une for loop : une fois pour chgaque eval_runs (n = 3; R)
							for(R in eval_runs) {
							
											# Corresponding points	
											run_pa_pts <- which(!is.na(calib_lines[,paste("_",R,sep=""),"_AllData"])) # Watchout for the added underscore
											run_pa_coord <- coordinates_with_PA[run_pa_pts,]
											head(run_pa_coord)
											### a ddf containing the points' coordinates ! 
							
											run_pa_resp <- resp_with_PA[run_pa_pts]
											### Whether point is a presence (1) or an absence (0)

											# So... for the cabib/validation points ofthe GLM and RF for the RUN1 and PA1 you can build
											# a kind of summary table like.. © D. Georges
											# data_run1_pa1 <- data.frame(run1_pa1_coord, ## the coordinates
											# response = run1_pa1_resp, ## response points associated
											# is_a_calib_pt = na.omit(calib_lines[,,"_AllData"]) )
											# Or :
											# data_run1_pa1 <- data.frame(run1_pa1_coord,
														#response = run1_pa1_resp,
														#is_a_calib_pt = na.omit(calib_lines[,"_RUN1","_AllData"]),
														#pred_glm_run1_pa1 = pred[,"GLM","RUN1","AllData"], 
														#pred_rf_run1_pa1 = pred[,"RF","RUN1","AllData"])
							
											### For each eval_run, create a data.frame where you will store : 
											# - calib_points coordinates
											# - whether they are presences or absence (response)
											# - whether they were used to calibrate or test the models
											# - the SDMs' predicted HSI values (0-1000)
											### ¡! watch out, some SDMs may fail (especially GBMs), therefore create pred_HSI columns accroding to 'pred' dimensions ¡!
											# Save data for eval_run #1
											data_run <- data.frame(x = run_pa_coord$x, 
																   y = run_pa_coord$y,
																   response = run_pa_resp,
																   is_calib = na.omit(calib_lines[,paste("_",R,sep=""),"_AllData"]),
																   pred_sre = pred[,"SRE",R,"AllData"],
																   pred_cta = pred[,"CTA",R,"AllData"],
																   pred_rf = pred[,"RF",R,"AllData"],
																   pred_gbm = pred[,"GBM",R,"AllData"],
																   pred_glm = pred[,"GLM",R,"AllData"],
																   pred_gam = pred[,"GAM",R,"AllData"],
																   pred_mars = pred[,"MARS",R,"AllData"],
																   pred_fda = pred[,"FDA",R,"AllData"],
																   pred_ann = pred[,"ANN",R,"AllData"],
																   pred_maxent = pred[,"MAXENT",R,"AllData"] ) 
													
											### Go to species directory, create a new folder for printing the calib data above
											### and save the 3 dataframes as Rdata : 
											setwd(paste(second.wd,"/",myRespName,"/", sep=""))
											if(R == "RUN1") {
													dir.create("eval_runs_data")
							   			 			setwd(paste(second.wd,"/",myRespName,"/","eval_runs_data","/", sep=""))
											} else {
							   			 			setwd(paste(second.wd,"/",myRespName,"/","eval_runs_data","/", sep=""))
											}
											# Save as Rdata
											save(data_run, file = paste("Data_eval","_",preval,"_run",r,"_",R,".Rdata", sep="") )
											# go back to second.wd
											setwd(second.wd)
							
								} # eo in R in eval_runs
								
						    ### Check content : 
							# setwd(paste(second.wd,"/",myRespName,"/","DATA_FOR_EVALRUNS","/", sep=""))
							# file <- load("Data_eval_RUN9_run1.Rdata")
							# file <- get(file)
							# summary(file)
							# head(file)
							# Seems ok !
							
							############################# 4°) Ensemble modelling :
  
 			   				### Building ensemble-models
  			  				myBiomodEM <- BIOMOD_EnsembleModeling(
    										modeling.output = myBiomodModelOut,
    										chosen.models = 'all',
    										em.by = 'all',  						
    										eval.metric = c('TSS'),
    										eval.metric.quality.threshold = c(0.0), 
    										prob.mean = TRUE,
    										prob.cv = FALSE,
    										prob.ci = FALSE,
    										#prob.ci.alpha = 0.05,
    										prob.median = FALSE,
    										committee.averaging = FALSE,
    										prob.mean.weight = FALSE,
    										prob.mean.weight.decay = FALSE)
											
  
  			  				### Make projections on current variable / T0
  			  				myBiomodProj <- BIOMOD_Projection(
    											modeling.output = myBiomodModelOut,
    											new.env = varT0[,3:5],
    											proj.name = paste(myRespName,"_",preval,"_run_",r,"_","T0", sep=""),
    											selected.models = 'all',
    											binary.meth = 'TSS',
    											compress = 'xz',
    											clamping.mask = FALSE )
							
  			  				### Make projections on current global environment : T0_glob
							### ¡! Takes most of the time ¡!
  			  				myBiomodProj <- BIOMOD_Projection(
    											modeling.output = myBiomodModelOut,
    											new.env = varT0_glob,
    											proj.name = paste(myRespName,"_",preval,"_run_",r,"_","T0_glob", sep=""),
    											selected.models = 'all',
    											binary.meth = 'TSS',
    											compress = 'xz',
    											clamping.mask = FALSE )
  
		
							############################# 5°) Projection des distributions pour T1 (2020-2049) : 
							
							#### Pour chaque scénario d'émission SRES (Adloff & al., 2015) :					  
							for(s in scenario) {
								
										setwd(initial.wd) # aller au WD de base ou sont storés les climatos futures
							
										if(s == "A2") { 
													vars.envT1 <- read.table("clims_A2_2050_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else if (s == "A2F"){ 
													vars.envT1 <- read.table("clims_A2F_2050_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else if(s == "A2RF"){ 
													vars.envT1 <- read.table("clims_A2RF_2050_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else if(s == "A2ARF"){ 
													vars.envT1 <- read.table("clims_A2ARF_2050_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else if(s == "A1BARF"){ 
													vars.envT1 <- read.table("clims_A1BARF_2050_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else { 
													vars.envT1 <- read.table("clims_B1ARF_2050_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
										} # eo else if
										
										setwd(second.wd) # revenir au second WD
				
										### Projection à T1
  			  							myBiomodProj <- BIOMOD_Projection(
    		  												modeling.output = myBiomodModelOut,
    														new.env = vars.envT1[,3:5],
    														proj.name = paste(myRespName,"_",preval,"_run_",r,"_","T1_",s, sep=""),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE)
  
  			  							### Make ensemble-models projections on T1 variables
  			  							myBiomodEF <- BIOMOD_EnsembleForecasting(
    		  												projection.output = myBiomodProj,
    														EM.output = myBiomodEM)
															
										rm(vars.envT1) # retirer les variables de T1
										gc()					

										# -------------------------------------------------------------------------------------------------------

										### Projection à T2
										setwd(initial.wd) # aller au WD de base ou sont storés les climatos futures
										
										if(s == "A2") { 
													vars.envT2 <- read.table("clims_A2_6998_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else if(s == "A2F"){ 
													vars.envT2 <- read.table("clims_A2F_6998_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else if(s == "A2RF"){ 
													vars.envT2 <- read.table("clims_A2RF_6998_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else if(s == "A2ARF"){ 
													vars.envT2 <- read.table("clims_A2ARF_6998_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else if(s == "A1BARF"){ 
													vars.envT2 <- read.table("clims_A1BARF_6998_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
												} else { 
													vars.envT2 <- read.table("clims_B1ARF_6998_final_Levin.txt", header=TRUE, sep="\t", dec=".") 
										} #eo if
										
										setwd(second.wd) # revenir au second WD
				
  			  			  				myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = vars.envT2[,3:5],
    														proj.name = paste(myRespName,"_",preval,"_run_",r,"_","T2_",s, sep=""),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE )
  
  										### Make ensemble-models projections on T3 variables
  			  							myBiomodEF <- BIOMOD_EnsembleForecasting(
    													projection.output = myBiomodProj,
    													EM.output = myBiomodEM )
										
										# Clear				
										gc()
											
							} #eo s in scenario
  
					} #eo r in runs
					
			} # eo preval in prevalences			
  
} #eo NicheModelling


### Parallel Computating of Niche Modelling

### Need to re-launch (because of brackets in species names):
### To find them quickly:
# sp.names <- str_replace_all(sp.names, "_", ".")
# species[match(sp.names, species)]
 
#	- 3 Acartia spp: Acartia_Acartia_danae, Acartia_Acartia_negligens, Acartia_Acartiura_clausi
#	- 3 Corycaeus spp: Corycaeus_Agetus_flaccus, Corycaeus_Agetus_limbatus, Corycaeus_Agetus_typicus
#	- 2 Temora spp: Temora_stylifera, Temora_longicornis
#	- Subeucalanus_monachus
#	- Rhincalanus_nasutus
#	- Pontellina_plumata
#	- Pseudocalanus_elongatus
#	- 3 Paracalanus: Paracalanus_denudatus, Paracalanus_nanus, Paracalanus_parvus
#	- Oithona_similis, Oithona_setigera and Oithona_plumifera
#	- Microsetella_rosea
#	- Monothula_subtilis
#	- Mormonilla_phasma
#	- Labidocera_wollastoni
#	- Lubbockia_squillimana
#	- Lucicutia_clausi
#	- Farranula_rostrata
#	- Euterpina_acutifrons
#	- Euchirella_rostrata
#	- Chiridius_poppei
#	- Clausocalanus_furcatus ; Clausocalanus_arcuicornis

# Total = 30/30

### 09/05/2016: When re-running some missing models projections
sp.names <- c("Calanus_helgolandicus")
#sp.names
runs <- c(10)
# ENMs list:
SDMs <- c("GAM")
# List of eval_runs:
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5") ## !! need an underscore because biomod2 adds one by default in pred and calib_lines !! 
# List of prevalence levels:
prevalences <- c("n1.5")

library("parallel")
#mclapply(sp.names, SpeciesNicheModelling, mc.cores = 3)
SpeciesNicheModelling("Calanus_helgolandicus")

### Works fine. Results are stored as follows :
# - one main file per species, containing
# - one file containing all models' ouputs (with sub-files for each run of pseudo-absences)
# - one file per PROJECTIONS (i.e for each run of pseudoAbs, each time period and each scenario : 130 per species then) which contain R objects for :
#		o projection outputs (for each run of each model type) ; can be extracted with apply with mean() ; 
#		  probabilities are rescaled bewteen 0 and 1000 : divide by 1000 to rescale bewteen 0 and 1
#		o TSS bins
#		o projection ouputs for ensemble modelling : committee averaging (CA) & model averaging (mean)
#		o Clamping Mask	


##### 09/05/2016 : To monitor the evolution of the above script: 
setwd(second.wd)
species <- dir() ; species
species <- species[-c(43,59,83)]

### With a lapply: 
res <- lapply(species, function(sp) {
		setwd(paste(second.wd,"/",sp,"/", sep= ""))
		#dir()[grep(".txt", dir())] 
		n <- length(dir()[grep(".txt", dir())] )
		return(c(sp, n))
}) # eo lapply

prog <- data.frame(do.call(rbind, res)) ; rm(res)
colnames(prog)[1:2] <- c("species", "nb")
prog

### Still that freaking problem with semicolons :()
# To check number of species that have finished:
nrow(prog[which(prog$nb == "50"),])
nrow(prog[which(prog$nb == "50"),]) / 106


### Crap...it's a mix of runs, prevalences, species...need to re-run these the models below: 

# Calanus.helgolandicus       	10       1.5 GAM
# Calanus.helgolandicus       	10       1.5 GAM
# Calanus.helgolandicus       	10       1.5 GAM
# Calanus.helgolandicus       	10       1.5 GAM

