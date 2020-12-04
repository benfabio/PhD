
##### 01/06/2015 - LOV - Fabio Benedetti
##### Script for : 
#					- creating a function that computes maximal TSS thresholds for all binomials runs (3000 per species...cf.Script 20.3)
#					- compute a similar function that will compute TSS maximal threshold from HSI values (should be the same as the ones given by 'biomod2') 
#					- compare and act consequently

### Last update : 12/10/2015

# ------------------------------------------------------------------------------------------------------------------------------------

library("biomod2")
library("stringr")
library("foreign")
library("SDMTools")
library("gtools")
library("reshape2")
library("sp")
library("parallel")
library("ggplot2")
library("RColorBrewer")
library("scales")
library("fields")

# ------------------------------------------------------------------------------------------------------------------------------------

WD <- getwd()

# For Global coordinates :
coords_glob <- as.data.frame(ras <- brick("clims_woa13_baseline_04.grd"), xy = TRUE)
coords_glob <- coords_glob[,c(1:2,5)] # !! Loaded coordinates twice !!
coords_glob <- na.omit(coords_glob)
coords_glob$ids <- paste(coords_glob$x, coords_glob$y, sep="_")

# Basic stuff to start with
setwd("./niche_modelling")
second.wd <- getwd()
setwd(second.wd)
sp.names <- dir() ; length(sp.names)
sp.names <- sp.names[-c(6:10,48,56,65,89,116)]
setwd(WD) ; length(sp.names)

# List of pseudo-absences runs : n = 10
runs <- as.character(c(1:10))
# List of evaluation runs : n = 3
RUNS <- c("RUN1","RUN2","RUN3")
# List of SDMs
SDMs <- c("SRE","CTA","RF","MARS","FDA","GLM","GAM","ANN","GBM","MAXENT")
# List of binomial experiments
experiments <- c("exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10")
gc()

### BONUS to create a data table containing all possible combinations of the above lists : expand.grid()
# D <- expand.grid(pabs_run = runs, eval_run = RUNS, SDM = SDMs, binom_exp = experiments, scores = c("TSS","AUC"), value = NA)
# dim(D)

# ------------------------------------------------------------------------------------------------------------------------------------

##### Stratégie à suivre : créer un Rdata (de classe 'multidimensional array') pour CHAQUE espèce contenant les dimensions suivantes : 
# - 10 pseudo-absence runs
# - 10 SDM algorithms
# - 3 evaluation runs 
# - 10 binomial experiments
		

### For testing :
sp <- "Rhincalanus.nasutus"
r <- 7
R <- "RUN1"
sdm <- "GLM"
exp <- "exp1"

### Create a functioj to compute TSS and AUC scores from a ddf thatyou will build within the function below
confusion.matrix.operator <- function(x) {
			# Build a confusion matrix from observed testing set (presences and pseudo-absences) 
			# All predicted values > 0.5 will be turned into 1, which does not matter because we already transformed them into 1 and zeros
			cmx <- SDMTools::confusion.matrix(obs = x$observed, x$predicted, threshold = 0.5)
			# Compute AUC :
			AUC <- SDMTools::auc(obs = x$observed, pred = x$predicted)
			# Compute TSS :
			sensit <- SDMTools::sensitivity(mat = cmx)
			specif <- SDMTools::specificity(mat = cmx)
			TSS <- (sensit + specif - 1)
			# Remove some stuff and return a table containing both evaluation scores
			rm(sensit, specif, cmx)
			return( data.frame(AUC = AUC, TSS = TSS) )
} # eo FUN


### Create function : 
binom_eval_scorer <- function(sp = sp.names) {

				# Message
				message(paste("Running evaluations for ", sp, sep=""))
				setwd(paste(second.wd,"/","Binom_eval_scores","/", sep=""))
				dir.create(sp)
				setwd(second.wd)

				### 1°) Create a multidimensional array where you will store the calculated TSS and AUC scores
				eval_marray <- array(data = NA, 
								dim = c(2, 10, 3, 10, 10), 
						        dimnames= list(
								c("TSS","AUC"),
								runs,
								RUNS, 
						        SDMs,
						        experiments))
										
				### 2°) For each pseudo-absence runs, for each evaluation runs, for each SDM algorithm and for each binomial experiment
				###     fill the multidimensional array
				
				for(r in runs) {
					
					for(R in RUNS) {
						
							### Load species testing set : 
							setwd( paste(second.wd, "/", sp, "/","eval_runs_data","/", sep="") )  # Watchout here, need to replace points by '_' in the species name : gsub(".","_",sp, fixed=TRUE)
							Data <- load(paste("Data_eval_run",r,"_",R,".Rdata", sep="")) 
							DataSpecies <- get(Data) ; rm(Data)
							# summary(DataSpecies)
							# dim(DataSpecies)
							# Restrict DataSpecies to testing set :
							Data <- DataSpecies[DataSpecies$is_calib == FALSE, c("x","y","response")]
							# nrow(Data[Data$response == 1,])
							# nrow(Data[Data$response == 0,])
							# 0.2*nrow(DataSpecies)
							# Gut gut

							# Give a cell id :
							Data$ids <- as.factor(paste(Data$x, Data$y, sep="_"))
							rm(DataSpecies)  
					
							### Load the marray of binomial ouputs :
							setwd( paste(second.wd, "/","Binomials_T0_glob","/",sp,"/", sep="") ) 
							res_bin <- load( paste(sp,"_run",r,"_T0_glob.Rdata", sep="") )
							marray <- get(res_bin)  
							rm(res_bin)
							# str(marray)
							# head(marray)
							gc()

							# Marray contains the binary outputs for all SDMs and all experiments, that is why these loops are started afterwards :
							for(sdm in SDMs) {
										
								 # Keep the outputs from the chosen SDM algorithm
								 marray2 <- marray[,1,1,sdm,]
									
										for(exp in experiments) {
									
												### Create a ddf containing the binary ouputs :
												ddf2 <- data.frame(x= coords_glob[,"x"], y= coords_glob[,"y"], pred = marray2[,exp], ids = coords_glob[,"ids"])
												# restric to matching cells : 
												matching <- ddf2[,c("ids")] %in% Data[,c("ids")]
												ddf2 <- ddf2[matching,]
												
												### ¡! Data$response and ddf2$pred have the same length, but do not share the same order ¡!
												### Need to re-order the latter according to the former : 
												# length(unique(ddf2$ids))
												# length(unique(Data$ids))
												ddf2 <- ddf2[order(ddf2$ids),] 
												Data <- Data[order(Data$ids),]
												### dim(Data[which(Data$response == 0),])
												# head(ddf2) ; head(Data)
												# tail(ddf2) ; tail(Data)

												### Bind response var and predicted var : 
												data_evaluation <- data.frame(ID = Data$ids, observed = Data$response, predicted = ddf2$pred)
												rm(ddf2, matching)

												### With package 'SDMTools', compute Area Under the Curve (AUC) and True Skill Statistics (TSS) scores :
												library("SDMTools")
												table <- confusion.matrix.operator(data_evaluation)
									
												### Supply scores' value to marray
												eval_marray["TSS", r, R, sdm, exp] <- table$TSS
												eval_marray["AUC", r, R, sdm, exp] <- table$AUC
												rm(AUC, TSS, data_evaluation)
												gc()
													
										 } # eo exp in experiments	
										
								 rm(marray2)
						
							} # eo sdm in SDMs 
								
						rm(marray)	
						
						### Go to scores directory, and print 'eval_marray' as Rdata :
						setwd(paste(second.wd,"/","Binom_eval_scores","/",sp,"/", sep=""))
						save(eval_marray, file = paste(sp,"_run",r,"_",R,".Rdata", sep=""))
						rm(eval_marray)
						setwd(WD)	
						
				} # eo R in RUNS				
							
		} # eo r in runs
				
} # eo FUN


### Apply function within a mclapply : 
mclapply(sp.names, binom_eval_scorer, mc.cores = 12)
### ¡! WAY TOO SLOW ¡!
### ¡! Scores way lower than thresholded ones ¡!

### With paralleled l_ply
#require("plyr")
#require("doParallel")
#registerDoParallel(cores = 12)
#l_ply(sp.names, binom_eval_scorer, .parallel = TRUE)


# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------



#### For testing :
# sp <- sp.names[6]
# r <- 1
# #sdm <- "GLM"
## R <- 1


###### Create a similar function as above, but only for HSI values : 
#binom_eval_scorer2 <- function(sp = sp.names) {

#				# Message
#				message(paste("Running evaluations for ",sp, sep=""))

#				### 1°) Create a multidimensional array where you will store the calculated TSS and AUC scores
#				eval_marray <- array(data = NA, 
#								dim = c(2, 10, 10, 3), 
#						        dimnames= list(
#								c("AUC","TSS"),
#								runs,
#								SDMs, 
#						        RUNS) )
#										
#				### 2°) For each pseudo-absence runs, for each evaluation runs, for each SDM algorithm and for each binomial experiment
#				###     Fill the multidimentsional array
#				
#				for(r in runs) {
#					
#					for(sdm in SDMs) {
#						
#						for(R in RUNS) {
#								
#									### Load species testing set : 
#									setwd( paste(second.wd, "/", sp, "/","eval_runs_data","/", sep="") )  # Watchout here, need to replace points by '_' in the species name : gsub(".","_",sp, fixed=TRUE)
#									Data <- load(paste("Data_eval_run",r,"_",R,".Rdata", sep="")) 
#									DataSpecies <- get(Data) ; rm(Data)
#									# summary(DataSpecies)
#									# dim(DataSpecies)

#									# Restrict DataSpecies to testing set :
#									Data <- DataSpecies[DataSpecies$is_calib == FALSE, c("x","y","response")]

#									# Give a cell id :
#									Data$ids <- paste(Data$x, Data$y, sep="_") 
#								
#									### Load HSI values :
#									if(sdm == "SRE") { 
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_sre")])/1000)
#									} else if (sdm == "CTA") { 
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_cta")])/1000)
#									} else if (sdm == "RF") {
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_rf")])/1000)
#									} else if (sdm == "GBM") {
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_gbm")])/1000)
#									} else if (sdm == "GLM") {
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_glm")])/1000)
#									} else if (sdm == "GAM") {
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_gam")])/1000)
#									} else if (sdm == "MARS") {
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_mars")])/1000)
#									} else if (sdm == "FDA") {
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_fda")])/1000)
#									} else if (sdm == "ANN") {
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_ann")])/1000)
#									} else if (sdm == "MAXENT") {
#											predi <- ((DataSpecies[DataSpecies$is_calib == FALSE, c("pred_maxent")])/1000)
#									} # eo if else # 2

#									### Create a ddf containing the HSI values + the coordinates :
#									ddf2 <- data.frame(	x= DataSpecies[DataSpecies$is_calib == FALSE, c("x")], 
#														y = DataSpecies[DataSpecies$is_calib == FALSE, c("y")], 
#													   	pred = predi)
#									
#									ddf2$ids <- paste(ddf2$x, ddf2$y, sep="_")
#									
#									# restric to matching cells : 
#									matching <- ddf2[,c("ids")] %in% Data[,c("ids")]
#									#sum(matching)
#									ddf2 <- ddf2[matching,]
#									# dim(ddf2)
#									# dim(Data)

#									### Bind response var and predicted var : 
#									data_evaluation <- data.frame(ID = Data$ids, observed = Data$response, predicted = ddf2$pred)
#									rm(ddf2, DataSpecies)

#									### With package 'SDMTools', compute maximal Area Under the Curve (AUC) and True Skill Statistics (TSS) scores 
#									### according to a varying hsi threshold value
#									library("SDMTools")
#									# Compute AUC first:
#									AUC <- SDMTools::auc(obs = data_evaluation$observed, pred = data_evaluation$predicted)
#									
#									# Build a ddf containing the threshold values a,d their corresponding scores : 
#									results <- sapply(seq(from= 0, to= 1, by= 0.01), function(t) {
#												cmx <- SDMTools::confusion.matrix(obs = data_evaluation$observed, data_evaluation$predicted, threshold = t)
#												sensit <- SDMTools::sensitivity(mat = cmx)
#									  			specif <- SDMTools::specificity(mat = cmx)
#									  			tss <- (sensit + specif - 1)  # Calcule du TSS
#									  			res <- c(t, tss, sensit, specif) # On stocke le threshold, TSS, sensitivity et specificity
#									}) # eo sapply
#									
#									results <- t(results) # On transpose pour avoir les thresholds et TSS en colonnes
#									results <- data.frame(results)
#									colnames(results) <- c("threshold","tss", "sensitivity", "specificity")

#									TSSmax <- max(results$tss) # Get maximum TSS value (should be the same as the one computed by 'biomod2')
#									
#									### Supply scores' value to marray
#									eval_marray["TSS",r,sdm,R] <- TSSmax
#									eval_marray["AUC",r,sdm,R] <- AUC
#									
#									rm(AUC, TSSmax, data_evaluation, results)
#						
#						} # eo R in RUNS	
#						
#					} # eo sdm in SDMs 
#					
#				} # eo r in runs
#				
#				### Go to scores directory, and print 'eval_marray' as Rdata :
#				setwd(paste(second.wd,"/","HSI_eval_scores","/", sep=""))
#				save(eval_marray, file = paste(sp,".Rdata", sep=""))
#				rm(eval_marray)
#				
#		setwd(WD)			

#} # eo FUN


#### Apply function within a mclapply : 
#mclapply(sp.names, binom_eval_scorer2, mc.cores = 12)
#### Works fine, and quickly.



###### Check : 
setwd(paste(second.wd,"/","Eval_scores_28_05_15","/", sep=""))
dir()
scores <- read.table("Eval_scores_Rhincalanus.nasutus.txt", h=TRUE)  
str(scores)
median(scores[which(scores$stat.type == "TSS"),"Testing.data"])
## melt multidimensional array : 
m_scores <- melt(scores)
colnames(m_scores) <- c("stat","pabs_run","sdm","eval_run","value")
rm(scores)
m_scores$sp.name <- "Calanus.helgolandicus"

#### Ok, within a lapply, rbind all species'data: 
results <- lapply(sp.names, function(sp) {	
		
				# Read each file
				setwd(paste(second.wd,"/","HSI_eval_scores","/", sep=""))
				res <- load(paste(sp,".Rdata",sep=""))
				scores <- get(res) ; rm(res)	
				# Melt
				m_scores <- melt(scores) ; rm(scores)			
				# Change colnames
				colnames(m_scores) <- c("stat","pabs_run","sdm","eval_run","value")
				# Add species name
				m_scores$sp.name <- sp

				return(m_scores)
				setwd(second.wd)
}) #eo lapply

scores <- do.call(rbind, results) ; rm(results)
dim(scores)
summary(scores)
scores[which(scores$value <= 0),"value"] <- 0

#### Compute mean TSS and AUC per species : 
library("dplyr")
avg_tss <- scores[scores$stat == "TSS",] %>%
  group_by(sdm) %>%
  summarise(average = mean(value), sd = sd(value))

avg_tss <- data.frame(avg_tss[order(avg_tss$average, decreasing = T),]) # sort by order

require("dplyr")
avg_auc <- scores[scores$stat == "AUC",] %>%
  group_by(sdm) %>%
  summarise(average = mean(value), sd = sd(value))
  
avg_auc <- data.frame(avg_auc[order(avg_auc$average, decreasing = T),]) # sort by order

#### Very similar to the scores given by 'biomod2'

### Boxplots : 



##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################


### 10/06/2015 : Examine results from binom_eval_scores : 
library("reshape2")
library("dplyr")

#files <- dir()
#res <- load(files[3]) # retrieving 3rd file, because it is the last one that was created
#scores <- get(res) ; rm(res)

scores <- read.csv("eval_scores_thresh_29_09_15.csv", h=TRUE, sep=";", dec=",")
head(scores) ; tail(scores)
class(scores)
str(scores)

#scores <- melt(scores)
#dim(scores)
#colnames(scores)[1:5] <- c("stat","pabs_run","eval_run","sdm","binom_exp")
#head(scores)
#scores[is.na(scores$value),] # OK no NAs

avg_tss <- scores[scores$stat == "TSS",] %>%
  	group_by(model) %>%
  	summarise(average = mean(Testing.data), sd = sd(Testing.data))

avg_tss <- data.frame(avg_tss[order(avg_tss$average, decreasing = T),]) # sort by order
avg_tss
mean(avg_tss$average)

write.csv2(avg_tss,"SDM_avgTSS_thresh_29_09_15.csv", row.names=FALSE)



avg_auc <- scores[scores$stat == "ROC",] %>%
 	 group_by(model) %>%
 	 summarise(average = mean(Testing.data), sd = sd(Testing.data))
  
avg_auc <- data.frame(avg_auc[order(avg_auc$average, decreasing = T),]) # sort by order
avg_auc
mean(avg_auc$average)
write.csv2(avg_auc,"SDM_avgAUC_thresh_29_09_15.csv",row.names=FALSE)




### Now, for all species : 
# Ok, within a lapply, rbind all species'data:
second.wd <- getwd()
sp.names <- dir()[93]
runs <- as.character(c(2:10))
R <- "RUN3"

results <- lapply(sp.names, function(sp) {	

	for(r in runs) {
			
		#for(R in RUNS) {
			
				# Read each file
				setwd(paste(second.wd,"/",sp,"/", sep=""))
				res <- load(paste(sp,"_run",r,"_",R,".Rdata",sep=""))
				scores <- get(res) ; rm(res)	
				# Melt
				m_scores <- melt(scores) ; rm(scores)			
				# Change colnames
				colnames(m_scores) <- c("stat","pabs_run","sdm","eval_run","value")
				# Add species name
				m_scores$sp.name <- sp
				
				return(m_scores)
				setwd(second.wd)
				
		   #}
		}
		
}) #eo lapply

scores <- do.call(rbind, results) ; rm(results)
colnames(scores)[1:7] <- c("stat","pabs_run","eval_run","sdm","experiment","score","sp.name") 
## dim(scores)
## summary(scores)
scores[which(scores$score <= 0),"score"] <- 0

scores[is.na(scores),c("pabs_run")]
scores <- na.omit(scores)

scores[which(scores$stat == "AUC"),]

write.csv2(scores,"eval_scores_proba_29_09_15.csv", row.names = FALSE)

scores <- read.csv("evaluation_scores_29_09_15.csv", h=TRUE, sep=";", dec=",")
ggplot(scores[which(scores$stat == "AUC"),], aes(factor(sdm), value)) + geom_boxplot() + xlab("SDM") + ylab("AUC")


### 12/10/2015 : Retrieving sepcies' and SDMs' thresholded TSS/AUC scores
library("reshape2")
library("dplyr")

result <- lapply(sp.names, function(sp) {
		# Go to species dir
		setwd(paste(second.wd,"/","Binom_eval_scores","/",sp,"/", sep=""))
		# Load/get .Rdata
		res <- load(paste(sp,"_eval_marray.Rdata", sep=""))
		res <- get(res)
		# Melt multidim array
		molten <- melt(res)
		molten$sp.name <- sp
		# Return object
		return(molten)

	} # eo FUN

) #eo lapply

scores <- do.call(rbind, result)
rm(result)
colnames(scores)[1:5] <- c("stat","pabs_run","eval_run","sdm","exp")
summary(scores)
scores[which(scores$value <= 0),"value"] <- 0

### Compute species' average TSS/AUC scores 
avg_tss <- scores[scores$stat == "TSS",] %>%
  	group_by(sdm) %>%
  	summarise(average = mean(value), sd = sd(value))

avg_tss <- data.frame(avg_tss[order(avg_tss$average, decreasing = T),]) # sort by order
(nrow(avg_tss[which(avg_tss$average >= 0.9),])/106)*100

#avg_tss
#mean(avg_tss$average)
setwd(WD)


write.csv2(avg_tss,"SDM_AUC_probabilistic_12_10_15.csv", row.names=FALSE)

avg_auc <- scores[scores$stat == "ROC",] %>%
 	 group_by(model) %>%
 	 summarise(average = mean(Testing.data), sd = sd(Testing.data))
  
avg_auc <- data.frame(avg_auc[order(avg_auc$average, decreasing = T),]) # sort by order
avg_auc
mean(avg_auc$average)
write.csv2(avg_auc,"SDM_avgAUC_thresh_29_09_15.csv",row.names=FALSE)



### 07/10/2015 : Re-computed R.nasutus evaluation scores (AUC/TSS):
scores <- load("Rhincalanus.nasutus_eval_marray.Rdata")
scores <- get(scores)
dim(scores)
str(scores)
summary(scores)
m_scores <- melt(scores) ; rm(scores)	
colnames(m_scores)[1:5] <- c("stat","pabs_run","eval_run","sdm","exp")
summary(m_scores) # No NAs gut gut :)

# Average and sd of species' TSS score
mean(m_scores[which(m_scores$stat == "TSS"),"value"])
sd(m_scores[which(m_scores$stat == "TSS"),"value"])

mean(m_scores[which(m_scores$stat == "AUC"),"value"])
sd(m_scores[which(m_scores$stat == "AUC"),"value"])



### 12/10/2015 : Plotting plotting
scores <- read.table("scores.csv", h=T, dec=".", sep=";")
dim(scores)
head(scores) ; tail(scores)
summary(scores)

### Relation between average TSS and species prevalence  
p <- ggplot(scores) + geom_point(aes(x=log(n), y=avg_TSS_t, size = sd_TSS_t), shape = 21, fill = "blue", colour = "black") +
				scale_size_continuous("Standard deviation\nof average TSS") + scale_y_continuous(limits=c(0.5,1)) +
				xlab("Number of presences (Log)") + ylab("Average TSS") + theme_bw()
				
pdf("plot1.pdf")
p
dev.off()
	
				
### Relation between average AUC and species prevalence  
p <- ggplot(scores) + geom_point(aes(x=log(n), y=avg_AUC_t, size = sd_AUC_t), shape = 21, fill = "blue", colour = "black") +
				scale_size_continuous("Standard deviation\nof average AUC") + scale_y_continuous(limits=c(0.8,1)) +
				xlab("Number of presences (Log)") + ylab("Average AUC") + theme_bw()
				
pdf("plot2.pdf")
p
dev.off()	


### Relation between average TSS (probabilistic) and species prevalence  
p <- ggplot(scores) + geom_point(aes(x=log(n), y=avg_TSS_p, size = sd_TSS_p), shape = 21, fill = "red", colour = "black") +
				scale_size_continuous("Standard deviation\nof average TSS") + scale_y_continuous(limits=c(0.3,1)) +
				xlab("Number of presences (Log)") + ylab("Average TSS") + theme_bw()
				
pdf("plot3.pdf")
p
dev.off()		
	
				
### Relation between average AUC (probabilistic) and species prevalence  
p <- ggplot(scores) + geom_point(aes(x=log(n), y=avg_AUC_p, size = sd_AUC_p), shape = 21, fill = "red", colour = "black") +
				scale_size_continuous("Standard deviation\nof average AUC") + scale_y_continuous(limits=c(0.7,1)) +
				xlab("Number of presences (Log)") + ylab("Average AUC") + theme_bw()
				
pdf("plot4.pdf")
p
dev.off()	


### Relation between probabilistic TSS and thresholded TSS
p <- ggplot(scores) + geom_point(aes(x=avg_TSS_t, y=avg_TSS_p), shape = 21, fill = "red", colour = "black") +
				scale_y_continuous(limits=c(0.4,1)) + scale_x_continuous(limits=c(0.6,1)) +
				xlab("TSS (threshold)") + ylab("TSS (probabilistic)") + theme_bw()
				
pdf("plot5.pdf")
p
dev.off()	


### Relation between probabilistic AUC and thresholded AUC
p <- ggplot(scores) + geom_point(aes(x=avg_AUC_t, y=avg_AUC_p), shape = 21, fill = "red", colour = "black") +
				scale_y_continuous(limits=c(0.7,1)) + scale_x_continuous(limits=c(0.8,1)) +
				xlab("AUC (threshold)") + ylab("AUC (probabilistic)") + theme_bw()
				
pdf("plot6.pdf")
p
dev.off()



### Relation between species prevalence and TSS range
p <- ggplot(scores) + geom_point(aes(x=log(n), y=sd_AUC_p), shape = 21, fill = "orange", colour = "black", size = 3) +
				scale_y_continuous(limits=c(0, 0.3)) +
				xlab("Number of presences (Log)") + ylab("Standard deviation of AUC (probabilistic)") + theme_bw()
				
pdf("plot9.pdf")
p
dev.off()

