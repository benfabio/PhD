
# ------------------------------------------------------------------------------------------------------------------------------------

library("stringr")
library("reshape2")
library("sp")
library("matrixStats")
library("ggplot2")
library("RColorBrewer")
library("scales")
library("fields")
library("parallel")
library("doParallel")

# ------------------------------------------------------------------------------------------------------------------------------------

### Basic stuff to start with
WD <- getwd()
# Vector containing cells'id : 
# List of pseudo-absences runs
runs <- c(1:10) ; runs
# List of SDMs
SDMs <- c("SRE","CTA","RF","MARS","FDA","GLM","GAM","ANN","GBM","MAXENT")
# List of scenarios 
scenarios <- c("A2","A2F","A2RF","A2ARF")
# List of prevalence levels :
prevalences <- c("n1","n1.5","n2","n10","n50")
# List of indices : 
indices <- c("SRT0","SRT2","∆SR","Jac","Jtu","Jne","Jratio")
# evaluation runs
eval_runs <- c("RUN1","RUN2","RUN3")

# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------

### 1) For ∆SR
# i <- 1

for(i in 1:10) {
	
		# Useless message :
		message(paste("Loading subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		table <- load(paste("ANOVA_subtable_all_T1_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# Register number of cores for parallel computing
		library("doParallel")
		registerDoParallel(cores = 20) ### BEWARE OF THE MACHINE YOU'RE USING ###
		# Compute means in parallel
		message(paste("Computing average predictions for subtable n°",i, sep=""))
		library("plyr")
		# Restrict 'table' object to columns of interest (make room)
		# table <- table[,-c(6,7,8,9,10)]
		# summary(table)
		
		means <- ddply(table, ~ id, function(X) {
					# Compute mean dSR/Jac for each SDM:
					avg <- mean( X[,"dSR"] )
					sd <- sd(X[,"dSR"])
					#coeff <- (sd/avg)*100  # coefficient of variation
					# Return a ddf containing all these values : 
					means <- data.frame(avg = avg, sd = sd)
					
		}, .parallel = TRUE)
		gc()
		
		### NOTE : 'means' object kept cells' id as a column
		# summary(means)
		# Save in dir:
		message(paste("Printing average ∆SR values for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		save(means, file = paste("avg_all_dSR_T1_subtable_",i,".Rdata", sep=""))
		
		# Make room:
		rm(means, table)
		gc()
		setwd(paste(WD,"/","ANOVA_tables","/", sep=""))
		
} # eo for loop	



### 2) For Jaccard index
#i <- 5

for(i in 1:10) {
	
		# Useless message :
		message(paste("Loading subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		table <- load(paste("ANOVA_subtable_all_T1_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		# dim(table)
		
		# Register number of cores for parallel computing
		library("doParallel")
		registerDoParallel(cores = 20) ### BEWARE OF THE MACHINE YOU'RE USING ###
		# Compute means in parallel
		message(paste("Computing average predictions for subtable n°",i, sep=""))
		library("plyr")
		# Restrict 'table' object to columns of interest (make room)
		# table <- table[,-c(6,7,8,9,10)]
		# summary(table)
		# summary(table[is.na(table$Jac),])
		### NEED TO FIX SRE Jaccard values to 0
		table[is.na(table$Jac ),"Jac"] <- 0
		### 08/06/16: What do we do with the 999 values ? (i.e. cases where non analog conditions to dont allow to compute Jaccard index)
		table[which(table$Jac == 999),"Jac"] <- NA
		
		means <- ddply(table, ~ id, function(X) {
					# Compute mean dSR/Jac for each SDM:
					avg <- mean( X[,"Jac"] )
					sd <- sd(X[,"Jac"])
					#coeff <- (sd/avg)*100  # coefficient of variation
					# Return a ddf containing all these values : 
					means <- data.frame(avg = avg, sd = sd)					
		}, .parallel = TRUE)
		gc()
		
		### NOTE : 'means' object kept cells' id as a column
		# summary(means)
		# Save in dir:
		message(paste("Printing SDMs average values for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		save(means, file = paste("avg_all_Jac_T1_subtable_",i,".Rdata", sep=""))
		
		# Make room:
		rm(means, table)
		gc()
		setwd(paste(WD,"/", sep=""))
		
} # eo for loop	



### 3) For Jaccard's component of Nestedness (beware of Bioclim NULL values)
#i <- 1

for(i in 1:10) {
	
		# Useless message :
		message(paste("Loading subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		table <- load(paste("ANOVA_subtable_all_T1_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# Register number of cores for parallel computing
		library("doParallel")
		registerDoParallel(cores = 20) ### BEWARE OF THE MACHINE YOU'RE USING ###
		# Compute means in parallel
		message(paste("Computing average predictions for subtable n°",i, sep=""))
		library("plyr")
		# Restrict 'table' object to columns of interest (make room)
		# table <- table[,-c(6,7,8,9,10)]
		### Fix SRE's Jne values to 1 !
		table[is.na(table$Jne ),"Jne"] <- 1
		table[is.na(table$Jtu ),"Jtu"] <- 0
		
		### 08/06/16: What do we do with the 999 values ? (i.e. cases where non analog conditions to dont allow to compute Jaccard index)
		table[which(table$Jac == 999),"Jne"] <- NA
		
		means <- ddply(table, ~ id, function(X) {
					# Compute mean dSR/Jac for each SDM:
					avg <- mean( X[,"Jne"] )
					sd <- sd(X[,"Jne"])
					#coeff <- (sd/avg)*100  # coefficient of variation
					# Return a ddf containing all these values : 
					means <- data.frame(avg = avg, sd = sd)
					
		}, .parallel = TRUE)
		gc()
		
		### NOTE : 'means' object kept cells' id as a column
		# summary(means)
		# Save in dir:
		message(paste("Printing SDMs average values for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		save(means, file = paste("avg_all_Jne_T1_subtable_",i,".Rdata", sep=""))
		
		# Make room:
		rm(means, table)
		gc()
		setwd(paste(WD,"/", sep=""))
		
} # eo for loop	





### 4) For Jaccard's component of Turn-over (beware of Bioclim NULL values)
i <- 1

for(i in 1:10) {
	
		# Useless message :
		message(paste("Loading subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		table <- load(paste("ANOVA_subtable_all_T1_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# Register number of cores for parallel computing
		library("doParallel")
		registerDoParallel(cores = 20) ### BEWARE OF THE MACHINE YOU'RE USING ###
		# Compute means in parallel
		message(paste("Computing average predictions for subtable n°",i, sep=""))
		library("plyr")
		# Restrict 'table' object to columns of interest (make room)
		# table <- table[,-c(6,7,8,9,10)]
		table[is.na(table$Jne ),"Jne"] <- 1
		table[is.na(table$Jtu ),"Jtu"] <- 0
		
		### 08/06/16: What do we do with the 999 values ? (i.e. cases where non analog conditions to dont allow to compute Jaccard index)
		table[which(table$Jac == 999),"Jtu"] <- NA
		
		means <- ddply(table, ~ id, function(X) {
					# Compute mean dSR/Jac for each SDM:
					avg <- mean( X[,"Jtu"] )
					sd <- sd(X[,"Jtu"])
					#coeff <- (sd/avg)*100  # coefficient of variation
					# Return a ddf containing all these values : 
					means <- data.frame(avg = avg, sd = sd)					
		}, .parallel = TRUE)
		gc()
		
		### NOTE : 'means' object kept cells' id as a column
		# summary(means)
		# Save in dir:
		message(paste("Printing SDMs average values for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		save(means, file = paste("avg_all_Jtu_T1_subtable_",i,".Rdata", sep=""))
		
		# Make room:
		rm(means, table)
		gc()
		setwd(paste(WD,"/", sep=""))
		
} # eo for loop	




### 5) For ßratio (Jne/Jac)
i <- 1

for(i in 1:10) {
	
		# Useless message :
		message(paste("Loading subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		table <- load(paste("ANOVA_subtable_all_T1_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# Register number of cores for parallel computing
		library("doParallel")
		registerDoParallel(cores = 20) ### BEWARE OF THE MACHINE YOU'RE USING ###
		# Compute means in parallel
		message(paste("Computing average predictions for subtable n°",i, sep=""))
		library("plyr")
		# Restrict 'table' object to columns of interest (make room)
		# table <- table[,-c(6,7,8,9,10)]
		table[is.na(table$Jne),"Jne"] <- 1
		table[is.na(table$Jtu),"Jtu"] <- 0
		### With some SRE projections: total loss of species --> generates NA in ßratio because 100% loss of species, fix it 
		table[which(is.na(table$Jratio) & table$sdm != "SRE"),"Jratio"] <- 0
		### With some SRE projections: total loss of species --> generates NA in ßratio because 0% loss of species, fix it
		table[which(is.na(table$Jratio) & table$sdm == "SRE"),"Jratio"] <- 1
		
		### 08/06/16: What do we do with the 999 values ? (i.e. cases where non analog conditions to dont allow to compute Jaccard index)
		table[which(table$Jac == 999),"Jratio"] <- NA
		
		means <- ddply(table, ~ id, function(X) {
					# Compute mean dSR/Jac for each SDM:
					avg <- mean( X[,"Jratio"] )
					sd <- sd(X[,"Jratio"])
					#coeff <- (sd/avg)*100  # coefficient of variation
					# Return a ddf containing all these values : 
					means <- data.frame(avg = avg, sd = sd)					
		}, .parallel = TRUE)
		gc()
		
		### NOTE : 'means' object kept cells' id as a column
		# summary(means)
		# Save in dir:
		message(paste("Printing SDMs average values for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables_wout_non_analogs","/", sep=""))
		save(means, file = paste("avg_all_Jratio_T1_subtable_",i,".Rdata", sep=""))
		
		# Make room:
		rm(means, table)
		gc()
		setwd(paste(WD,"/", sep=""))
		
} # eo for loop	


####################################################################################################################################

