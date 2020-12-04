
##### 25/06/2015 - LOV - Fabio Benedetti 
##### Script for : 
#		- loading the ANOVA table (big)
#		- testing three-way variance analysis (ANOVA) with replication (because of 100 communities)
#		- divide the big ANOVA table into 10 fragments of 2649 cells
#		- in a for loop, extract each of ten fragments, and run aov with ddply, by paralleling on ten cores
#		- test potential impact of community as a source of uncertainty

### Last update : 15/12/2015

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
## For Mediterranean coordinates :
varT0 <- read.table("climatos_woa13_6594_Levin.txt", header=T, sep="\t", dec=".")
# Vector containing cells'id : 
ids <- paste(varT0[,"x"], y = varT0[,"y"], sep="_")
# length(ids)
# List of pseudo-absences runs
runs <- as.character(c(1:10)) ; runs
# List of SDMs
SDMs <- c("SRE","CTA","RF","MARS","FDA","GLM","GAM","ANN","GBM","MAXENT")
# List of binomial experiments 
experiments <- c("exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10")
# try sample() function 
# exp <- sample(experiments, 1) ; exp
# List of scenarios 
scenarios <- c("A2","A2F","A2RF","A2ARF")
# List of potential communities :
communities <- c(1:100)
# List of indices : 
indices <- c("SRT0","SRT2","∆SR","Jac","Jtu","Jne","Jratio")



### For Model_forcing :
for(i in 1:10) {
	
		# Useless message :
		message(paste("Doing ANOVA for subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables","/","Model_forcing","/", sep=""))
		table <- load(paste("ANOVA_subtable_T1_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# Restrict to columns of interest
		table <- table[,-c(6,7,9,10,11)]
		table$run <- factor(table$run)
		table$community <- factor(table$community)
		
		### Get rid of SRE projections !
		table <- table[which(table$sdm != "SRE"),]
		# dim(table)
		
		require("lme4")
		require("doParallel")
		require("plyr")
		registerDoParallel(cores = 12)
		
		sum_squares <- daply(table, ~ id, function(X) {
					# Do ANOVA (without interactions)
					m <- lmer(formula =  dSR ~ (sdm + scenar)^2 + (1 | community/run ), data= X[, c(2:6)] )
					# Retrieve SSQ
					ssq <- anova(m)[,2]
					names(ssq) <- rownames(anova(m))
					ssq <- t(as.data.frame(ssq))
					
					return(ssq)
					
					rm(m)

		}, .parallel=TRUE)
		gc()
		
		sum_squares <- as.data.frame(sum_squares)
		sum_squares$pixel <- factor(rownames(sum_squares))
		sum_squares <- na.omit(sum_squares)
		
		
		# Save in dir:
		message(paste("Printing F-ratios for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables","/","Model_forcing","/","SSQ_15_01_16","/", sep="")) 
		save(sum_squares, file = paste("ssq_dSR_subtable_T1_",i,".Rdata", sep=""))
		
		# Make room:
		rm(table, sum_squares)
		gc()
		setwd(paste(WD,"/","ANOVA_tables","/", sep=""))
		
} # eo for loop





### For Emission_scenario :
for(i in 1:10) {
	
		# Useless message :
		message(paste("Doing ANOVA for subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables","/","Emission_scenarios","/", sep=""))
		table <- load(paste("ANOVA_subtable2_T1_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# restrict to columns of interest
		table <- table[,-c(6,7,9,10,11)]
		table$run <- factor(table$run)
		table$community <- factor(table$community)
		
		### Get rid of SRE projections !
		table <- table[which(table$sdm != "SRE"),]
		# dim(table)
		
		require("lme4")
		require("doParallel")
		require("plyr")
		registerDoParallel(cores = 12)
		
		sum_squares <- daply(table, ~ id, function(X) {
					# Do ANOVA (without interactions)
					m <- lmer(formula =  dSR ~ (sdm + scenar)^2 + (1 | community/run ), data= X[, c(2:6)] )
					# Retrieve SSQ
					ssq <- anova(m)[,2]
					names(ssq) <- rownames(anova(m))
					ssq <- t(as.data.frame(ssq))
					
					return(ssq)
					
					rm(m)

		}, .parallel = TRUE)
		gc()
		
		sum_squares <- as.data.frame(sum_squares)
		sum_squares$pixel <- factor(rownames(sum_squares))
		sum_squares <- na.omit(sum_squares)
		
		# Save in dir:
		message(paste("Printing Sum of Squares for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables","/","Emission_scenarios","/","SSQ_15_01_16","/", sep="")) 
		save(sum_squares, file = paste("ssq_dSR_subtable2_T1_",i,".Rdata", sep=""))
		
		# Make room:
		rm(table, sum_squares)
		gc()
		setwd(paste(WD,"/","ANOVA_tables","/", sep=""))
		
} # eo for loop		

