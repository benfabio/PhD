
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

# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------


### Load the big ANOVA table : 
setwd(paste(WD,"/","ANOVA_tables","/", sep=""))
table <- load("ANOVA_table.Rdata") ### Takes a while...
table <- get(table)
# dim(table) # 1 059 600 000 rows for 11 columns
# colnames(table)

### 06/07/15 : Dividing ANOVA table into 10 fragments of 2649 cells each. Print them as .Rdata in dir
table_1 <- subset(table, ids %in% ids[1:2649])
# dim(table_1)
# length(unique(table_1$id))
save(table_1, file = paste("ANOVA_subtable_1.Rdata", sep=""))

table_2 <- subset(table, ids %in% ids[2650:5298])
# dim(table_2)
# length(unique(table_2$id))
save(table_2, file = paste("ANOVA_subtable_2.Rdata", sep=""))

table_3 <- subset(table, ids %in% ids[5299:7947])
# dim(table_3)
# length(unique(table_3$id))
save(table_3, file = paste("ANOVA_subtable_3.Rdata", sep=""))

table_4 <- subset(table, ids %in% ids[7948:10596])
# dim(table_4)
# length(unique(table_4$id))
save(table_4, file = paste("ANOVA_subtable_4.Rdata", sep=""))

table_5 <- subset(table, ids %in% ids[10597:13245])
# dim(table_5)
# length(unique(table_5$id))
save(table_5, file = paste("ANOVA_subtable_5.Rdata", sep=""))

table_6 <- subset(table, ids %in% ids[13246:15894])
# dim(table_6)
# length(unique(table_6$id))
save(table_6, file = paste("ANOVA_subtable_6.Rdata", sep=""))

table_7 <- subset(table, ids %in% ids[15895:18543])
# dim(table_7)
# length(unique(table_7$id))
save(table_7, file = paste("ANOVA_subtable_7.Rdata", sep=""))

table_8 <- subset(table, ids %in% ids[18544:21192])
# dim(table_8)
# length(unique(table_8$id))
save(table_8, file = paste("ANOVA_subtable_8.Rdata", sep=""))

table_9 <- subset(table, ids %in% ids[21193:23841])
# dim(table_9)
# length(unique(table_9$id))
save(table_9, file = paste("ANOVA_subtable_9.Rdata", sep=""))

table_10 <- subset(table, ids %in% ids[23842:26490]) 
# dim(table_10)
# length(unique(table_10$id))
save(table_10, file = paste("ANOVA_subtable_10.Rdata", sep=""))

# Check whether there are redundancies betwwen the 2 datasets: 
# which(unique(table_10$id) %in% unique(table_7$id)) # No redundancies

rm(table)

# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------

### Test aov() function
# ?aov

# Subsample a cell based on a id : 
#length(unique(table$id))
#ids <- unique(table$id)
#i <- ids[10000]

require("dplyr")
t <- subset(table, id == unique(table$id)[1] )
dim(t)
head(t)

# Test aov
aov <- aov(formula = dSR ~ run + sdm + scenar + run*sdm+ run*scenar + sdm*scenar, data = t[,c(2,3,4,5,8)])
summary(aov)
# Extract sum of squares for each factor + interaction
summary(aov)[[1]][1:6,2]
# Gut gut

# When turning all explanatory variables into factors : 
#t$run <- as.factor(t$run)
#t$community <- as.factor(t$community)

aov <- aov(formula = dSR ~ run + sdm + scenar + run*sdm+ run*scenar + sdm*scenar, data = t[2,c(2,3,4,5,8)])
summary(aov)
# Extract sum of squares for each factor + interaction
summary(aov)[[1]][1:6,4]
summary(aov)[[1]]$F # To extract Fratio
summary(aov)[[1]]$S 


### 06/07/15 : Test code from ©J-O. Irisson
# fake data
pixel <- 1:2000
run <- letters[1:10]
sdm <- letters[1:10]
scenar <- letters[1:4]
commu <- 1:100
d <- expand.grid(pixel=pixel, run=run, sdm=sdm, scenar=scenar, commu=commu)
d$x <- runif(nrow(d))

# parallel computing
library("doParallel")
registerDoParallel(cores=10)

# compute ANOVA in parallel
library("plyr")
system.time(fratio <- daply(d, ~ pixel, function(X) {
	m <- aov(x ~ run + sdm + scenar, data=X)
 	fratio <- summary(m)[[1]]$F
}, .parallel=TRUE))

# handle the output
fratio <- as.data.frame(fratio[,-4])
names(fratio) <- c("run", "sdm", "scenar")
fratio$pixel <- as.numeric(row.names(fratio))

### Test it on one of your ANOVA_subtables you created above:
head(table_4)
# parallel computing
library("doParallel")
registerDoParallel(cores = 12)

# Compute ANOVA in parallel
library("plyr")
system.time(fratio <- daply(table_4, ~ id, function(X) {
	m <- aov(formula = dSR ~ run + sdm + scenar, data= X[,c(2,3,4,5,8)])
 	#fratio <- summary(m)[[1]]$F
	squares <- summary(m)[[1]]$F
}, .parallel=TRUE))

fratio <- as.data.frame(fratio[,-4])
names(fratio) <- c("run", "sdm", "scenar")
fratio$pixel <- as.factor(row.names(fratio))
head(fratio)
# OK...


# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------


### If everything seems ok, make a for loop (i in 1:10) that will read ANOVA_subtable and compute anova for each cell
### by paralleling with daply, like above, and print results in a Rdata. Final table should have 4 columns and 2649 rows


for(i in 1:10) {
	
		# Useless message :
		message(paste("Doing ANOVA for subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables","/","Model_forcing","/", sep=""))
		table <- load(paste("ANOVA_subtable_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# restrict to columns of interest
		table <- table[,-c(6,7,9,10,11)]
		table$run <- factor(table$run)
		table$community <- factor(table$community)
		
		# Register number of cores for parallel computing
		#library("doParallel")
		#registerDoParallel(cores = 12)
		# Compute ANOVA in parallel
		#library("plyr")
		#fratio <- daply(table, ~ id, function(X) {
					# Do ANOVA (without interactions)
					#m <- aov(formula =  dSR ~ run + sdm + scenar + sdm*scenar + sdm*run, data= X[,c(2:6)])
					# Retrieve F-ratio
		 			#fratio <- summary(m)[[1]]$F
		#}, .parallel=TRUE)
		#gc()
		
		############################################################## With LMER:
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
		
		#t <- table[which(table$id == "18.9940447244757_30.2989752360427"),c(2:6)]
		
		#m <- lmer(formula =  dSR ~ (sdm + scenar)^2 + (1 | community/run ), data= t)
		#fixef(m) # effets fixes
		#ranef(m) # effets aléatoires
		#summary(m)
		#ssq <- anova(m)[,2]
		#names(ssq) <- rownames(anova(m))
		
		#sos <- as.data.frame(ssq)
		sum_squares$pixel <- rownames(sum_squares)
		sum_squares <- na.omit(sum_squares)
		
		####################################################################################################################################
		 
		## Get rid of extra column
		#fratio <- as.data.frame(fratio[,-6])
		## Rename columns
		#names(fratio) <- c("run", "sdm", "scenar", "sdm:scenar", "sdm:run")
		# Add a column specifying cells' id
		#fratio$pixel <- as.factor(row.names(fratio))
		#dim(fratio)
		# Get rid of NAs
		#fratio <- na.omit(fratio)
		# dim(fratio)
		
		# Save in dir:
		message(paste("Printing F-ratios for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables","/","Model_forcing","/","SSQ_15/12","/", sep="")) 
		save(fratio, file = paste("Fratios_dSR_subtable2_",i,".Rdata", sep=""))
		
		# Make room:
		rm(table, fratio)
		gc()
		setwd(paste(WD,"/","ANOVA_tables","/", sep=""))

} # eo for loop


# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------

##### 15/12/15 : Re-computing sum of squares by using mixed models
#####			 the pseudo-absence realisation and the community randomly - probabilistically assembled cannot be investigated as 
#####			 uncertainty sources, like the ENM and the forcing configurations. We must use mixed models, with the "lme4" package,
#####			 to identify these two factors as "random effects", meaning they will be identified as pseudo-replicates.
#####			 Code remains the same, except that you will add a line in the daply


### For Model_forcing :
for(i in 1:10) {
	
		# Useless message :
		message(paste("Doing ANOVA for subtable n°",i, sep=""))
		setwd(paste(WD,"/","ANOVA_tables","/","Model_forcing","/", sep=""))
		table <- load(paste("ANOVA_subtable_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# Restrict to columns of interest
		table <- table[,-c(6,7,9,10,11)]
		table$run <- factor(table$run)
		table$community <- factor(table$community)
		
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
		setwd(paste(WD,"/","ANOVA_tables","/","Model_forcing","/","SSQ_15_12","/", sep="")) 
		save(sum_squares, file = paste("ssq_dSR_subtable_",i,".Rdata", sep=""))
		
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
		table <- load(paste("ANOVA_subtable2_",i,".Rdata",sep=""))
		table <- get(table)
		gc()
		
		# restrict to columns of interest
		table <- table[,-c(6,7,9,10,11)]
		table$run <- factor(table$run)
		table$community <- factor(table$community)
		
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
		message(paste("Printing Sum of Squares for subtable n°",i, sep="")) 
		setwd(paste(WD,"/","ANOVA_tables","/","Emission_scenarios","/","SSQ_15_12","/", sep="")) 
		save(sum_squares, file = paste("ssq_dSR_subtable2_",i,".Rdata", sep=""))
		
		# Make room:
		rm(table, sum_squares)
		gc()
		setwd(paste(WD,"/","ANOVA_tables","/", sep=""))
		
} # eo for loop		




