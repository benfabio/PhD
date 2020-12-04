
##### 24/02/2017 - LOV - Fabio Benedetti
##### Script for : 
#		- Assessing main results' sentivity to varying parameters : data omission, ENM choice and smoothing parameter choice by:
#		- Loading the several niche traits tables + functional traits table
#		- Plotting species in niche space, based on the selected set of niche traits
#		- Plotting FG in the same niche space and test PC1 variation between them
#		- Also check average TSS variations


### Last update : 24/02/2017

# -----------------------------------------------------------------------------------------------------------------------------

library("sp")
library("dplyr")
library("stringr")
library("reshape2")
library("fields")
library("scales")
library("ggplot2")
library("RColorBrewer")
library("corrgram")
library("FactoMineR")
library("autoplot")

WD <- getwd()

jet.colors <- c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")

# -----------------------------------------------------------------------------------------------------------------------------


##### A) Load functional traits data
ddf <- read.csv("table_traits_biogeography.csv", h = T, dec = ",", sep = ";")
# summary(ddf)
ddf$fun_group <- factor(ddf$fun_group)
# Give size classes
distance <- dist((ddf[,c("length_max")]), method = "euclidean") # distance matrix
fit <- hclust(distance, method= "ward") 
# clustering on the distance matrix computed from the CA, aggregation method = Ward's
# plot(fit, hang= -1) # display dendogram
groups <- cutree(fit, k = 4) # cut tree into 4 clusters
# groups
ddf$size_class <- factor(groups)
#																		Median	 Mean
# summary(ddf[which(ddf$size_class == 1),"length_max"]) # 1.300   1.407   1.555   1.548   1.660   1.800
# summary(ddf[which(ddf$size_class == 2),"length_max"]) # 1.890   2.000   2.240   2.304   2.550   3.030
# summary(ddf[which(ddf$size_class == 3),"length_max"]) # 3.400   3.925   4.550   4.824   5.198   8.250
# summary(ddf[which(ddf$size_class == 4),"length_max"]) # 0.5700  0.7600  0.8600  0.8619  0.9600  1.2000

### Restrict to functional traits columns 
ddf <- ddf[,c(1,65:72)]




##### B) Use MCA to re-define functional groups
require("FactoMineR")
ddf$C <- (ddf$trophism == "Carnivore")|(ddf$trophism == "Omnivore-Carnivore")
ddf$O <- (ddf$trophism == "Omnivore")|(ddf$trophism == "Omnivore-Carnivore")|(ddf$trophism == "Omnivore-Detritivore")|(ddf$trophism == "Omnivore-Herbivore")
ddf$D <-(ddf$trophism == "Omnivore-Detritivore")
ddf$H <-(ddf$trophism == "Omnivore-Herbivore")

ddf$na <- is.na(ddf$trophism) + is.na(ddf$feeding) + is.na(ddf$spawning) + is.na(ddf$DVM) + is.na(ddf$layer)
ddf$na3 <- is.na(ddf$trophism) + is.na(ddf$feeding) + is.na(ddf$spawning)
known <- which(ddf$na3 == 0)
unknown <- which(ddf$na3 > 0)
rownames(ddf) <- ddf$sp.name 
dataf.MCA <- MCA(ddf[,c("size_class","feeding","spawning","C","O","D","H")], ncp = 4, na.method="Average", ind.sup = unknown)
#dataf.MCA <- MCA(ddf[,c("size_class","feeding","C","O","D","H")], ncp=4, na.method="Average", ind.sup = unknown)
#dataf.MCA <- MCA(ddf[,c("size_class","C","O","D","H")], ncp=3, na.method="Average", ind.sup = unknown)
#dataf.MCA <- MCA(ddf[,c("size_class","spawning","C","O","D","H")], ncp=4, na.method="Average", ind.sup = unknown)
summary(dataf.MCA)
dataf.MCA$var$contrib
#                     Dim 1        Dim 2      Dim 3        Dim 4
# size_class.1   0.63665401  1.958099075  2.8521745 44.047430304
# size_class.2   2.61019770  4.617668794  7.9832014 12.775615025
# size_class.3   0.00795111  5.814286099 31.5767475  0.005465652
# size_class.4   1.02939458  9.276813739  0.9163114 11.690656358
# Active Ambush  8.48616339  1.555657911 12.1373869  2.203213220
# Cruise         1.58223467 12.196609260  4.5412704  0.511903040
# Filter         7.83297694  3.798000321 11.7128782  0.076995848
# Mixed          3.09582820  0.283006874 20.3020808  2.788275612
# Broadcaster   11.31091457  3.888249416  0.1757640  2.774856968
# Sac-spawner    9.31487082  3.202087755  0.1447468  2.285176326
# C_FALSE        3.28957833  3.079889871  0.2044773  0.066022705
# C_TRUE        13.70657636 12.832874464  0.8519886  0.275094603
# O_FALSE       13.70657636 12.832874464  0.8519886  0.275094603
# O_TRUE         3.28957833  3.079889871  0.2044773  0.066022705
# D_FALSE        0.22243279  2.436201567  0.1851246  0.974402828
# D_TRUE         1.74768620 19.141583738  1.4545507  7.656022220
# H_FALSE        7.89549053  0.002702953  1.7004910  5.020150057
# H_TRUE        10.23489513  0.003503828  2.2043401  6.507601926

# Looking at the eigenvalues
eig <- data.frame(prop=dataf.MCA$eig$"percentage of variance", nb= c(1:11))
quartz()
ggplot(eig) + geom_bar(aes(x=nb, y=prop), stat="identity") + geom_line(aes(x=nb, y=mean(prop) )) + xlab("Eigenvalue Number") +  ylab("Value") + theme_bw()
### Keep first 4 MCA axes

# Getting the results
mca1 <- paste0("MCA 1 (",floor(eig$prop[1]*100)/100,"%)")
mca2 <- paste0("MCA 2 (",floor(eig$prop[2]*100)/100,"%)")
mca3 <- paste0("MCA 3 (",floor(eig$prop[3]*100)/100,"%)")
mca4 <- paste0("MCA 4 (",floor(eig$prop[4]*100)/100,"%)")

#mca_var <- data.frame(dataf.MCA$var$coord, Variable = rep(names(cats), cats),Class=row.names(dataf.MCA$var$coord))
mca_sp <- data.frame(dataf.MCA$ind$coord)
mca_sp_sup <- data.frame(dataf.MCA$ind.sup$coord)

# colnames(mca_var) <- c("MCA1","MCA2","MCA3","MCA4","MCA5","MCA6",#"MCA7","MCA8","MCA9", "Variable","Class")
colnames(mca_sp) <- c("MCA1","MCA2","MCA3","MCA4")
colnames(mca_sp_sup) <- c("MCA1","MCA2","MCA3","MCA4")

### Gather all spp (3 were left as suppl objects)
mca_all_temp <- rbind(mca_sp, mca_sp_sup)
mca_all_sp <- mca_all_temp[ order(row.names(mca_all_temp)), ]
dim(mca_all_sp) ; dim(mca_all_temp)
  
dist_fct <- dist(mca_all_sp, method = "euclidean")
# mca_distance <- dist(mca_sp[,1:DIM], method = "euclidean") # distance matrix
fit_mca <- hclust(dist_fct, method= "ward") 
quartz()
plot(fit_mca, hang= -1, labels = rownames(mca_all_sp)) # display dendogram with names

#phylo <- as.phylo(fit_mca)
#plot(phylo)
#save(phylo, file = "fct_dendro.Rdata")

### Vector of new functional groups
groups <- cutree(fit_mca, k = 7) # cut tree into 4 clusters
# Supply to 'ddf'
ddf$fun_group_2 <- factor(groups)
str(ddf)

# New column: 
ddf$FG <- NA
ddf[which(ddf$fun_group_2 == 1), "FG"] <- 4
ddf[which(ddf$fun_group_2 == 2), "FG"] <- 2
ddf[which(ddf$fun_group_2 == 3), "FG"] <- 3
ddf[which(ddf$fun_group_2 == 4), "FG"] <- 7
ddf[which(ddf$fun_group_2 == 5), "FG"] <- 1
ddf[which(ddf$fun_group_2 == 6), "FG"] <- 5
ddf[which(ddf$fun_group_2 == 7), "FG"] <- 6



##### C) Load environmental niche traits data (8 tables) + use PCA to make niche space + plot spp. + FG in niche space 
tables <- dir()[grep("23_02_17", dir())]

# t <- tables[8] ### for testing

# Perform analyses and produce plots for each of the 8 tables
for(t in tables) {
	
		message(paste("ANALYZING table ", t, sep = ""))
		
		# Load niche traits table
		niche <- read.table(t, h = T, sep = ";", dec = ".")
		
		# Keep an index to inform which parameter you are dealing with:
		para <- strsplit(t, "_")[[1]][3]
		
		# Change rownames
		rownames(niche) <- niche[,"sp.name"]
	
		# If t == "niche_traits_BRT_23_02_17.txt" or "niche_traits_MAXENT_23_02_17.txt", get rid of the few spp you couldn't model with logChla or MLPAR1
		if( t == "niche_traits_BRT_23_02_17.txt") {
				niche <- niche[-c(17,37,70,85,87,94),]
		} else if (t == "niche_traits_MAXENT_23_02_17.txt") {
				niche <- niche[-c(85),]
		} # eo if else loop
	
		# Select niche traits
		niche <- niche[,c("sp.name","center_SST","breadth_SST","upper_dSST","lower_dSST","center_SSS","center_MLD1","center_MLPAR1","breadth_MLPAR1","upper_MLPAR1","center_logChla","breadth_logChla")]	
		
		# Perform PCA
		res.pca <- PCA(niche[,-c(1)], scale.unit = TRUE, ncp = 3, graph = T, axes = c(1,2))
		# res.pca <- PCA(table, scale.unit = TRUE, ncp = 3, graph = T, axes = c(2,3))
		# summary(res.pca)
		
		### Plot a nice PCA space (variables)
		RES <- augment(res.pca, dimensions= c(1,2), which = "col")
		# str(RES)
		# RES
		colnames(RES)[1:6] <- c("vars", "type", "PC1", "PC2", "cos2", "contrib")
		### For drawing correlation circle: ! 
		# http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
		circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){		   
		    	r = diameter / 2
				tt <- seq(0,2*pi,length.out = npoints)
		    	xx <- center[1] + r * cos(tt)
		    	yy <- center[2] + r * sin(tt)
		    	return(data.frame(x = xx, y = yy))		
		} # eo circleFun
		dat <- circleFun(c(0,0), diameter = 2, npoints = 100)

		plot1 <- ggplot() + 
				geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
				geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2, label= vars), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour= "firebrick3", data= RES) +
				geom_path(data = dat, aes(x,y) ) + geom_text(data = RES, aes(x = PC1, y = PC2, label = vars)) + 
				xlab( paste("PC1 (", round(res.pca$eig$per[1], 2), "%)", sep = "") ) + ylab( paste("PC2 (", round(res.pca$eig$per[2], 2), "%)", sep = "") ) + theme_linedraw()

		# quartz()
		# plot1
		ggsave(paste("plot_pca_vars_1&2_",para,".pdf", sep = ""), plot = plot1, dpi = 300, width = 7, height = 7)
		gc() ; rm(plot1, dat, RES)
		
		### Plot a nice PCA space (objects)
		niche$PC1_ind <- res.pca$ind$coord[,1]
		niche$PC2_ind <- res.pca$ind$coord[,2]
		niche$PC3_ind <- res.pca$ind$coord[,3]
		# niche$index <- c(1:106)
		
		library("ggrepel")
		plot1 <- ggplot() + geom_point(aes(x = PC1_ind, y = PC2_ind, fill = PC1_ind), colour = "black", size = 3, data = niche, shape = 21) + 
			 		#geom_segment(aes(x=0, y=0, xend= PC1*5, yend= PC2*5, label= vars), arrow = arrow(length = unit(0.4, "cm"),type = "closed"), colour = "grey50", data = vars) + 
					scale_fill_gradient2(name = "Coordinate on PC1", guide = "colorbar", low = "blue", high = "red") + 
					#scale_shape_manual(name = "Niche traits\ngroups", values = c(21,22,23,24)) + 
			 		geom_text_repel( aes(x = PC1_ind, y = PC2_ind, label = sp.name), size = 2.75, data = niche) +
					#geom_text(aes(x= PC1*5, y= PC2*5, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
		 		   	geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
					xlab( paste("PC1 (", round(res.pca$eig$per[1], 2), "%)", sep = "") ) + ylab( paste("PC2 (", round(res.pca$eig$per[2], 2), "%)", sep = "") ) + 
					theme_bw()
		# quartz()
		# plot1
		ggsave(paste("plot_pca_sp_1&2_", para, ".pdf", sep = ""), plot = plot1, dpi = 300, width = 15, height = 7)
		gc() ; rm(plot1)
		
		### Now, plot FG in niche space
		# Supplying functional groups to niche space 
		niche$FG <- NA
		for(i in rownames(niche)) {
			niche[i,"FG"] <- ddf[i,"FG"]
		}
		# niche
		
		# ANOVA:
		### Variation of functional groups across PC1 & 2
		fit <- aov(PC1_ind ~ factor(FG), data = niche)
		# shapiro.test(fit$residuals) # p-value = 1.839e-08 ; not normally distributed
		# summary(fit) # p-value = 0.00247 ! ANOVA shows significant variation of fun group in PC1
		### If fligner$p.value is > 0.05 and KW$p.value is < 0.05, then save last plot
		fligner_pval <- fligner.test(PC1_ind ~ factor(FG),  data = niche)$p.value
		KW_pval <- kruskal.test(PC1_ind ~ factor(FG), data = niche)$p.value
		# Save them somewhere
		writeLines(as.character(fligner_pval), paste("fligner_pval_", para, ".txt", sep = ""))
		writeLines(as.character(KW_pval), paste("KW_pval_", para, ".txt", sep = ""))
			
		### Prepare dataframes for plotting FG in niche space --> need to compute groiups' niche centroids
		bary <- data.frame(FG = unique(ddf$fun_group_2), PC1 = NA, PC2 = NA, PC3 = NA)
		#
		bary[which(bary$FG == 1), "PC1"] <- mean(niche[which(niche$FG == 1), "PC1_ind"] )
		bary[which(bary$FG == 1), "PC2"] <- mean(niche[which(niche$FG == 1), "PC2_ind"] )
		bary[which(bary$FG == 1), "PC3"] <- mean(niche[which(niche$FG == 1), "PC3_ind"] )
		#
		bary[which(bary$FG == 2), "PC1"] <- mean(niche[which(niche$FG == 2), "PC1_ind"] )
		bary[which(bary$FG == 2), "PC2"] <- mean(niche[which(niche$FG == 2), "PC2_ind"] )
		bary[which(bary$FG == 2), "PC3"] <- mean(niche[which(niche$FG == 2), "PC3_ind"] )
		#
		bary[which(bary$FG == 3), "PC1"] <- mean(niche[which(niche$FG == 3), "PC1_ind"] )
		bary[which(bary$FG == 3), "PC2"] <- mean(niche[which(niche$FG == 3), "PC2_ind"] )
		bary[which(bary$FG == 3), "PC3"] <- mean(niche[which(niche$FG == 3), "PC3_ind"] )
		#
		bary[which(bary$FG == 4), "PC1"] <- mean(niche[which(niche$FG == 4), "PC1_ind"] )
		bary[which(bary$FG == 4), "PC2"] <- mean(niche[which(niche$FG == 4), "PC2_ind"] )
		bary[which(bary$FG == 4), "PC3"] <- mean(niche[which(niche$FG == 4), "PC3_ind"] )
		#
		bary[which(bary$FG == 5), "PC1"] <- mean(niche[which(niche$FG == 5), "PC1_ind"] )
		bary[which(bary$FG == 5), "PC2"] <- mean(niche[which(niche$FG == 5), "PC2_ind"] )
		bary[which(bary$FG == 5), "PC3"] <- mean(niche[which(niche$FG == 5), "PC3_ind"] )
		#
		bary[which(bary$FG == 6), "PC1"] <- mean(niche[which(niche$FG == 6), "PC1_ind"] )
		bary[which(bary$FG == 6), "PC2"] <- mean(niche[which(niche$FG == 6), "PC2_ind"] )
		bary[which(bary$FG == 6), "PC3"] <- mean(niche[which(niche$FG == 6), "PC3_ind"] )
		#
		bary[which(bary$FG == 7), "PC1"] <- mean(niche[which(niche$FG == 7),"PC1_ind"] )
		bary[which(bary$FG == 7), "PC2"] <- mean(niche[which(niche$FG == 7),"PC2_ind"] )
		bary[which(bary$FG == 7), "PC3"] <- mean(niche[which(niche$FG == 7),"PC3_ind"] )
				
		### Plot FG in niche space
		plot1 <- ggplot() + 
					geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
					geom_point(aes(x = PC1_ind, y = PC2_ind, shape= factor(FG), fill = factor(FG)), data = niche, colour = "black", size = 2) + 
					#geom_segment(aes(x=0, y=0, xend=PC1*4, yend=PC2*4, label= vars), arrow = arrow(length= unit(0.1, "cm"), type= "closed"), colour= "black", data= vars) +
					geom_point(aes(x = PC1, y = PC2, shape= factor(FG), fill = factor(FG)), colour = "black", data = bary, size = 5) + 
					#scale_fill_manual(name = "Functional\ngroups", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")   ) +
					scale_fill_manual(name = "Functional group", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")  ) +
					scale_shape_manual(name= "Functional group", values = c(21,22,23,24,25,21,22)) + 
					#geom_text(aes(x= PC1*4, y= PC2*4, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
					xlab( paste("PC1 (",round(res.pca$eig$per[1],2), "%)", sep = "") ) + ylab( paste("PC2 (",round(res.pca$eig$per[2],2), "%)", sep = "") ) + 
					theme_bw()
						
		# quartz()
		# plot1
				
		ggsave(paste("plot_pca_FG_1&2_", para, ".pdf", sep = ""), plot = plot1, dpi = 300, width = 10, height = 7)
		
		plot2 <- ggplot() + 
					geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
					geom_point(aes(x= PC2_ind, y= PC3_ind, shape= factor(FG), fill= factor(FG)), data = niche, colour = "black", size = 2) + 
					#geom_segment(aes(x=0, y=0, xend=PC2*4, yend=PC3*4, label= vars), arrow = arrow(length= unit(0.1, "cm"), type= "closed"), colour= "black", data= vars) +
					geom_point(aes(x=PC2, y= PC3, shape= factor(FG), fill= factor(FG)), colour= "black", data = bary, size = 5) + 
					#scale_fill_manual(name = "Functional\n groups", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")   ) +
					scale_fill_manual(name = "Functional group", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")  ) +
					scale_shape_manual(name= "Functional group", values = c(21,22,23,24,25,21,22)) + 
					#geom_text(aes(x= PC2*4, y= PC3*4, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
					xlab( paste("PC2 (",round(res.pca$eig$per[2],2), "%)", sep = "") ) + ylab( paste("PC3 (",round(res.pca$eig$per[3],2), "%)", sep = "") ) + 
					theme_bw()

		# quartz()
		# plot2
				
		ggsave(paste("plot_pca_FG_2&3_", para, ".pdf", sep = ""), plot = plot2, dpi = 300, width = 10, height = 7)
				
		gc() 
		rm(bary, plot1, plot2, niche, KW_pval, fligner_pval, fit)
		# setwd(WD)


} # eo FOR LOOP





