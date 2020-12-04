
##### 16/05/2016 - ETHZ - Fabio Benedetti
##### Script for : 
#			- Plot species in niche space, based on the 11 dimensions selected from the correlogram
#			- Re-run analysis of funct. groups in this niche space
#			- Assess sensitivity of functional groups definition by removing trauts from MCA space
#			- Re-run cophenetic correlation


### Last update : 10/11/2017

# -------------------------------------------------------------------------------------------------------------

library("sp")
#library("plyr")
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
library("ape")
library("picante")

WD <- getwd()

jet.colors <- c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")

# -------------------------------------------------------------------------------------------------------------

##### A) Correlogram & PCA to build environmental/niche space 
ddf <- read.csv("table_traits_biogeography.csv", h=T, dec=",", sep=";")
# summary(ddf)
ddf$fun_group <- factor(ddf$fun_group)
# Give size classes
distance <- dist((ddf[,c("length_max")]), method = "euclidean") # distance matrix
fit <- hclust(distance, method= "ward") 
# clustering on the distance matrix computed from the CA, aggregation method = Ward's
# plot(fit, hang= -1) # display dendogram
groups <- cutree(fit, k=4) # cut tree into 4 clusters
# groups
ddf$size_class <- factor(groups)
#																		Median	 Mean
summary(ddf[which(ddf$size_class == 1),"length_max"]) # 1.300   1.407   1.555   1.548   1.660   1.800
summary(ddf[which(ddf$size_class == 2),"length_max"]) # 1.890   2.000   2.240   2.304   2.550   3.030
summary(ddf[which(ddf$size_class == 3),"length_max"]) # 3.400   3.925   4.550   4.824   5.198   8.250
summary(ddf[which(ddf$size_class == 4),"length_max"]) # 0.5700  0.7600  0.8600  0.8619  0.9600  1.2000

# 4 < 1 < 2 < 3

### Spearman's correlation rank
# quartz()
# corrgram(ddf[,c(2:5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65)], 
# order= TRUE, lower.panel= panel.shade, upper.panel= panel.pie, text.panel= panel.txt, cor.method = "spearman") 
### SEE SCRIPT9.9 FOR BETTER LOOKING CORRELOGRAM

### Based on spearman' correlations, we're keeping the follwoing niche traits to build the niche space:
table <- ddf[,c("center_SST","breadth_SST","upper_dSST","lower_dSST","center_SSS","center_MLD1",
				"center_MLPAR1","breadth_MLPAR1","upper_MLPAR1","center_logChla","breadth_logChla")]				
rownames(table) <- ddf[,"sp.name"]

# Principal Component Analysis
#quartz()
res.pca <- PCA(table, scale.unit = TRUE, ncp = 3, graph = T, axes = c(1,2))
#res.pca <- PCA(table, scale.unit = TRUE, ncp = 3, graph = T, axes = c(2,3))
summary(res.pca)

### NO DIFFERENCES WHEN EXCLUDING SPP WITH NO TROPHISM (n = 5)
# [1] Copilia_quadrata      Lubbockia_squillimana 	Monothula_subtilis   
# [4] Neomormonilla_minor   Phaenna_spinifera 

### Variables' contribution
# res.pca$var$contrib
#                       Dim.1     	Dim.2      	Dim.3
# center_SST      15.436528266  	2.712200  	2.080729224
# breadth_SST      3.344354912 		26.316745  	2.290959041
# upper_dSST      10.741381482  	3.432873  	5.709559321
# lower_dSST       0.001886183 		19.594769 	14.394307005
# center_SSS      20.001916074  	1.187824  	0.009354675
# center_MLD1      5.526586606  	2.079869  	5.298041669
# center_MLPAR1   19.943395687  	3.853846  	0.389301606
# breadth_MLPAR1   0.280270740  	3.968663 	48.536780243
# upper_MLPAR1    12.497602948  	7.598016 	16.703132861
# center_logChla  12.098392948  	1.331801  	1.898366910
# breadth_logChla  0.127684154 		27.923394  	2.689467446

### PC1 : center_SSS; center_SST; center_MLPAR; upper_MLPAR VERSUS center_logChla; upper_dSST; center_MLD to a lesser extent
### --> trophic gradient; from warm stable oceanographic conditions, to colder, fresher, more seasonally variable and more productive environments

### PC2 : breadth_logChla; breadth_SST VERSUS lower_dSST
### --> Species with broader thermal and Chla niches VERSUS species exhibiting narrower niches for these 2 dimensions 
#      (whether they are cold-water, warm-water, productive waters- affiliated)

### PC3 : breadth_MLPAR1; upper_MLPAR1; lower_dSST

### Better PCA plots ? 
library("autoplot")
RES <- augment(res.pca, dimensions = c(1,2), which = "col")
str(RES)
RES
colnames(RES)[1:6] <- c("vars", "type", "PC1", "PC2", "cos2", "contrib")
### For drawing correlation circle: ! 
# http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
dat <- circleFun(c(0,0), diameter = 2, npoints = 100)

plot1 <- ggplot() + theme_light() + 
		geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
		geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2 ), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour= "firebrick3", data = RES) +
		geom_path(data = dat, aes(x,y) ) + xlab("PC1 (35.7%)") + ylab("PC2 (21.7%)")

quartz()
plot1
ggsave("Fig.3.pdf", plot = plot1, dpi = 300, width = 5, height = 5)




### get species coordinates in niche space
res.pca$ind$coord
table$PC1_ind <- res.pca$ind$coord[,1]
table$PC2_ind <- res.pca$ind$coord[,2]
table$PC3_ind <- res.pca$ind$coord[,3]

### get variables coordinates in another table
res.pca$var$coord
vars <- data.frame(vars = colnames(table)[1:11],
				PC1 = res.pca$var$coord[,1],
				PC2 = res.pca$var$coord[,2],
				PC3 = res.pca$var$coord[,3])


# HAC from the species coordinates (PC1-3)
distance <- dist((table[,12:14]), method = "euclidean") # distance matrix
fit <- hclust(distance, method= "ward") 

names <- stringr::str_replace_all(rownames(table), "_", " ")
table$names <- names

##### For plotting a dendrogram of niche space...
library("ggdendro")
dendr    <- dendro_data(fit, type="rectangle") # convert for ggplot
clust    <- cutree(fit, k=4)
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
# dendr[["labels"]] has the labels, merge with clust.df based on label column
dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by="label")
dendr$label$label <- names

dendro <- ggplot() + geom_segment(data= segment(dendr), aes(x= x, y= y, xend= xend, yend= yend)) + 
  		   geom_text(data= label(dendr), aes(x, y, label= str_replace_all(label, "_", " "), hjust= 0), size= 3, face= "italic") +
  		   coord_flip() + scale_y_reverse(expand= c(0.2, 0)) + 
  		   theme_bw()
quartz()
dendro
# ggsave(filename= "dendro_niche_groups.pdf", plot = dendro, dpi = 300)		   


table$index <- c(1:106)
table$sp.names <- ddf$sp.name 
table$sp.names <- str_replace_all(table$sp.names, "_", " ")

library("ggrepel")

plot1 <- ggplot() + geom_point(aes(x= PC1_ind, y= PC2_ind, fill= PC1_ind), colour = "black", size = 3, data = table, shape = 21) + 
	 		#geom_segment(aes(x=0, y=0, xend= PC1*5, yend= PC2*5, label= vars), arrow = arrow(length = unit(0.4, "cm"),type = "closed"), colour = "grey50", data = vars) + 
			scale_fill_gradient2(name="Coordinate on PC1", guide = "colorbar", low = "blue", high = "red") + 
			#scale_shape_manual(name = "Niche traits\ngroups", values = c(21,22,23,24)) + 
	 		geom_text_repel( aes(x= PC1_ind, y= PC2_ind, label= sp.names), size = 2.75, data = table) +
			#geom_text(aes(x= PC1*5, y= PC2*5, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
 		   	geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
			xlab("PC1 (35.68%)") + ylab("PC2 (21.74%)") + 
			theme_bw()
quartz()
plot1


plot2 <- ggplot() + geom_point(aes(x= PC2_ind, y= PC3_ind, fill= PC2_ind), colour = "black", size = 3, data = table, shape = 21) + 
	 		#geom_segment(aes(x=0, y=0, xend= PC1*5, yend= PC2*5, label= vars), arrow = arrow(length = unit(0.4, "cm"),type = "closed"), colour = "grey50", data = vars) + 
			scale_fill_gradient2(name="Coordinate on PC2", guide = "colorbar", low = "blue", high = "red") + 
			#scale_shape_manual(name = "Niche traits\ngroups", values = c(21,22,23,24)) + 
	 		geom_text_repel( aes(x= PC2_ind, y= PC3_ind, label= sp.names), size = 2.75, data = table) +
			#geom_text(aes(x= PC1*5, y= PC2*5, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
 		   	geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
			xlab("PC2 (21.74%)") + ylab("PC2 (14.85)%)") + 
			theme_bw()
quartz()
plot2
					
			
ggsave(filename="pca_biplot_1.pdf", plot = plot1, dpi = 300)
ggsave(filename="pca_biplot_2.pdf", plot = plot2, dpi = 300)





##### B) MCA to build functional space and define functional groups
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
ggplot(eig) + geom_bar(aes(x=nb, y=prop), stat="identity") + geom_line(aes(x=nb, y=mean(prop) )) + xlab("Eigenvalue Number") +  ylab("Value")
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
groups <- cutree(fit_mca, k = 7) # cut tree into 7 clusters
# Supply to 'ddf'
ddf$fun_group_2 <- groups

### Rename fun_groups properly : 
head( ddf[which(ddf$fun_group_2 == 1),"sp.name"] ) # 1 --> FG 4
head( ddf[which(ddf$fun_group_2 == 2),"sp.name"] ) # 2 --> FG 2
head( ddf[which(ddf$fun_group_2 == 3),"sp.name"] ) # 3 --> FG 3
head( ddf[which(ddf$fun_group_2 == 4),"sp.name"] ) # 4 --> FG 7
head( ddf[which(ddf$fun_group_2 == 5),"sp.name"] ) # 5 --> FG 1
head( ddf[which(ddf$fun_group_2 == 6),"sp.name"] ) # 6 --> FG 5
head( ddf[which(ddf$fun_group_2 == 7),"sp.name"] ) # 7 --> FG 6

# New column: 
ddf$FG <- NA
ddf[which(ddf$fun_group_2 == 1), "FG"] <- 4
ddf[which(ddf$fun_group_2 == 2), "FG"] <- 2
ddf[which(ddf$fun_group_2 == 3), "FG"] <- 3
ddf[which(ddf$fun_group_2 == 4), "FG"] <- 7
ddf[which(ddf$fun_group_2 == 5), "FG"] <- 1
ddf[which(ddf$fun_group_2 == 6), "FG"] <- 5
ddf[which(ddf$fun_group_2 == 7), "FG"] <- 6

# Niche characteristics
mean(ddf[which(ddf$FG == 1), "upper_logChla"])
mean(ddf[which(ddf$FG == 2), "upper_logChla"])	
mean(ddf[which(ddf$FG == 3), "upper_logChla"])	
mean(ddf[which(ddf$FG == 4), "breadth_logChla"])		
mean(ddf[which(ddf$FG == 5), "upper_logChla"])		
mean(ddf[which(ddf$FG == 6), "upper_logChla"])		
mean(ddf[which(ddf$FG == 7), "upper_logChla"])		


### Do another dendrogram, like above, with ggdendro package
library("ggdendro")
dendr2    <- dendro_data(fit_mca, type="rectangle") # convert for ggplot
clust2    <- cutree(fit_mca, k=7)
clust2.df <- data.frame(label = names(clust2), cluster = factor(clust2), FG = NA)
clust2.df[which(clust2.df$cluster == 1), "FG"] <- 4
clust2.df[which(clust2.df$cluster == 2), "FG"] <- 2
clust2.df[which(clust2.df$cluster == 3), "FG"] <- 3
clust2.df[which(clust2.df$cluster == 4), "FG"] <- 7
clust2.df[which(clust2.df$cluster == 5), "FG"] <- 1
clust2.df[which(clust2.df$cluster == 6), "FG"] <- 5
clust2.df[which(clust2.df$cluster == 7), "FG"] <- 6

# dendr[["labels"]] has the labels, merge with clust.df based on label column
dendr2[["labels"]] <- merge(dendr2[["labels"]], clust2.df, by="label")
dendr2$label <- ddf$sp.names

dendro2 <- ggplot() + geom_segment(data=segment(dendr2), aes(x= x, y= y, xend= xend, yend= yend)) + 
  		   geom_text(data= label(dendr2), aes(x=x, y=y, label = str_replace_all(label, "_", " "), hjust=0, color= factor(FG)), size = 3 , fontface = "italic") +
		   scale_colour_manual(name = "Functional\ngroups", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63") ) +
  		   coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + geom_hline(y= 7, type = "dashed") + 
  		   theme_classic() 

quartz()
dendro2		   
ggsave(filename= "dendro_fct_groups_k7.pdf", plot = dendro2, height= 16, width= 14)		   


### 31/05/2016: Examine the species traits within the 7 functional groups (FG)
#colnames(ddf)
ddf[which(ddf$FG == 1),c(65:70,72,80)]  # n = 10 ; 9.4% ; mean(ddf[which(ddf$FG == 1),c("length_max")])
# Large Carnivores (+2 Omni-Carni) ; mainly sac-spawning ; no leading feeding strategy
# Genera: Euchaeta + Pareuchaeta (Euchaetidae), 1 Sapphirina, 2 Haloptilus, 1 Heterorhabdus, 1 Candacia
# Indicative: Both strong and weak migrants, no driving water layer...

ddf[which(ddf$FG == 2),c(65:70,72,80)] # n = 15 ; 14% des spp
# Smaller and strictly Carnivores, mainly ambush feeding (but 50% NaNs) and 50/50 between sac-spawning and broadcasting,
# main Taxa: Candacia and Corycaeus, Aetideus, 1 Heterorhabdus, 1 Augaptilus
# Indicative: Almost all epi-mesopelagics, generally weak and non-migrants (but C. furcifer)

ddf[which(ddf$FG == 3),c(65:70,72,80)] # n = 8 ; 7.5% des spp
# Large (all from largest SC) filter-feeding Herbivores, 5 Broadcasters (3 NaNs)
# Taxa: only large calanoids: Calanus, Eucalanus, Mesocalanus, Neocalanus, 1 Pleuromamma, Subeucalanus...open ocean taxa, Atlantic origins...
# Indicative: Either strong or non-migrant (lol), epi to deeper layers

ddf[which(ddf$FG == 4),c(65:70,72,80)] # n = 27 ; 25% des espèces
# Smaller (SC 4 to 2, but mostly intermediate SC) filter feeding Herbivores + mixed feeders that are Omnivores (Acartia spp + Centropages spp) ; ALL broadcasters (8 NaNs)
# Taxa: Calanoids, with some that are numerically important in the Med: Acartia + Centropages, Calocalanus spp., 1 Clausocalanus, 1 Ctenocalanus, 2 Pleuromamma, Lucicutia spp., 
# 2 Paracalanus, T. stylifera, Subeucalanus monachus, Nannocalanus minor...
# Indicative traits : 50% of Epipelagic species that are weak and non-migrants...plus some that can cross different layers, 3 strong migrants (Pleuromamma and Subeucalanus)

ddf[which(ddf$FG == 5),c(65:70,72,80)] # n = 21 (19%)
# Small cruise feeding and sac-spawning (+ 4 filter feeding) Detritivores + some herbivores that were clustered with them because of Cruise feeding
# Taxa: both calanoids (Clausocalanus, Scaphocalanus, and Scolecithrichidae) and cyclopoids (Oncaea, Triconia), and Harpacticoids (Microsetella spp.)
# Indicative traits: mostly weak migrants that inhabit the epimesopelagic

ddf[which(ddf$FG == 6),c(65:70,72,80)] # n = 12 (11%)
# Oithonids: small active ambush omnivores, sac-spawners. + Species with many NaNs (Isias, Lubbockia and 1 lost Haloptilus)

ddf[which(ddf$FG == 7),c(65:70,72,80)] # n = 13 (12%)
# Small (SC == 1 only) sac-spawning herbivores, mostly filter feeding (+ 4 cruise feeding Clausocalanus spp.) + 4 NaNs in terms of feeding strat
# Similar to FG n°4, but separated because of egg-spawning strat... 
# Taxa: Clausocalanus spp., Paracalanus parvus, Temora longicornis, Pseudocalanus elongatus...species with many NaNs: Neomormonilla minor, M.clausi
# Indicative traits: mostly epipelagic species, weak and non-migrants





##### C) Plotting the 7 functional groups in defined niche space + ANOVA to assess variance of FGs' coordinates along PCs

### Supplying functional groups to niche space 
table$fun_group_2 <- NA
table$fun_group_2 <- ddf[rownames(table), "FG"]

### Variation of functional groups across PC1 & 2
fit <- aov(PC1_ind ~ factor(fun_group_2), data = table)
shapiro.test(fit$residuals) # p-value = 1.839e-08 ; not normally distributed
summary(fit) # p-value = 0.00247 ! ANOVA shows significant variation of fun group in PC1
fligner.test(PC1_ind ~ factor(fun_group_2),  data = table)
# p-value = 0.6405 ; homogeneity of variances
kruskal.test(PC1_ind ~ factor(fun_group_2), data = table)
# p-value = 0.007553 --> Significant variations of FG of PC1 !!

														 #   Min.					Median		 Mean					  	Max
summary(table[which(table$fun_group_2 == 1),"PC1_ind"])	 # -2.1000  	0.1027  	0.6812  	0.5564  	1.2630  		2.2540 
summary(table[which(table$fun_group_2 == 2),"PC1_ind"])	 # -2.3550  	0.8165  	1.2000  	1.0280  	1.9600  		2.5620 
summary(table[which(table$fun_group_2 == 3),"PC1_ind"])	 # -3.6190 		-1.9820 	-1.1140 	-1.3110 	-0.2379  		0.4937
summary(table[which(table$fun_group_2 == 4),"PC1_ind"])	 # -4.7720 		-0.8160  	0.5698  	0.1783  	1.3070  		2.7120 
summary(table[which(table$fun_group_2 == 5),"PC1_ind"])  # -5.4610 		-0.2155  	0.8673  	0.2255  	1.5340  		2.0950 
summary(table[which(table$fun_group_2 == 6),"PC1_ind"])  # -7.1320 		-1.2090  	0.6522 		-0.4408  	0.9524  		1.8330
summary(table[which(table$fun_group_2 == 7),"PC1_ind"])  # -7.1780 		-1.6480 	-0.4661 	-1.1350  	0.9861  		1.3670
### Like before: Carnivores and Omni-Carni have warmer/ saltier niches
# pairwise.wilcox.test(table$PC1_ind, factor(table$trophism), p.adjust.method = "BH" )


###### Nothing going on with PC2 though.


### But: PC3
fligner.test(PC3_ind ~ factor(fun_group_2),  data = table)
# p-value = 0.1421 ; homogeneity of variances
kruskal.test(PC3_ind ~ factor(fun_group_2), data = table)
# p-value = 0.007329 ! 
														 #   Min.					Median		 Mean					  	Max
summary(table[which(table$fun_group_2 == 1),"PC3_ind"])	 # -1.62600 	-1.12000  	0.06963 	-0.04926  	0.67630  		2.09000 
summary(table[which(table$fun_group_2 == 2),"PC3_ind"])	 # -1.4000 		-0.1034  	0.6038  	0.5815  	1.0060  		2.8490 
summary(table[which(table$fun_group_2 == 3),"PC3_ind"])	 # -0.2670  	0.0998  	1.2020  	1.1770  	2.0930  		2.7650 
summary(table[which(table$fun_group_2 == 4),"PC3_ind"])	 # -2.5370 		-0.2478  	0.3475  	0.2262  	0.7383  		2.5960 
summary(table[which(table$fun_group_2 == 5),"PC3_ind"])  # -1.6720 		-0.7143 	-0.1424 	-0.1844  	0.3591  		1.1190 
summary(table[which(table$fun_group_2 == 6),"PC3_ind"])  # -3.1440 		-1.8850 	-0.9610 	-0.9842  	0.1266  		1.1690 
summary(table[which(table$fun_group_2 == 7),"PC3_ind"])  # -3.0940 		-1.4850 	-0.9063 	-0.6207  	0.3463  		2.1380 


### Now, functional traits values: 
table$trophism <- NA
table$trophism <- ddf[rownames(table), "trophism"]
table$feeding <- NA
table$feeding <- ddf[rownames(table), "feeding"]
table$spawning <- NA
table$spawning <- ddf[rownames(table), "spawning"]
table$size_class <- NA
table$size_class <- ddf[rownames(table), "size_class"]
table$length_max <- NA
table$length_max <- ddf[rownames(table), "length_max"]

### Trophism
# PC1:
fit <- aov(PC1_ind ~ factor(trophism), data= table)
shapiro.test(fit$residuals) # p-value = 7.815e-08 ; NOT normally distributed
fligner.test(PC1_ind ~ factor(trophism),  data= table) # p-value = 0.3436 ; Variances are homogeneous
kruskal.test(PC1_ind ~ factor(trophism), data= table)
# p-value = 0.02258 --> Significant variations of trophism on PC1		
																		 #   Min.		Median		 Mean		Max
summary(table[which(table$trophism == "Carnivore"),"PC1_ind"])			 # -2.3550    	1.0460  	0.7630  	2.5620 
summary(table[which(table$trophism == "Omnivore"),"PC1_ind"])			 # -7.1320 	 	-0.5122 	-0.9504    	1.1830 
summary(table[which(table$trophism == "Omnivore-Carnivore"),"PC1_ind"])	 # 0.5818  	  	1.0500  	1.2950    	2.2540 
summary(table[which(table$trophism == "Omnivore-Herbivore"),"PC1_ind"])	 # -7.1780 	  	0.4833 		-0.1950    	2.7120 
summary(table[which(table$trophism == "Omnivore-Detritivore"),"PC1_ind"])# -5.46100  	0.70030 	-0.03836   	1.85100 
### Like before: Carnivores and Omni-Carni have warmer/ saltier niches
# pairwise.wilcox.test(table$PC1_ind, factor(table$trophism), p.adjust.method = "BH" )


# PC2: 
fligner.test(PC2_ind ~ factor(trophism),  data= table) # p-value = 0.2769 ; Variances are homogeneous
kruskal.test(PC2_ind ~ factor(trophism), data= table)
# p-value = 0.05281 ; NOPE.

# PC3: 
fligner.test(PC3_ind ~ factor(trophism),  data= table) # p-value = 0.06105; Variances are homogeneous
kruskal.test(PC3_ind ~ factor(trophism), data= table)
# p-value = 0.182 ; NOPE.


### Feeding strategy
# PC1:
fligner.test(PC1_ind ~ factor(feeding),  data= table)
# p-value = 0.5753 ; homogeneity of variances
kruskal.test(PC1_ind ~ factor(feeding), data= table)
# p-value = 0.4479 ; non significant

# PC2:
fligner.test(PC2_ind ~ factor(feeding),  data= table)
# p-value = 0.4502 ; homogeneity of variances
kruskal.test(PC2_ind ~ factor(feeding), data= table)
# p-value = 0.9234 ; non significant

# PC3:
fligner.test(PC3_ind ~ factor(feeding),  data= table)
# p-value = 0.132 ; homogeneity of variances
kruskal.test(PC3_ind ~ factor(feeding), data= table)
# p-value = 0.5105 ; non significant


# Spawning strategy...you know it's useless

# Size class
fligner.test(PC1_ind ~ factor(size_class),  data= table)
# p-value = 0.2301 ; homogeneity of variances
kruskal.test(PC1_ind ~ factor(size_class), data= table)
# p-value = 0.01975 ;  significant !..					#   Min.					Median		 Mean					Max
summary(table[which(table$size_class == 1),"PC1_ind"])	# -7.1780 		-1.2980 	-0.1284 	-0.4243  	1.1440  	2.5620 
summary(table[which(table$size_class == 2),"PC1_ind"])	# -4.7720 		-0.5122 	1.0400  	0.4692  	1.4690  	2.7120 
summary(table[which(table$size_class == 3),"PC1_ind"])	# -3.6190 		-1.5520 	-0.1812 	-0.5848  	0.5704  	2.2540
summary(table[which(table$size_class == 4),"PC1_ind"])	# -7.13200  	0.03294  	0.91980  	0.28080  	1.53400  	1.98600 
### But doesn't really make sense does it ?

# PC2:
fligner.test(PC2_ind ~ factor(size_class),  data= table)
# p-value = 0.4551 ; homogeneity of variances
kruskal.test(PC2_ind ~ factor(size_class), data= table)
# p-value = 0.2091 ; non significant

# PC3:
fligner.test(PC3_ind ~ factor(size_class),  data= table)
# p-value = 0.09341 ; homogeneity of variances
kruskal.test(PC3_ind ~ factor(size_class), data= table)
# p-value = 3.738e-05 ; quite significant




### Variation of length_max across PC1/2
cor(table$PC1_ind, table$length_max, method= "spearman") # -0.1339733
cor(table$PC2_ind, table$length_max, method= "spearman") # -0.1951245
cor(table$PC3_ind, table$length_max, method= "spearman") # 0.4106379

summary( lm(length_max ~ PC1_ind, data = table) ) # p-value: 0.415
summary( lm(length_max ~ PC2_ind, data = table) ) # p-value: 0.3145
summary( lm(length_max ~ PC3_ind, data = table) ) # p-value: 0.002197

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

plot <- ggplot(table) + geom_vline(xintercept = 0, linetype = "dashed") +
		geom_point(aes(x=PC1_ind, y=length_max, fill= PC1_ind), shape= 21, size= 4, colour= "black") + 
		scale_fill_gradient2(guide = "colorbar", name = "Coordinate on PC1\n(35.68%)", low="blue", high = "red") + 
		xlab("Coordinate on PC1 (35.68%)") + ylab("Maximal body length (mm)") + theme_bw()	
quartz()
plot
ggsave("plot_size_classes_PC1.pdf", plot= plot, height= 6, width= 7)

		
plot <- ggplot(table) + geom_vline(xintercept = 0, linetype = "dashed") +
		geom_point(aes(x=PC2_ind, y=length_max, fill= PC2_ind), shape= 21, size= 4, colour= "black") + 
		scale_fill_gradient2(guide = "colorbar", name = "Coordinate on PC2\n(21.74%)", low="blue", high = "red") + 
		xlab("Coordinate on PC2 (21.74%)") + ylab("Maximal body length (mm)") + theme_bw()	
quartz()
plot
ggsave("plot_size_classes_PC2.pdf", plot= plot, height= 6, width= 7)


plot <- ggplot(table) + geom_vline(xintercept = 0, linetype = "dashed") +
		geom_point(aes(x=PC3_ind, y=length_max, fill= PC3_ind), shape= 21, size= 4, colour= "black") + 
		scale_fill_gradient2(guide = "colorbar", name = "Coordinate on PC3\n(14.85%)", low="blue", high = "red") + 
		xlab("Coordinate on PC3 (14.85%)") + ylab("Maximal body length (mm)") + theme_bw()	
quartz()
plot
ggsave("plot_size_classes_PC3.pdf", plot= plot, height= 6, width= 7)
### Slight correlation between PC3 (higher MLPAR niches) and size...could be due to predators ? NOT the case.

# Build a new dataframe combining the first two
ns <- data.frame(species = ddf[,"sp.name"],
 				FG = ddf[,"FG"], 
				PC1 = table$PC1_ind, 
				PC2 = table$PC2_ind, 
				PC3 = table$PC3_ind)
				
rownames(ns) <- table$names

# Gravity centers of each size class:
bary <- data.frame(FG = unique(ddf$FG), PC1 = NA, PC2 = NA, PC3 = NA)
#
bary[which(bary$FG == 1), "PC1"] <- mean(ns[which(ns$FG == 1), "PC1"] )
bary[which(bary$FG == 1), "PC2"] <- mean(ns[which(ns$FG == 1), "PC2"] )
bary[which(bary$FG == 1), "PC3"] <- mean(ns[which(ns$FG == 1), "PC3"] )
#
bary[which(bary$FG == 2), "PC1"] <- mean(ns[which(ns$FG == 2), "PC1"] )
bary[which(bary$FG == 2), "PC2"] <- mean(ns[which(ns$FG == 2), "PC2"] )
bary[which(bary$FG == 2), "PC3"] <- mean(ns[which(ns$FG == 2), "PC3"] )
#
bary[which(bary$FG == 3), "PC1"] <- mean(ns[which(ns$FG == 3), "PC1"] )
bary[which(bary$FG == 3), "PC2"] <- mean(ns[which(ns$FG == 3), "PC2"] )
bary[which(bary$FG == 3), "PC3"] <- mean(ns[which(ns$FG == 3), "PC3"] )
#
bary[which(bary$FG == 4), "PC1"] <- mean(ns[which(ns$FG == 4), "PC1"] )
bary[which(bary$FG == 4), "PC2"] <- mean(ns[which(ns$FG == 4), "PC2"] )
bary[which(bary$FG == 4), "PC3"] <- mean(ns[which(ns$FG == 4), "PC3"] )
#
bary[which(bary$FG == 5), "PC1"] <- mean(ns[which(ns$FG == 5), "PC1"] )
bary[which(bary$FG == 5), "PC2"] <- mean(ns[which(ns$FG == 5), "PC2"] )
bary[which(bary$FG == 5), "PC3"] <- mean(ns[which(ns$FG == 5), "PC3"] )
#
bary[which(bary$FG == 6), "PC1"] <- mean(ns[which(ns$FG == 6), "PC1"] )
bary[which(bary$FG == 6), "PC2"] <- mean(ns[which(ns$FG == 6), "PC2"] )
bary[which(bary$FG == 6), "PC3"] <- mean(ns[which(ns$FG == 6), "PC3"] )
#
bary[which(bary$FG == 7), "PC1"] <- mean(ns[which(ns$FG == 7),"PC1"] )
bary[which(bary$FG == 7), "PC2"] <- mean(ns[which(ns$FG == 7),"PC2"] )
bary[which(bary$FG == 7), "PC3"] <- mean(ns[which(ns$FG == 7),"PC3"] )

nrow(ns[which(ns$FG == 1),]) # 10
nrow(ns[which(ns$FG == 2),]) # 15
nrow(ns[which(ns$FG == 3),]) # 8
nrow(ns[which(ns$FG == 4),]) # 22
nrow(ns[which(ns$FG == 5),]) # 20
nrow(ns[which(ns$FG == 6),]) # 10
nrow(ns[which(ns$FG == 7),]) # 10


plot1 <- ggplot() + 
		geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
		geom_point(aes(x= PC1, y= PC2, shape= factor(FG), fill= factor(FG)), data = ns, colour = "black", size = 2) + 
		#geom_segment(aes(x=0, y=0, xend=PC1*4, yend=PC2*4, label= vars), arrow = arrow(length= unit(0.1, "cm"), type= "closed"), colour= "black", data= vars) +
		geom_point(aes(x=PC1, y= PC2, shape= factor(FG), fill= factor(FG)), colour= "black", data = bary, size = 5) + 
	    #scale_fill_manual(name = "Functional\ngroups", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")   ) +
		scale_fill_manual(name = "Functional group", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")  ) +
		scale_shape_manual(name= "Functional group", values = c(21,22,23,24,25,21,22)) + 
		#geom_text(aes(x= PC1*4, y= PC2*4, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
		xlab("PC1 (35.68%)") + ylab("PC2 (21.74%)") + theme_bw()
quartz()
plot1
ggsave("plot_pca_FG_1&2.pdf", plot = plot1, dpi = 300)
		
plot2 <- ggplot() + 
		geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
		geom_point(aes(x= PC2, y= PC3, shape= factor(FG), fill= factor(FG)), data = ns, colour = "black", size = 2) + 
		#geom_segment(aes(x=0, y=0, xend=PC2*4, yend=PC3*4, label= vars), arrow = arrow(length= unit(0.1, "cm"), type= "closed"), colour= "black", data= vars) +
		geom_point(aes(x=PC2, y= PC3, shape= factor(FG), fill= factor(FG)), colour= "black", data = bary, size = 5) + 
	   	#scale_fill_manual(name = "Functional\n groups", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")   ) +
		scale_fill_manual(name = "Functional group", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")  ) +
		scale_shape_manual(name= "Functional group", values = c(21,22,23,24,25,21,22)) + 
		#geom_text(aes(x= PC2*4, y= PC3*4, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
		xlab("PC2 (21.74%)") + ylab("PC3 (14.85%)") + theme_bw()

quartz()
plot2
ggsave("plot_pca_FG_2&3.pdf", plot = plot2, dpi = 300)


### Do the same with trophism on PC1 + PC2 !
# Build a new dataframe combining the first two
ns <- data.frame(species = ddf[!is.na(ddf$trophism),"sp.name"],
 				TG = ddf[!is.na(ddf$trophism),"trophism"], 
				PC1 = table$PC1_ind, 
				PC2 = table$PC2_ind, 
				PC3 = table$PC3_ind)
				
rownames(ns) <- table$names

# Gravity centers of each size class:
bary <- data.frame(TG = unique(ddf$trophism), PC1 = NA, PC2 = NA, PC3 = NA)
#
bary[which(bary$TG == "Carnivore"), "PC1"] <- mean(ns[which(ns$TG == "Carnivore"), "PC1"] )
bary[which(bary$TG == "Carnivore"), "PC2"] <- mean(ns[which(ns$TG == "Carnivore"), "PC2"] )
bary[which(bary$TG == "Carnivore"), "PC3"] <- mean(ns[which(ns$TG == "Carnivore"), "PC3"] )
#
bary[which(bary$TG == "Omnivore"), "PC1"] <- mean(ns[which(ns$TG == "Omnivore"), "PC1"] )
bary[which(bary$TG == "Omnivore"), "PC2"] <- mean(ns[which(ns$TG == "Omnivore"), "PC2"] )
bary[which(bary$TG == "Omnivore"), "PC3"] <- mean(ns[which(ns$TG == "Omnivore"), "PC3"] )
#
bary[which(bary$TG == "Omnivore-Herbivore"), "PC1"] <- mean(ns[which(ns$TG == "Omnivore-Herbivore"), "PC1"] )
bary[which(bary$TG == "Omnivore-Herbivore"), "PC2"] <- mean(ns[which(ns$TG == "Omnivore-Herbivore"), "PC2"] )
bary[which(bary$TG == "Omnivore-Herbivore"), "PC3"] <- mean(ns[which(ns$TG == "Omnivore-Herbivore"), "PC3"] )

bary[which(bary$TG == "Omnivore-Detritivore"), "PC1"] <- mean(ns[which(ns$TG == "Omnivore-Detritivore"), "PC1"] )
bary[which(bary$TG == "Omnivore-Detritivore"), "PC2"] <- mean(ns[which(ns$TG == "Omnivore-Detritivore"), "PC2"] )
bary[which(bary$TG == "Omnivore-Detritivore"), "PC3"] <- mean(ns[which(ns$TG == "Omnivore-Detritivore"), "PC3"] )
#
bary[which(bary$TG == "Omnivore-Carnivore"), "PC1"] <- mean(ns[which(ns$TG == "Omnivore-Carnivore"), "PC1"] )
bary[which(bary$TG == "Omnivore-Carnivore"), "PC2"] <- mean(ns[which(ns$TG == "Omnivore-Carnivore"), "PC2"] )
bary[which(bary$TG == "Omnivore-Carnivore"), "PC3"] <- mean(ns[which(ns$TG == "Omnivore-Carnivore"), "PC3"] )
#
#bary[which(bary$TG == 6), "PC1"] <- mean(ns[which(ns$TG == 6), "PC1"] )
#bary[which(bary$TG == 6), "PC2"] <- mean(ns[which(ns$TG == 6), "PC2"] )
#bary[which(bary$TG == 6), "PC3"] <- mean(ns[which(ns$TG == 6), "PC3"] )
#
#bary[which(bary$TG == 7), "PC1"] <- mean(ns[which(ns$TG == 7),"PC1"] )
#bary[which(bary$TG == 7), "PC2"] <- mean(ns[which(ns$TG == 7),"PC2"] )
#bary[which(bary$TG == 7), "PC3"] <- mean(ns[which(ns$TG == 7),"PC3"] )


### plots in niche qpace
plot1 <- ggplot() + 
		geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
		geom_point(aes(x= PC1, y= PC2, shape= factor(TG), fill= factor(TG)), data = ns, colour = "black", size = 2) + 
		#geom_segment(aes(x=0, y=0, xend=PC1*4, yend=PC2*4, label= vars), arrow = arrow(length= unit(0.1, "cm"), type= "closed"), colour= "black", data= vars) +
		geom_point(aes(x=PC1, y= PC2, shape= factor(TG), fill= factor(TG)), colour= "black", data = bary, size = 5) + 
	    #scale_fill_manual(name = "Functional\ngroups", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")   ) +
		scale_fill_manual(name = "Trophism", values = c("#a50026","#f46d43","#d73027","#3288bd","#66c2a5","#66bd63")  ) +
		scale_shape_manual(name= "Trophism", values = c(21,22,23,24,25,21,22)) + 
		#geom_text(aes(x= PC1*4, y= PC2*4, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
		xlab("PC1 (35.78%)") + ylab("PC2 (22.05%)") + theme_bw()
quartz()
plot1
ggsave("plot_pca_TG_1&2.pdf", plot = plot1, dpi = 300)
		
plot2 <- ggplot() + 
		geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
		geom_point(aes(x= PC2, y= PC3, shape= factor(TG), fill= factor(TG)), data = ns, colour = "black", size = 2) + 
		#geom_segment(aes(x=0, y=0, xend=PC2*4, yend=PC3*4, label= vars), arrow = arrow(length= unit(0.1, "cm"), type= "closed"), colour= "black", data= vars) +
		geom_point(aes(x=PC2, y= PC3, shape= factor(TG), fill= factor(TG)), colour= "black", data = bary, size = 5) + 
	   	#scale_fill_manual(name = "Functional\n groups", values = c("#a50026","#d73027","#006837","#66c2a5","#f46d43","#3288bd","#66bd63")   ) +
		scale_fill_manual(name = "Trophism", values = c("#a50026","#f46d43","#d73027","#3288bd","#66c2a5","#66bd63")  ) +
		scale_shape_manual(name= "Trophism", values = c(21,22,23,24,25,21,22)) + 
		#geom_text(aes(x= PC2*4, y= PC3*4, label= vars), hjust=-0.1, vjust=1, size = 2.5, data = vars) + 
		xlab("PC2 (22.05%)") + ylab("PC3 (15.15%)") + theme_bw()

quartz()
plot2
ggsave("plot_pca_TG_2&3.pdf", plot = plot2, dpi = 300)






##### D) Check variation of missing niche traits according to FGs 
### 
fligner.test(center_MLPAR1 ~ factor(spawning),  data= ddf) 
kruskal.test(center_MLPAR1 ~ factor(spawning), data= ddf)  





##### E) Computing correlation between two cophenetic distance matrices of the two trees: functional and niche
library("dendextend")
dend1 <- as.dendrogram(fit)     # dendrogram of niche space
dend2 <- as.dendrogram(fit_mca) # dendrogram of functional space
match_order_by_labels(dend1, dend2)
length(labels(dend1)) ; length(labels(dend2))

cophenetic(dend1)
cophenetic(dend2)

cor_cophenetic(dend1, dend2, "spearman")
cor(cophenetic(fit), cophenetic(fit_mca), method = "spearman") # 0.06912228
cor(cophenetic(dend1), cophenetic(dend2), method = "spearman") # 0.63 !! Gut gut
cor_cophenetic(dendlist(dend1, dend2), method = "spearman")


#cor(cophenetic(dend1), cophenetic(dend2), method = "kendall") # 0.5524801 (for info)
#cor(cophenetic(dend1), cophenetic(dend2), method = "pearson") # 0.5490363 (for info)

### For exact p-value one should result to a permutation test. 
# One such option will be to permute over the labels of one tree many times, and calculating the distribution under the null hypothesis (keeping the trees topologies constant).
# ...

# Randomize tips around dend1 n times (n = 999)
n <- 499
library("parallel")
null_distrib <- mclapply(1:n, function(i) {
						# Randomize tips
						rand_dend1 <- dend1
						#rand_dend2 <- dend2
						labels(rand_dend1) <- sample(labels(dend1), size = length(labels(dend1)))
						#labels(rand_dend2) <- sample(labels(dend2), size = length(labels(dend2)))
						# labels(rand_dend1) ; labels(dend1)
						#cor <- cor(cophenetic(rand_dend1), cophenetic(dend2), method = "spearman")
						cor <- cor_cophenetic(rand_dend1, dend2, "spearman")
						return(cor)
}, mc.cores = 3)

length(null_distrib)

null_cor <- do.call(rbind, null_distrib)

#obs_cor <- cor(cophenetic(dend1), cophenetic(dend2), method = "spearman") # 0.63 !! Gut gut
obs_cor <- cor_cophenetic(dend1, dend2, "spearman")
# final p-value should be computed as follows: 
p.value <- sum(null_cor >= obs_cor) / (n+1)
# or standardized effect size: 
ses <- (obs_cor - mean(null_cor)) / sd(null_cor)


### Other stuff maybe ? 
# The all.equal.dendrogram function makes a global comparison of two or more dendrograms trees.
all.equal(dend1, dend1) # TRUE
all.equal(dend1, dend2) # "Difference in branch heights -  Mean relative difference: 0.5739686 "
all.equal(dends_15_51) # Same as above

### Baker’s Gamma Index :
# Baker’s Gamma Index (see baker’s paper from 1974) is a measure of association (similarity) between two trees of Hierarchical clustering (dendrograms). 
# It is defined as the rank correlation between the stages at which pairs of objects combine in each of the two trees.
# Or more detailed: It is calculated by taking two items, and see what is the highest possible level of k 
# (number of cluster groups created when cutting the tree) for which the two item still belongs to the same tree. 
# That k is returned, and the same is done for these two items for the second tree. 
# There are n over 2 combinations of such pairs of items from the items in the tree, and all of these numbers are calculated for each of the two trees. 
# Then, these two sets of numbers (a set for the items in each tree) are paired according to the pairs of items compared, and a Spearman correlation is calculated.
# The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar. 
# For exact p-value one should use a permutation test. 
# One such option will be to permute over the labels of one tree many times, calculating the distribution under the null hypothesis (keeping the trees topologies constant).
# Notice that this measure is not affected by the height of a branch but only of its relative position compared with other branches.

#cor_bakers_gamma(dend1, dend2, use_labels_not_values = F, to_plot = T)
#dend1 <- match_order_by_labels(dend1, dend2) # if you are not sure
#cor_bakers_gamma(dend1, dend2, use_labels_not_values = T)  
### Not working for some reason...




##### F) Finally, check the species/ fucntionl groups and traits for which MLPAR/logChla models are performing best
ddf[which(ddf$avgTSS_logChla >= 0.3),c("trophism","feeding","FG")]

ddf[which(ddf$avgTSS_MLPAR1 >= 0.3),c("trophism","feeding","FG")]

ddf[which(ddf$avgTSS_MLD1 >= 0.3),c("trophism","feeding","FG")]




##### G) 20/05/16 - Assess differences in average univarite TSS values (for part 1 of Results)
# May need to melt first

m <- melt(ddf[,c(5,15,25,35,45,55)])
fit <- aov(value ~ factor(variable), data= m)
shapiro.test(fit$residuals) # p-value = 0.4193 ; normally distributed residuals 
bartlett.test(value ~ factor(variable),  data= m) # p-value = 2.412e-10 ; Variances are not homogeneous...

# T.test ? 
shapiro.test(ddf$avgTSS_logChla) # Normally distributed: MLD, PAR, SST (p-value > 0.05)
fligner.test(ddf$avgTSS_SSS)

wilcox.test(ddf$avgTSS_SSS, ddf$avgTSS_dSST) # p-value = 0.8763
wilcox.test(ddf$avgTSS_SST, ddf$avgTSS_dSST) # p-value = 0.003142
wilcox.test(ddf$avgTSS_SST, ddf$avgTSS_SSS) # p-value = 0.02277

wilcox.test(ddf$avgTSS_SST, ddf$avgTSS_MLD1) # p-value < 2.2e-16
wilcox.test(ddf$avgTSS_SST, ddf$avgTSS_logChla) # p-value < 2.2e-16
wilcox.test(ddf$avgTSS_SST, ddf$avgTSS_MLPAR) # p-value < 2.2e-16

wilcox.test(ddf$avgTSS_dSST, ddf$avgTSS_MLD1) # p-value < 2.2e-16
wilcox.test(ddf$avgTSS_dSST, ddf$avgTSS_logChla) # p-value < 2.2e-16
wilcox.test(ddf$avgTSS_dSST, ddf$avgTSS_MLPAR) # p-value < 2.2e-16

wilcox.test(ddf$avgTSS_SSS, ddf$avgTSS_MLD1) # p-value < 2.2e-16
wilcox.test(ddf$avgTSS_SSS, ddf$avgTSS_logChla) # p-value < 2.2e-16
wilcox.test(ddf$avgTSS_SSS, ddf$avgTSS_MLPAR) # p-value < 2.2e-16


wilcox.test(ddf$avgTSS_MLPAR, ddf$avgTSS_logChla) # p-value = 0.3657
wilcox.test(ddf$avgTSS_MLD1, ddf$avgTSS_logChla) # p-value = 0.0002704
wilcox.test(ddf$avgTSS_MLD1, ddf$avgTSS_MLPAR) # p-value = 3.396e-08
wilcox.test(ddf$avgTSS_SSS, ddf$avgTSS_logChla) # p-value < 2.2e-16

# ∆SST - SSS ; SST ; MLD ; MLPAR - logChla



### Plot species' univarite ncihes ! (Fig. 3 and 4 of Brun et al., 2015): each point being a species, plot centers with breadth
ddf[ddf$center_logChla > 0.6,]
summary(lm(center_MLD1 ~ breadth_MLD1, data = ddf))

# 1) SST
summary(lm(breadth_SST ~ center_SST, data = ddf)) # p-value: 0.3371
cor(ddf$breadth_SST, ddf$center_SST, method = "spearman") # -0.07539613

plot <- ggplot(ddf) + geom_point(aes(x= center_SST, y= breadth_SST), shape = 21, colour = "black", fill= "red") + 
	xlab("SST center (°C)") + ylab("SST breadth (°C)") + theme_linedraw()

quartz()
plot
ggsave("plot_SSTniches.pdf", plot= plot, dpi = 300)

# 2) ∆SST
summary(lm(breadth_dSST ~ center_dSST, data = ddf)) # p-value: 2.087e-14
cor(ddf$breadth_dSST, ddf$center_dSST, method = "spearman") # 0.6364058
lm <- lm(breadth_dSST ~ center_dSST, data = ddf)
plot <- ggplot(ddf) + geom_point(aes(x= center_dSST, y= breadth_dSST), shape = 21, colour = "black", fill= "red") + 
	geom_path(aes(x= center_dSST, y= (1.449536*center_dSST) - 5.605651), colour = "black", size = 1) +
	geom_ribbon(aes(x= center_dSST, ymin= predict(lm, interval="confidence")[,2], ymax= predict(lm, interval="confidence")[,3]), alpha= 0.2) + 
	xlab("∆SST center (°C)") + ylab("∆SST breadth (°C)") + theme_linedraw()

quartz()
plot
ggsave("plot_dSSTniches.pdf", plot= plot, dpi = 300)

# 3) SSS
summary(lm(breadth_SSS ~ center_SSS, data = ddf)) # p-value: 4.39e-16
cor(ddf$breadth_SSS, ddf$center_SSS, method = "spearman") # -0.8236038
lm <- lm(breadth_SSS ~ center_SSS, data = ddf)
plot <- ggplot(ddf) + geom_point(aes(x= center_SSS, y= breadth_SSS), shape = 21, colour = "black", fill= "red") + 
	geom_path(aes(x= center_SSS, y= (-1.235774*center_SSS) + 51.260820), colour = "black", size = 1) +
	geom_ribbon(aes(x= center_SSS, ymin= predict(lm, interval="confidence")[,2], ymax= predict(lm, interval="confidence")[,3]), alpha= 0.2) + 
	xlab("SSS center") + ylab("SSS breadth") + theme_linedraw()

quartz()
plot
ggsave("plot_SSSniches.pdf", plot= plot, dpi = 300)

# 4) MLD
summary(lm(breadth_MLD1 ~ center_MLD1, data = ddf)) # p-value < 2.2e-16
cor(ddf$breadth_MLD1, ddf$center_MLD1, method = "spearman") # 0.9079628
lm <- lm(breadth_MLD1 ~ center_MLD1, data = ddf)
coefficients(lm)
plot <- ggplot(ddf) + geom_point(aes(x= center_MLD1, y= breadth_MLD1), shape = 21, colour = "black", fill= "red") + 
	geom_path(aes(x= center_MLD1, y= (1.805067*center_MLD1) + 29.434573), colour = "black", size = 1) +
	geom_ribbon(aes(x= center_MLD1, ymin= predict(lm, interval="confidence")[,2], ymax= predict(lm, interval="confidence")[,3]), alpha= 0.2) + 
	xlab("MLD center (m)") + ylab("MLD breadth (m)") + theme_linedraw()

quartz()
plot
ggsave("plot_MLDniches.pdf", plot= plot, dpi = 300)

# 5) MLPAR
summary(lm(breadth_MLPAR1 ~ center_MLPAR1, data = ddf)) # p-value = 0.04547 
cor(ddf$breadth_MLPAR1, ddf$center_MLPAR1, method = "spearman") # -0.01459056
lm <- lm(breadth_MLPAR1 ~ center_MLPAR1, data = ddf)
coefficients(lm)
plot <- ggplot(ddf) + geom_point(aes(x= center_MLPAR1, y= breadth_MLPAR1), shape = 21, colour = "black", fill= "red") + 
	geom_path(aes(x= center_MLPAR1, y= (0.1341208*center_MLPAR1) + 31.3353295), colour = "black", size = 1) +
	geom_ribbon(aes(x= center_MLPAR1, ymin= predict(lm, interval="confidence")[,2], ymax= predict(lm, interval="confidence")[,3]), alpha= 0.2) + 
	xlab("MLPAR center (E/m2.d)") + ylab("MLPAR breadth (E/m2.d)") + theme_linedraw()

quartz()
plot
ggsave("plot_MLPARniches.pdf", plot= plot, dpi = 300)

# 6) logChla 
summary(lm(breadth_logChla ~ center_logChla, data = ddf)) # p-value = 0.06852
cor(ddf$breadth_logChla, ddf$center_logChla, method = "spearman") # 0.1931682
lm <- lm(breadth_logChla ~ center_logChla, data = ddf)
coefficients(lm)
plot <- ggplot(ddf) + geom_point(aes(x= center_logChla, y= breadth_logChla), shape = 21, colour = "black", fill= "red") + 
	xlab("logChla center (mg/m3)") + ylab("logChla breadth (mg/m3)") + theme_linedraw()

quartz()
plot
ggsave("plot_Chlaniches.pdf", plot= plot, dpi = 300)



##### 01/06/2016: Re-checking significance tests for univariate TSS mean values between FG (they were all actually OK)

fligner.test(avgTSS_SST ~ factor(FG),  data= ddf) # p-value = 0.1741 ; Variances are homogeneous
kruskal.test(avgTSS_SST ~ factor(FG), data= ddf) # p-value = 0.03212
pairwise.wilcox.test(ddf$avgTSS_SST, factor(ddf$FG), p.adjust.method = "BH")
#  1     2     3     4     5     6    
#2 0.336 -     -     -     -     -    
#3 0.739 0.523 -     -     -     -    
#4 0.961 0.421 0.847 -     -     -    
#5 0.858 0.523 0.961 0.872 -     -    
#6 0.160 0.044 0.132 0.132 0.132 -    
#7 0.523 0.101 0.160 0.421 0.315 0.354
### FG 2 better modelled than FG 6 with SST

fligner.test(avgTSS_dSST ~ factor(FG),  data= ddf) # p-value = 0.04209 ; Variances are homogeneous
### Cannot test


fligner.test(avgTSS_SSS ~ factor(FG),  data= ddf) # p-value = 0.1852 ; Variances are homogeneous
kruskal.test(avgTSS_SSS ~ factor(FG), data= ddf) # p-value = 0.003431
pairwise.wilcox.test(ddf$avgTSS_SSS, factor(ddf$FG), p.adjust.method = "BH")
#   1     2     3     4     5     6    
#2 0.892 -     -     -     -     -    
#3 0.031 0.056 -     -     -     -    
#4 0.277 0.310 0.457 -     -     -    
#5 0.690 0.890 0.031 0.121 -     -    
#6 0.031 0.056 0.690 0.199 0.019 -    
#7 0.199 0.199 0.598 0.508 0.077 0.690
### FG 6 > FG 5 
### FG 6 > FG 1
### FG 5 > FG 3

fligner.test(avgTSS_MLD1 ~ factor(FG),  data= ddf) # p-value = 0.04701 
### Cannot test

fligner.test(avgTSS_MLPAR1 ~ factor(FG),  data= ddf) # p-value = 0.2497 ; Variances are homogeneous
kruskal.test(avgTSS_MLPAR1 ~ factor(FG), data= ddf) # p-value = 0.02099
pairwise.wilcox.test(ddf$avgTSS_MLPAR1, factor(ddf$FG), p.adjust.method = "BH")
#   1     2     3     4     5     6    
#2 0.891 -     -     -     -     -    
#3 0.103 0.103 -     -     -     -    
#4 0.686 0.612 0.103 -     -     -    
#5 0.917 0.917 0.031 0.418 -     -    
#6 0.652 0.301 0.031 0.119 0.333 -    
#7 0.624 0.612 0.103 0.805 0.242 0.103
### 

fligner.test(avgTSS_logChla ~ factor(FG),  data= ddf) # p-value = 0.5877
kruskal.test(avgTSS_logChla ~ factor(FG), data= ddf) # p-value = 0.4765
### No diff


fligner.test(avgTSS_MLD1 ~ factor(spawning),  data= ddf) # p-value = 0.3913 ; Variances are homogeneous
kruskal.test(avgTSS_MLD1 ~ factor(spawning), data= ddf) # p-value = 0.01168
pairwise.wilcox.test(ddf$avgTSS_MLD1, factor(ddf$spawning), p.adjust.method = "BH")

summary(ddf[which(ddf$spawning == "Broadcaster"),"avgTSS_MLD1"])
summary(ddf[which(ddf$spawning == "Sac-spawner"),"avgTSS_MLD1"])

summary(ddf[which(ddf$spawning == "Broadcaster"),"avgTSS_MLPAR1"])
summary(ddf[which(ddf$spawning == "Sac-spawner"),"avgTSS_MLPAR1"])


##### 03/03/2017: Work on Fig.2 (panel showing variations of niche characteristics among species) before submission to L&O
head(ddf)

# SST
summary(lm(breadth_SST ~ center_SST, data = ddf)) # p-value: 0.3371
#quartz()
plot <- ggplot() + geom_point(aes(x = center_SST, y = breadth_SST), data = ddf, shape = 21, colour = "black", fill = "#d53e4f", size = 3) + 
			xlab("SST niche center (°C)") + ylab("SST niche breadth (°C)") + theme_bw()
ggsave(filename = "niche_traits_sst.pdf", plot = plot, height = 5, width = 5, dpi = 300)

# SSS
summary(lm(breadth_SSS ~ center_SSS, data = ddf)) # p-value: 4.39e-16
lm <- lm(breadth_SSS ~ center_SSS, data = ddf)
lm_ddf <- data.frame(predict(lm, interval="confidence"))
lm_ddf$raw <- ddf$center_SSS
#quartz()
plot <- ggplot() + geom_point(aes(x = center_SSS, y = breadth_SSS), data = ddf, shape = 21, colour = "black", fill = "#f46d43", size = 3) + 
			geom_path(aes(x = raw, y = fit), data = lm_ddf) + geom_ribbon(aes(x = raw, ymax = upr, ymin = lwr), data = lm_ddf, alpha = 0.5, fill = "grey70") + 
			xlab("SSS niche center") + ylab("SSS niche breadth") + theme_bw()
ggsave(filename = "niche_traits_sss.pdf", plot = plot, height = 5, width = 5, dpi = 300)

# dSST
summary(lm(breadth_dSST ~ center_dSST, data = ddf)) # p-value: 2.087e-14
lm <- lm(breadth_dSST ~ center_dSST, data = ddf)
lm_ddf <- data.frame(predict(lm, interval="confidence"))
lm_ddf$raw <- ddf$center_dSST
#quartz()
plot <- ggplot() + geom_point(aes(x = center_dSST, y = breadth_dSST), data = ddf, shape = 21, colour = "black", fill = "#fdae61", size = 3) + 
			geom_path(aes(x = raw, y = fit), data = lm_ddf) + geom_ribbon(aes(x = raw, ymax = upr, ymin = lwr), data = lm_ddf, alpha = 0.5, fill = "grey70") + 
			xlab("∆SST niche center") + ylab("∆SST niche breadth") + theme_bw()
ggsave(filename = "niche_traits_dsst.pdf", plot = plot, height = 5, width = 5, dpi = 300)

# MLD
summary(lm(breadth_MLD1 ~ center_MLD1, data = ddf)) # p-value: < 2.2e-16
lm <- lm(breadth_MLD1 ~ center_MLD1, data = ddf)
lm_ddf <- data.frame(predict(lm, interval="confidence"))
lm_ddf$raw <- ddf$center_MLD1
#quartz()
plot <- ggplot() + geom_point(aes(x = center_MLD1, y = breadth_MLD1), data = ddf, shape = 21, colour = "black", fill = "#3288bd", size = 3) + 
			geom_path(aes(x = raw, y = fit), data = lm_ddf) + geom_ribbon(aes(x = raw, ymax = upr, ymin = lwr), data = lm_ddf, alpha = 0.5, fill = "grey70") + 
			xlab("MLD niche center (m)") + ylab("MLD niche breadth (m)") + theme_bw()
ggsave(filename = "niche_traits_mld.pdf", plot = plot, height = 5, width = 5, dpi = 300)

# MLPAR
summary(lm(breadth_MLPAR1 ~ center_MLPAR1, data = ddf)) # p-value: 0.04547
lm <- lm(breadth_MLPAR1 ~ center_MLPAR1, data = ddf)
lm_ddf <- data.frame(predict(lm, interval="confidence"))
lm_ddf$raw <- ddf$center_MLPAR1
#quartz()
plot <- ggplot() + geom_point(aes(x = center_MLPAR1, y = breadth_MLPAR1), data = ddf, shape = 21, colour = "black", fill = "#e6f598", size = 3) + 
			#geom_path(aes(x = raw, y = fit), data = lm_ddf) + geom_ribbon(aes(x = raw, ymax = upr, ymin = lwr), data = lm_ddf, alpha = 0.5, fill = "grey70") + 
			xlab("MLPAR niche center (E/m2.d)") + ylab("MLPAR niche breadth (E/m2.d)") + theme_bw()
			
ggsave(filename = "niche_traits_mlpar.pdf", plot = plot, height = 5, width = 5, dpi = 300)

# logChla
summary(lm(breadth_logChla ~ center_logChla, data = ddf)) # p-value: 0.06852
#lm <- lm(breadth_logChla ~ center_logChla, data = ddf)
#lm_ddf <- data.frame(predict(lm, interval="confidence"))
#lm_ddf$raw <- ddf$center_logChla
#quartz()
plot <- ggplot() + geom_point(aes(x = center_logChla, y = breadth_logChla), data = ddf, shape = 21, colour = "black", fill = "#66c2a5", size = 3) + 
			#geom_path(aes(x = raw, y = fit), data = lm_ddf) + geom_ribbon(aes(x = raw, ymax = upr, ymin = lwr), data = lm_ddf, alpha = 0.5, fill = "grey70") + 
			xlab("logChla niche center (mg/m3)") + ylab("logChla niche breadth (mg/m3)") + theme_bw()
			
ggsave(filename = "niche_traits_logChla.pdf", plot = plot, height = 5, width = 5, dpi = 300)			
			
			
			
