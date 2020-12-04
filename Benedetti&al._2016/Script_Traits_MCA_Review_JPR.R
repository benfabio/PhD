# Tests sur la diversité des communautés de copépodes entre Biorégions
# Sakina, 30 April 2015
# Script R
#
#              FUNCTIONAL GROUPS OF MEDITERRANEAN COPEPODS
#
# Proteus4 : List of 187 copepod species in the Med Sea
#
# Estimating the number of functional groups of Mediterranean Copepod -> MCA -> 3 groups!
#
# Last update : 2015 May 06 (new trait table -> Proteus5 and 5 traits)
#               2015 May 07 (new trait table with less NA in reprod)
#               2015 May 18 (new trait table from 2015 May 12)
#               2015 May 26 (without Siphono*)
#               2015 June 03 (correction of trophic regime)
#               2015 August 03: new trait table for JPR review
#               2015 August 04: choice of analysis with Fabio for the revision
#
####################################################################################


# ------------------------------------------------------------------------------
setwd("~/Documents/Travail_LOV/Recherche_LOV/Med-Perseus-Mermex/Fabio/Bioregions_copepod_div/SCRIPTS_divers")
library(reshape2)
library(ggplot2)
library(dplyr)
library(plyr)
library(dendextend)
#library(MASS)
library(FactoMineR)
# ------------------------------------------------------------------------------
# install.packages(c("ggplot2","plyr","FactoMineR","BAT","FD","ade4","heatmap3",
# "ape","geiger","phytools"))


####################################
#         TRAIT DATA TABLE
####################################

# Traits description : 12 variables incl. (5 traits) for the proteus5 species list,
# whitout Siphonostomatoida (two semi-parasitic species)
# 193 species

traits<-read.csv("../DATA_csv_xls_tre/Proteus_species_traits_V2.csv",sep=";")
# temp <- read.csv("../DATA_csv_xls_tre/Proteus5_Species_traits_2015_06_03_191sp.csv",sep=";")
# traits[!(traits$sp.name%in%temp$sp.name),]
# #                     sp.name             order     superfamily          family        genus length_min length_max       trophism feeding spawning  DVM layer
# # 155 Pontoeciella abyssicola Siphonostomatoida Pontoeciellidae Pontoeciellidae Pontoeciella        0.7       1.65 semi-parasitic    <NA>     <NA>   No  <NA>
# # 157           Ratania flava Siphonostomatoida      Rataniidae      Rataniidae      Ratania        1.0       1.40 semi-parasitic    <NA>     <NA> <NA>  <NA>
# filter(traits,trophism=="semi-parasitic")
traits <- filter(traits,order!="Siphonostomatoida")

# 191 species (without semi-parasitic)

colnames(traits)
# [1] "sp.name"     "order"       "superfamily" "family"     
# [5] "genus"       "length_min"  "length_max"  "trophism"   
# [9] "feeding"     "spawning"    "DVM"         "layer" 
head(traits)
summary(traits[,-1])
# order               superfamily             family              genus    
# Calanoida        :122   Clausocalanoidea:36   Sapphirinidae  : 18   Calocalanus  : 12  
# Cyclopoida       : 61   Diaptomoidea    :31   Paracalanidae  : 16   Oithona      : 12  
# Harpacticoida    :  6   Calanoidea      :22   Oncaeidae      : 13   Sapphirina   : 12  
# Mormonilloida    :  2   Sapphirinidae   :18   Corycaeus      : 12   Candacia     : 10  
# Siphonostomatoida:  0   Arietelloidea   :14   Oithonidae     : 12   Clausocalanus:  8  
# Corycaeidae     :14   Clausocalanidae: 11   Centropages  :  7  
# (Other)         :56   (Other)        :109   (Other)      :130  
#
# length_min      length_max                     trophism           feeding          spawning 
# Min.   :0.300   Min.   : 0.500   Carnivore           :54   Active Ambush:28   Broadcaster:52  
# 1st Qu.:0.680   1st Qu.: 1.300   Omnivore            :33   Cruise       :34   Sac-spawner:84  
# Median :1.100   Median : 2.000   Omnivore-Carnivore  : 8   Filter       :55   NA's       :55  
# Mean   :1.328   Mean   : 2.517   Omnivore-Detritivore:28   Mixed        :13                   
# 3rd Qu.:1.730   3rd Qu.: 3.175   Omnivore-Herbivore  :51   NA's         :61                   
# Max.   :5.450   Max.   :11.000   NA's                :17                                       
#
#       DVM                  layer   
#  No     :57   Epipelagic      :77  
#  Reverse: 5   Epimesopelagic  :63  
#  Strong :13   Epibathypelagic :33  
#  Weak   :50   Mesobathypelagic:11  
#  NA's   :66   Mesopelagic     : 3  
#               (Other)         : 2  
#               NA's            : 2 

length(unique(traits$order))          # 4 orders
length(unique(traits$superfamily))    # 19 superfamilies
length(unique(traits$family))         # 34 families
length(unique(traits$genus))          # 74 genus



#############################################
#    NOT REMOVING SPECIES WITH too many NAs 
############################################

# ####################################
# #    REMOVING SPECIES WITH NAs ?
# ####################################
# 
# # Number of NAs in traits
traits$na <- is.na(traits$trophism)+is.na(traits$feeding)+ is.na(traits$spawning) + is.na(traits$DVM) +is.na(traits$layer)
traits$na3 <- is.na(traits$trophism)+is.na(traits$feeding)+ is.na(traits$spawning)#+is.na(traits$layer)

# traits_red <- filter(traits,na3==0) # 99 species with 4 fully described traits for MCA
# 
# #traits <- 
# filter(traits,nad<=1 ) # 154 species
# filter(traits,nad<1 )    # 98 species

# filter(traits_,is.na(trophism))$sp.name # 17 species with unknow trophism
# 
# filter(traits_,na>2)$sp.name # 22 species with 3 or 4 NAs in traits
# # [1] Archescolecithrix auropecten  Augaptilus spinifrons         Candacia giesbrechti         
# # [4] Copilia mediterranea          Copilia vitrea                Diaixis pygmaea              
# # [7] Disco minutus                 Distioculus minor             Homeognathia brevis          
# # [10] Monothula subtilis            Mormonilla phasma             Neomormonilla minor          
# # [13] Pachos punctatum              Paracartia latisetosa         Parapontella brevicornis     
# # [16] Phaenna spinifera             Pontoeciella abyssicola       Ratania flava                
# # [19] Vettoria granulosa            Vettoria longifurca           Vettoria parva               
# # [22] Xanthocalanus agilis 
# 
# traits <- filter(traits_,na<3) # 171 species
# filter(traits,is.na(trophism))$sp.name # 3 species remaining with unknown trophism
# 
# traits <- filter(traits_,na<3 & !is.na(trophism)) # 168 species
# 
# # 168 species !  
# summary(traits[,-1])
# # order               superfamily             family             genus    
# # Calanoida        :113   Clausocalanoidea:32   Paracalanidae  :16   Calocalanus  : 12  
# # Cyclopoida       : 50   Diaptomoidea    :28   Corycaeus      :12   Oithona      : 12  
# # Harpacticoida    :  5   Calanoidea      :22   Oithonidae     :12   Sapphirina   : 12  
# # Mormonilloida    :  0   Corycaeidae     :14   Oncaeidae      :12   Candacia     :  9  
# # Siphonostomatoida:  0   Arietelloidea   :12   Sapphirinidae  :12   Clausocalanus:  8  
# #  Oithonidae      :12     Clausocalanidae:11     Centropages  :  7  
# #  (Other)         :48     (Other)        :93     (Other)      :108  
# #
# # length_min      length_max                     trophism           feeding          spawning 
# #  Min.   :0.300   Min.   : 0.500   Carnivore           :51   Active Ambush:28   Broadcaster:51  
# #  1st Qu.:0.665   1st Qu.: 1.308   Omnivore            :32   Cruise       :34   Sac-spawner:72  
# #  Median :1.100   Median : 1.990   Omnivore-Carnivore  : 8   Filter       :55   NA's       :45  
# #  Mean   :1.328   Mean   : 2.516   Omnivore-Detritivore:27   Mixed        :13                   
# #  3rd Qu.:1.715   3rd Qu.: 3.300   Omnivore-Herbivore  :50   NA's         :38                   
# #  Max.   :5.450   Max.   :11.000   semi-parasitic      : 0                                      
# # 
# # DVM                  layer   
# # No     :53   Epipelagic      :72  
# # Reverse: 4   Epimesopelagic  :53  
# # Strong :13   Epibathypelagic :32  
# # Weak   :47   Mesobathypelagic: 7  
# # NA's   :51   Mesopelagic     : 3  
# #              (Other)         : 0  
# #              NA's            : 1  
# 
# length(unique(traits$order))          # 3 orders (not 5)
# length(unique(traits$superfamily))    # 16 superfamilies (not 21)
# length(unique(traits$family))         # 27 families not (36)
# length(unique(traits$genus))          # 58 genus (not 76)
# 
# 
# # ORDER SUPER FAMILLY BY ORDER
# 
# traits$superfamily2 <-   factor(traits$superfamily, levels = c(
#   #"Mormonillidae","Pontoeciellidae","Rataniidae",
#   "Thalestroidea","Euterpinidae","Clytemnestridae","Ectinosomatidae",
#   "Sapphirinidae","Corycaeidae", #"Lubbockiidae",
#   "Oncaeidae","Oithonidae",#"incertaesedis(Cyclopoida)",
#   "Diaptomoidea","Clausocalanoidea","Arietelloidea","Calanoidea","Eucalanoidea",
#   "Lucicutiidae","Spinocalanoidea","Metridinidae"))
# traits$order2 <-   factor(traits$order, levels = c( #"Mormonilloida","Siphonostomatoida",
#   "Harpacticoida","Cyclopoida","Calanoida"))


# 7 Functional Traits -> length_min length_max trophism feeding  DVM  spawning layer
#-----------------------------------------------------------------------------------------------

#traits <- traits_red
data <-traits[,c(6,7,9,10,11)]

row.names(data) <- traits$sp.name

datab <-data
dim(datab)
# 191 7

#########################################
#    MCA on species qualitative traits
#########################################

# Defining size classes based on dendro on length -> 4 classes
#---------------------------------------------------------------

l_distance <- dist(traits$length_max, method = "euclidean") # distance matrix
hcl <- hclust(l_distance, method= "ward.D") 
plot(hcl, hang= -1, labels =row.names(data)) # display dendogram with species names
nn <- dim(datab)
par(cex=0.6); plot(hcl,labels = c(1:nn[1]),hang= -1) # display dendogram with species numbers
s_groups <- cutree(hcl, k = 4) # cut tree into 4 clusters
rect.hclust(hcl, 4)
datab$size_class <- as.factor(s_groups)
levels(datab$size_class )<- c("Size_1","Size_2","Size_3","Size_4")#,"Size_5")

# Checking the size ranges
c(min(filter(datab,size_class=="Size_1")[,2]),max(filter(datab,size_class=="Size_1")[,2]))  # [1] 0.5 1.8
c(min(filter(datab,size_class=="Size_1")[,2]),max(filter(datab,size_class=="Size_1")[,2]))  # [1] 0.5 1.8
c(min(filter(datab,size_class=="Size_2")[,2]),max(filter(datab,size_class=="Size_2")[,2]))  # [1] 1.89 2.85
c(min(filter(datab,size_class=="Size_3")[,2]),max(filter(datab,size_class=="Size_3")[,2]))  # [1] 3.0 5.7
c(min(filter(datab,size_class=="Size_4")[,2]),max(filter(datab,size_class=="Size_4")[,2]))  # [1]  6.1 11.0

# Looking at the size ditrib per class
templ <- data.frame(lm = datab$length_max,lc=datab$size_class)
ggplot(templ,aes(x=lm,fill=lc))+geom_bar(stat="bin",binwidth = 0.1)
ggplot(templ,aes(x=lm,fill=lc))+geom_bar(stat="bin",binwidth = 0.3)


# # Convert size into 5 quantiles
# #----------------------------------
# 
# # lminq <- quantile(datab$length_min,seq(0,1,.2))
# # lmin <- cut(datab$length_min,unique(lminq),include.lowest=TRUE)
# lmaxq <- quantile(datab$length_max,seq(0,1,.2))
# lmax <- cut(datab$length_max,unique(lmaxq),include.lowest=TRUE)
# #dataf<-cbind(datab[,-c(1,2)],lmin,lmax)
# #head(dataf)
# dataf<-cbind(datab[,-c(1,2)],lmax)
# head(dataf)



# Convert layers into new classes
#----------------------------------
summary(traits$layer)
datab$epi <- (traits$layer=="Epibathypelagic")|(traits$layer=="Epimesopelagic")|(traits$layer=="Epipelagic")
datab$meso <- (traits$layer=="Epibathypelagic")|(traits$layer=="Epimesopelagic")|(traits$layer=="Mesobathypelagic")|(traits$layer=="Mesopelagic")
datab$bathy <-(traits$layer=="Bathypelagic")|(traits$layer=="Epibathypelagic")|(traits$layer=="Hyperbenthic")|(traits$layer=="Mesobathypelagic")


# Convert trophism into new classes
#----------------------------------
summary(traits$trophism)
datab$C <- (traits$trophism=="Carnivore")|(traits$trophism=="Omnivore-Carnivore")
datab$O <- (traits$trophism=="Omnivore")|(traits$trophism=="Omnivore-Carnivore")|(traits$trophism=="Omnivore-Detritivore")|(traits$trophism=="Omnivore-Herbivore")
datab$D <-(traits$trophism=="Omnivore-Detritivore")
datab$H <-(traits$trophism=="Omnivore-Herbivore")



# Rename levels
#----------------------------------
#levels(datab$trophism) <- list("C"="Carnivore", "O"="Omnivore", "O-C"="Omnivore-Carnivore",
#                               "O-D"="Omnivore-Detritivore","O-H"="Omnivore-Herbivore")
#levels(datab$feeding)<-list("Af"="Active Ambush","Cf"="Cruise","Ff"="Filter","Mf"="Mixed")
#levels(datab$DVM) <-   list("no_DVM"="No", "weak_DVM"="Weak", "strong_DVM"="Strong", "reverse_DVM"="Reverse")
#levels(datab$lmin) <- c("s_1","s_2","s_3","s_4")#,"s_5")

# Not using min size
#------------------------------------------------

#?cor
# cor(datab$length_min,datab$length_max)
# # 0.8661344 (191 species)

# (with 99 species : 0.8714852)

#######################################
# Using 4 qualitative traits
########################################

head(datab)
dataf<-datab[,c(3:13)]
head(dataf)

summary(dataf)
#          feeding         spawning          DVM   ize_class       epi             meso           bathy        
# Active Ambush:28   Broadcaster:52   No     :57   Size_1:86   Mode :logical   Mode :logical   Mode :logical  
# Cruise       :34   Sac-spawner:84   Reverse: 5   Size_2:51   FALSE:16        FALSE:79        FALSE:143      
# Filter       :55   NA's       :55   Strong :13   Size_3:42   TRUE :173       TRUE :110       TRUE :46       
# Mixed        :13                    Weak   :50   Size_4:12   NA's :2         NA's :2         NA's :2        
# NA's         :61                    NA's   :66                                                              
#
# C               O               D               H          
# Mode :logical   Mode :logical   Mode :logical   Mode :logical  
# FALSE:112       FALSE:54        FALSE:146       FALSE:123      
# TRUE :62        TRUE :120       TRUE :28        TRUE :51       
# NA's :17        NA's :17        NA's :17        NA's :17



###############################################################################
# MCA computation
# With function MCA of FactoMineR Keep only 4 qualitative traits 
###############################################################################

# Number of MCA dimensions
#-----------------------------------
cats <- apply(dataf[,-c(3,5:7)], 2, function(x) nlevels(as.factor(x)))
cats
# feeding   spawning size_class          C          O          D          H 
# 4          2          4          2          2          2          2 
nb_axes <- sum(cats)-length(cats)
nb_axes
# 11 axes


# Computing the MCA with FactoMineR With 92 speices as supl.
#------------------------------------------------------------


known <-which(traits$na3==0)
unknown <-which(traits$na3>0)
dataf.MCA <- MCA(dataf[,1:11], ncp=6, na.method="Average",quali.sup=c(3,5:7),ind.sup =unknown)


# Looking at the eigenvalues
#---------------------------
eig <- data.frame(prop=dataf.MCA$eig$"percentage of variance",nb=c(1:nb_axes))
ggplot(eig) + geom_bar(aes(x=nb,y=prop),stat="identity") +
  geom_line(aes(x=nb,y=mean(prop))) +
  xlab("Eigenvalue Number") +  ylab("Value")

#DIM <- 5 # 5 axes -> 80.66 %
DIM <- 4 # 4 axes -> 70.77 %

summary(dataf.MCA)
dimdesc(dataf.MCA)
head(dataf.MCA$eig)
plot(dataf.MCA)
dataf.MCA$eig

# Getting the results
#---------------------
mca1 <-paste0("MCA 1 (",floor(eig$prop[1]*100)/100,"%)")
mca2 <-paste0("MCA 2 (",floor(eig$prop[2]*100)/100,"%)")
mca3 <-paste0("MCA 3 (",floor(eig$prop[3]*100)/100,"%)")
mca4 <-paste0("MCA 4 (",floor(eig$prop[4]*100)/100,"%)")

mca_var <- data.frame(dataf.MCA$var$coord, Variable = rep(names(cats), cats),Class=row.names(dataf.MCA$var$coord))
mca_sp <- data.frame(dataf.MCA$ind$coord)
mca_sp_sup <- data.frame(dataf.MCA$ind.sup$coord)

colnames(mca_var) <- c("MCA1","MCA2","MCA3","MCA4","MCA5","MCA6",#"MCA7","MCA8","MCA9",
                       "Variable","Class")
colnames(mca_sp) <- c("MCA1","MCA2","MCA3","MCA4","MCA5","MCA6")
colnames(mca_sp_sup) <- c("MCA1","MCA2","MCA3","MCA4","MCA5","MCA6")

# Explaining the MCA axes
#----------------------------
 temp1 <- melt(dataf.MCA$var$eta2)
 ggplot(temp1,aes(x=Var2,y=value,color=Var1))+geom_point(size=10)
# temp4 <- melt(data.frame(dataf.MCA$var$contrib, Variable = rep(names(cats), cats),Class=row.names(dataf.MCA$var$coord)),id=c("Class","Variable"))
# ggplot(temp4,aes(x=variable,y=value,color=Variable))+geom_point(size=4)+
#   geom_text(aes(label=Class),hjust=-0.2,vjust=0.2)
temp3 <- melt(mca_var,id=c("Class","Variable"))
ggplot(temp3,aes(x=variable,y=value,color=Variable))+geom_point(size=4)+
  geom_text(aes(label=Class),hjust=-0.2,vjust=0.2)
ggplot(temp3,aes(x=variable,y=value,color=Variable))+
  geom_text(aes(label=Class))

# plot of variable categories
ggplot(data=mca_var, aes(x=MCA1, y=MCA2)) +   theme_bw()  +
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(label=Class,colour=Variable)) +
  geom_point(data=mca_sp,aes(x=MCA1, y=MCA2))+
  geom_point(data=mca_sp_sup,aes(x=MCA1, y=MCA2),color="grey")+
  ggtitle("MCA plot of variables using R package FactoMineR")+xlab(mca1)+ylab(mca2)

ggplot(data=mca_var, aes(x=MCA3, y=MCA4)) +   theme_bw()  +
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(label=Class,colour=Variable)) +
  geom_point(data=mca_sp,aes(x=MCA3, y=MCA4))+
  geom_point(data=mca_sp_sup,aes(x=MCA3, y=MCA4),color="grey")+
  ggtitle("MCA plot of variables using R package FactoMineR")+xlab(mca3)+ylab(mca4)


# cluster sur MCA coordinates on 4 axes and 99 species 
#########################################################
color_cluster <- c("1"="purple", "2"="blue", "3"="green", "4"="grey","5"="orange","6"="red","7"="black","8"="darkred","9"="darkblue","10"="darkgreen","11"="grey80","12"="grey50")

dataf_99 <- dataf[known,]
mca_distance_99 <- dist(mca_sp[,1:DIM], method = "euclidean") # distance matrix
fit_mca_99 <- hclust(mca_distance_99, method= "ward.D") 
plot(fit_mca_99, hang= -1, labels = rownames(mca_sp)) # display dendogram
CLUSTER <- 6
rect.hclust(fit_mca_99, CLUSTER)
mca_groups_99 <- cutree(fit_mca_99, k = CLUSTER) # cut tree into 4 clusters
mca_groups_99 <- revalue(as.factor(mca_groups_99), c("5"="1","2"="2", "1"="3","3"="4", "6"="5", "4"="6"))
mca_groups_99 <- relevel(relevel(relevel(relevel(relevel(relevel(mca_groups_99, "6"), "5"), "4"), "3"), "2"), "1")
dataf_99$cluster_mca <- as.factor(mca_groups_99)


dend_99 <- as.dendrogram(fit_mca_99)
labels_colors(dend_99) <- color_cluster[dataf_99$cluster_mca][order.dendrogram(dend_99)]
plot(dend_99)
table(dataf_99[,1],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,2],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,4],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,8],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,9],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,10],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,11],mca_groups_99,useNA = c("ifany"))


#########################################################
# cluster sur MCA coordinates on 4 axes and 191 species 
#########################################################
color_cluster <- c("1"="purple", "2"="blue", "3"="green", "4"="grey","5"="orange",
                   "6"="red","7"="black","8"="darkred","9"="darkblue","10"="darkgreen","11"="grey80","12"="grey50")

mca_all_temp <- rbind(mca_sp[,1:DIM],mca_sp_sup[,1:DIM])
mca_all_sp <- mca_all_temp[ order(row.names(mca_all_temp)), ]

mca_distance <- dist(mca_all_sp, method = "euclidean")
# mca_distance <- dist(mca_sp[,1:DIM], method = "euclidean") # distance matrix
fit_mca <- hclust(mca_distance, method= "ward.D") 
par(cex=0.6); plot(fit_mca, hang= -1, label=c(1:nn[1])) # display dendogram with numbers
plot(fit_mca, hang= -1, labels = traits$species) # display dendogram with names
CLUSTER <- 6
rect.hclust(fit_mca, CLUSTER)

# Functional groups
mca_groups <- cutree(fit_mca, k = CLUSTER) # cut tree into clusters
mca_groups <- revalue(as.factor(mca_groups), c("2"="1","3"="2", "1"="3","4"="4", "6"="5", "5"="6"))
mca_groups <- relevel(relevel(relevel(relevel(relevel(relevel(mca_groups, "6"), "5"), "4"), "3"), "2"), "1")
dataf$cluster_mca <- as.factor(mca_groups)
traits$cluster_mca <- as.factor(mca_groups)

# write.csv(dataf,"../DATA_csv_xls_tre/Sakina_dataf.csv")
# write.csv(traits,"../DATA_csv_xls_tre/Sakina_traits.csv")

# Plot dendro with color
dend <- as.dendrogram(fit_mca)
labels_colors(dend) <- color_cluster[dataf$cluster_mca][order.dendrogram(dend)]
plot(dend)



# Characterize functional groups
#-------------------------------
table(dataf[,1],mca_groups,useNA = c("ifany"))
table(dataf[,2],mca_groups,useNA = c("ifany"))
table(dataf[,4],mca_groups,useNA = c("ifany"))
table(dataf[,8],mca_groups,useNA = c("ifany"))
table(dataf[,9],mca_groups,useNA = c("ifany"))
table(dataf[,10],mca_groups,useNA = c("ifany"))
table(dataf[,11],mca_groups,useNA = c("ifany"))



# AVERAGE method

?HCPC
# Performs an agglomerative hierarchical clustering on results from a factor analysis.
# It is possible to cut the tree by clicking at the suggested (or an other) level.
# Results include paragons, description of the clusters, graphics.
dataf.MCA4 <- MCA(dataf[,1:11], ncp=4, na.method="Average",quali.sup=c(3,5:7),ind.sup =unknown)
h1 <- HCPC(dataf.MCA4,metric="euclidean",method="ward")
h2 <- HCPC(dataf.MCA4,metric="euclidean",method="average")

dim(h1$data.clust)



fit_mca2 <- hclust(mca_distance, method= "average") 
par(cex=0.6); plot(fit_mca2, hang= -1, label=c(1:nn[1])) # display dendogram with numbers
plot(fit_mca2, hang= -1, labels = traits$species) # display dendogram with names
CLUSTER <- 6
rect.hclust(fit_mca2, CLUSTER)

# Functional groups
mca_groups2 <- cutree(fit_mca2, k = CLUSTER) # cut tree into clusters
mca_groups2 <- revalue(as.factor(mca_groups2), c("2"="1","3"="2", "1"="3","4"="4", "6"="5", "5"="6"))
mca_groups2 <- relevel(relevel(relevel(relevel(relevel(relevel(mca_groups2, "6"), "5"), "4"), "3"), "2"), "1")
dataf$cluster_mca2 <- as.factor(mca_groups2)
traits$cluster_mca2 <- as.factor(mca_groups2)

# write.csv(dataf,"../DATA_csv_xls_tre/Sakina_dataf.csv")
# write.csv(traits,"../DATA_csv_xls_tre/Sakina_traits.csv")

# Plot dendro with color
dend2 <- as.dendrogram(fit_mca2)
labels_colors(dend2) <- color_cluster[dataf$cluster_mca2][order.dendrogram(dend2)]
plot(dend2)



# Characterize functional groups
#-------------------------------
table(dataf[,1],mca_groups,useNA = c("ifany"))
table(dataf[,2],mca_groups,useNA = c("ifany"))
table(dataf[,4],mca_groups,useNA = c("ifany"))
table(dataf[,8],mca_groups,useNA = c("ifany"))
table(dataf[,9],mca_groups,useNA = c("ifany"))
table(dataf[,10],mca_groups,useNA = c("ifany"))
table(dataf[,11],mca_groups,useNA = c("ifany"))



#-------------------------------

color_troph <- c("O"="yellow", "C"="red", "O-H"="green", "O-D"="blue","O-C"="orange","NA"="black")
labels_colors(dend) <- color_troph[dataf$trophism][order.dendrogram(dend)]
plot(dend)
rect.hclust(fit_mca, CLUSTER)

color_feeding <- c("Active Ambush"="red", "Cruise"="green", "Filter"="blue","Mixed"="orange")
labels_colors(dend) <- color_troph[traits$feeding][order.dendrogram(dend)]
plot(dend)

color_size <- c("Size_1"="green", "Size_2"="yellow", "Size_3"="orange", "Size_4"="red","Size_5"="black")
labels_colors(dend) <- color_size[dataf$lmax][order.dendrogram(dend)]
plot(dend)

color_sp <- c("Broadcaster"="green", "spawner"="yellow")
labels_colors(dend) <- color_size[dataf$spawning][order.dendrogram(dend)]
plot(dend)


# Characterize functional groups
#-------------------------------
table(dataf_99[,1],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,2],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,4],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,8],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,9],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,10],mca_groups_99,useNA = c("ifany"))
table(dataf_99[,11],mca_groups_99,useNA = c("ifany"))



plotphy <- as.phylo(fit_mca)
plot(plotphy, show.tip.label=TRUE)
tiplabels(plotphy$tip.label,pch=22, col=(dataf$trophism), bg=(dataf$trophism), cex=0.2) 

plot(plotphy, show.tip.label=TRUE)
tiplabels(plotphy$tip.label,pch=22, col=(dataf$lmax), bg=(dataf$lmax), cex=0.2) 

plot(plotphy, show.tip.label=TRUE)
tiplabels(plotphy$tip.label,pch=22, col=(dataf$spawning), bg=(dataf$spawning), cex=0.2) 

plot(plotphy, show.tip.label=TRUE)
tiplabels(plotphy$tip.label,pch=22, col=(dataf$feeding), bg=(dataf$feeding), cex=0.2) 

plot(plotphy, show.tip.label=TRUE)
tiplabels(plotphy$tip.label,pch=22, col=(dataf$layer), bg=(dataf$layer), cex=0.2) 


par(mfrow=c(2,3))
for (i in 1:nrow(my.sample)) {
  plot(my.phylo, show.tip.label=FALSE, main=paste("Phylogeny - Bioregion",i))
  tiplabels(pch=22, col=(my.sample[i,]>0)+1, bg=(my.sample[i,]>0), cex=0.2) 
}
par(mfrow=c(1,1))



# Caracterizing each functional group
#------------------------------------
summary(filter(dataf,cluster_mca==1)[,c(1,2,4,8:11)])
# feeding          spawning   size_class    C               O               D               H          
# Active Ambush: 0   Broadcaster: 7   Size_1: 0   Mode:logical   Mode :logical   Mode :logical   Mode :logical  
# Cruise       :11   Sac-spawner:16   Size_2: 0   TRUE:31        FALSE:25        FALSE:31        FALSE:31       
# Filter       : 4   NA's       :10   Size_3:30   NA's:2         TRUE :6         NA's :2         NA's :2        
# Mixed        : 0                    Size_4: 3                  NA's :2                                        
# NA's         :18                                                                                              
summary(filter(dataf,cluster_mca==2)[,c(1,2,4,8:11)])
# feeding          spawning   size_class    C               O               D               H          
# Active Ambush:16   Broadcaster: 7   Size_1:12   Mode:logical   Mode :logical   Mode :logical   Mode :logical  
# Cruise       : 5   Sac-spawner:20   Size_2:20   TRUE:30        FALSE:29        FALSE:30        FALSE:30       
# Filter       : 0   NA's       : 6   Size_3: 1   NA's:3         TRUE :1         NA's :3         NA's :3        
# Mixed        : 0                    Size_4: 0                  NA's :3                                        
#  NA's         :12                                                                                              
summary(filter(dataf,cluster_mca==3)[,c(1,2,4,8:11)])
# feeding          spawning   size_class     C              O               D               H          
# Active Ambush: 0   Broadcaster:11   Size_1:6    Mode :logical   Mode:logical   Mode :logical   Mode :logical  
# Cruise       : 0   Sac-spawner: 0   Size_2:6    FALSE:10        TRUE:11        FALSE:11        FALSE:7        
# Filter       : 0   NA's       : 1   Size_3:0    TRUE :1         NA's:1         NA's :1         TRUE :4        
#  Mixed        :11                    Size_4:0    NA's :1                                        NA's :1        
#  NA's         : 1                                                                                              
summary(filter(dataf,cluster_mca==4)[,c(1,2,4,8:11)])
# feeding          spawning   size_class     C              O               D               H          
# Active Ambush: 0   Broadcaster:24   Size_1:16   Mode :logical   Mode:logical   Mode :logical   Mode :logical  
# Cruise       : 3   Sac-spawner: 3   Size_2:17   FALSE:48        TRUE:48        FALSE:48        FALSE:12       
# Filter       :43   NA's       :23   Size_3:10   NA's :2         NA's:2         NA's :2         TRUE :36       
# Mixed        : 0                    Size_4: 7                                                  NA's :2        
#  NA's         : 4                                                                                              
summary(filter(dataf,cluster_mca==5)[,c(1,2,4,8:11)])
# feeding          spawning   size_class     C              O               D               H          
# Active Ambush:11   Broadcaster: 0   Size_1:20   Mode :logical   Mode:logical   Mode :logical   Mode :logical  
# Cruise       : 0   Sac-spawner:17   Size_2: 0   FALSE:11        TRUE:11        FALSE:11        FALSE:11       
# Filter       : 0   NA's       : 3   Size_3: 0   NA's :9         NA's:9         NA's :9         NA's :9        
#  Mixed        : 0                    Size_4: 0                                                                 
#  NA's         : 9                                                                                              
summary(filter(dataf,cluster_mca==6)[,c(1,2,4,8:11)])
# feeding          spawning   size_class     C              O               D               H          
# Active Ambush: 1   Broadcaster: 3   Size_1:32   Mode :logical   Mode:logical   Mode :logical   Mode :logical  
# Cruise       :15   Sac-spawner:28   Size_2: 8   FALSE:43        TRUE:43        FALSE:15        FALSE:32       
# Filter       : 8   NA's       :12   Size_3: 1   NA's :0         NA's:0         TRUE :28        TRUE :11       
#  Mixed        : 2                    Size_4: 2                                  NA's :0         NA's :0        
#  NA's         :17                                                                                              



summary(filter(traits,cluster_mca==1)[,c(7,8,11,12)])
summary(filter(traits,cluster_mca==2)[,c(7,8,11,12)])
summary(filter(traits,cluster_mca==3)[,c(7,8,11,12)])
summary(filter(traits,cluster_mca==4)[,c(7,8,11,12)])
summary(filter(traits,cluster_mca==5)[,c(7,8,11,12)])
summary(filter(traits,cluster_mca==6)[,c(7,8,11,12)])

################################################
# Final Plot by cluster on 4 MCA axes
###################################################



ggplot(mca_sp,aes(x=MCA1,y=MCA2)) +   theme_bw()  +
  geom_point(aes(color=dataf$cluster_mca,shape=traits$order),size=5)+
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(label=c(1:nn[1]),color=dataf$cluster_mca),
            hjust=-0.4, vjust=-0.4,size=5)+
  geom_text(data=mca_var,aes(x=MCA1,y=MCA2,label=Class),size=7,family="Times",fontface="bold") +
  xlab(mca1) + ylab(mca2) + 
  labs(title="MCA on 4 copepod traits",color="Functional group",shape="Order") +
  coord_cartesian()


ggplot(mca_sp,aes(x=MCA3,y=MCA4)) +   theme_bw()  +
  geom_point(aes(color=dataf$cluster_mca,shape=traits$order),size=6)+
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(label=c(1:nn[1]),color=dataf$cluster_mca),
            hjust=-0.4, vjust=-0.4,size=6)+
  geom_text(data=mca_var,aes(x=MCA3,y=MCA4,label=Class),size=9,family="Times",fontface="bold") +
  xlab(mca3) + ylab(mca4) + 
  labs(title="MCA on 4 copepod traits",color="Functional group",shape="Order")       


ggplot(mca_sp,aes(x=MCA1,y=MCA2)) +   theme_bw()  +
  geom_point(aes(color=dataf$cluster_mca,shape=traits$order),size=6)+
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(label=c(1:nn[1]),color=dataf$cluster_mca),
            hjust=-0.4, vjust=-0.4,size=6)+
  geom_text(data=mca_var,aes(x=MCA1,y=MCA2,label=Class),size=9,family="Times",fontface="bold") +
  xlab(mca1) + ylab(mca2) + 
  labs(title="MCA on 4 copepod traits",color="Functional group",shape="Order")       


ggplot(mca_sp,aes(x=MCA3,y=MCA4)) +   theme_bw()  +
  geom_point(aes(color=dataf$cluster_mca,shape=traits$order),size=6)+
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(label=c(1:nn[1]),color=dataf$cluster_mca),
            hjust=-0.4, vjust=-0.4,size=6)+
  geom_text(data=mca_var,aes(x=MCA3,y=MCA4,label=Class),size=9,family="Times",fontface="bold") +
  xlab(mca3) + ylab(mca4) + 
  labs(title="MCA on 4 copepod traits",color="Functional group",shape="Order")       


# List species by functional groups
####################################

traits$nb <- c(1:nn[1])
write.csv(traits,"../DATA_csv_xls_tre/Proteus_species_traits_groups.csv")
          
# MCA with function dudi.acm of ade4
####################################
library(ade4)
mca3 = dudi.acm(dataf[,c(1:4,6)], scannf = FALSE, nf = 5)
# eigenvalues
head(mca3$eig/sum(mca3$eig)*100)


# With function mca of MASS
#############################

dataf.mca <- mca(dataf[,c(1:4,6)], abbrev=TRUE,nf=19)   # To use the abbrevations of the factors
dataf.mca
# Multiple correspondence analysis of 168 cases of 5 factors
# Correlations 0.647 0.574 0.560 0.527  cumulative % explained 16.17 30.53 44.54 57.71 
# Axis1: 16.17
# Axis2: 14.36
# Axis3: 14.01
# Axis4: 13.17

# eigenvalues
dataf.mca$d^2
head(dataf.mca$d^2*100/sum(dataf.mca$d^2))

mca1 <- "MCA 1 (16.17%)"
mca2 <- "MCA 2 (14.36%)"
mca3 <- "MCA 3 (14.01%)"
mca4 <- "MCA 4 (13.17%)"

  
# plot(dataf.mca)
mcf_sp <- data.frame(dataf.mca$rs)
mcf_var <- data.frame(dataf.mca$cs)
colnames(mcf_sp) <- c("MCA1","MCA2","MCA3","MCA4")
colnames(mcf_var) <- c("MCA1","MCA2","MCA3","MCA4")
bli <- data.frame(mcf_var)
rownames(mcf_var)

# cluster sur MCA coordinates on 4 axes: 

mca_distance2 <- dist(mcf_sp, method = "euclidean") # distance matrix
fit_mca <- hclust(mca_distance2, method= "ward.D") 
# clustering on the distance matrix computed from the MCA, aggregation method = Ward's
plot(fit_mca, hang= -1, labels = traits$species) # display dendogram
mca_groups <- cutree(fit_mca, k = 5) # cut tree into 5 clusters
mca_groups[which(mca_groups==1)] <- 5
mca_groups[which(mca_groups==4)] <- 1
mca_groups[which(mca_groups==3)] <- 4
mca_groups[which(mca_groups==2)] <- 3
mca_groups[which(mca_groups==5)] <-2
traits$cluster_mca <- mca_groups
dataf$cluster_mca <- mca_groups

nn <- dim(dataf.mca$rs)
par(cex=0.6); plot(fit_mca, hang= -1, label=c(1:nn[1])) # display dendogram
par(cex=0.4); plot(fit_mca, hang= -1, label=c(1:nn[1])) # display dendogram

rect.hclust(fit_mca, 5)


# by cluster on 4 MCA axes
ggplot(mcf_sp,aes(x=MCA1,y=MCA2)) + 
  geom_point(aes(color=as.factor(dataf$cluster_mca),shape=traits$order),size=6)+
  #xlim(-0.04,0.035)+
  geom_text(aes(label=c(1:nn[1]),color=as.factor(dataf$cluster_mca)),
            hjust=-0.4, vjust=-0.4,size=7)+
  geom_text(data=bli,aes(label=rownames(bli)),size=12,family="Times",fontface="bold") +
  theme_bw()  +
  xlab(mca1) + ylab(mca2) + labs(title="MCA on 4 copepod traits")       

ggplot(mcf_sp,aes(x=MCA3,y=MCA4)) + 
  geom_point(aes(color=as.factor(dataf$cluster_mca),shape=traits$order),size=6)+
  #xlim(-0.04,0.035)+
  geom_text(aes(label=c(1:nn[1]),color=as.factor(dataf$cluster_mca)),
            hjust=-0.4, vjust=-0.4,size=4)+
  geom_text(data=bli,aes(label=rownames(bli)),size=8) + theme_bw()  +
  xlab(mca3) + ylab(mca4) + labs(title="MCA on 4 copepod traits")       


# 2 axes, by superfamily
ggplot(mcf_sp) + geom_point(aes(x=dataf.mca$rs[,1],y=dataf.mca$rs[,2],color=traits$superfamily,shape=traits$order2),size=5)+
  xlim(-0.03,0.04)+
  geom_text(aes(x=dataf.mca$rs[,1],y=dataf.mca$rs[,2],label=rownames(dataf.mca$rs),color=as.factor(traits$superfamily)),hjust=-0.2, vjust=-0.2,size=7)+
  #geom_text(data=bli,aes(x=bli$X1,y=bli$X2,label=rownames(bli)),size=10,family="Times",fontface="bold") + 
  theme_bw()  +
  xlab("MCA1 (21.42%)") +
  ylab("MCA2 (19.43%)") + labs(title="MCA on 4 copepod traits")  



# Caracterizing each functional group :

 
summary(dataf[which(dataf$cluster_mca==1),1:4])
# feeding   migration_type        spawning   lmax   
# C   :46   m_n : 8        Broadcaster: 1   S_1:22  
# O   : 1   m_w :11        Sac-spawner:44   S_2:16  
# FF  : 1   m_s : 0        NA's       : 5   S_3:12  
#  NA's: 2   m_r : 5                         S_4: 0  
# NA's:26                         S_5: 0  
summary(dataf[which(dataf$cluster_mca==2),1:4])
#  feeding migration_type        spawning   lmax   
#  C : 5   m_n :21        Broadcaster:29   S_1:17  
#  O :22   m_w :28        Sac-spawner:12   S_2:21  
#  FF:39   m_s : 0        NA's       :25   S_3:25  
# m_r : 0                         S_4: 0  
# NA's:17                         S_5: 3  
summary(dataf[which(dataf$cluster_mca==3),1:4])
#  feeding   migration_type        spawning   lmax   
#  C   :23   m_n :13        Broadcaster:15   S_1: 0  
#  O   : 9   m_w :11        Sac-spawner: 9   S_2: 0  
#  FF  : 5   m_s : 4        NA's       :14   S_3: 0  
# NA's: 1   m_r : 0                         S_4:38  
#            NA's:10                         S_5: 0  
summary(dataf[which(dataf$cluster_mca==4),1:4])
# feeding migration_type        spawning   lmax   
# C :23   m_n :15        Broadcaster: 7   S_1: 0  
# O : 6   m_w : 0        Sac-spawner:19   S_2: 1  
# FF: 8   m_s : 9        NA's       :11   S_3: 1  
#          m_r : 0                         S_4: 0  
#          NA's:13                         S_5:35  



summary(traits[which(traits$cluster_mca==1),1:5])
summary(traits[which(traits$cluster_mca==2),1:5])
summary(traits[which(traits$cluster_mca==3),1:5])
summary(traits[which(traits$cluster_mca==4),1:5])














summary(dataf[which(dataf$cluster_mca6==1),1:5])
summary(dataf[which(dataf$cluster_mca6==2),1:5])
summary(dataf[which(dataf$cluster_mca6==3),1:5])
summary(dataf[which(dataf$cluster_mca6==4),1:5])
summary(dataf[which(dataf$cluster_mca6==5),1:5])
summary(dataf[which(dataf$cluster_mca6==6),1:5])
summary(dataf[which(dataf$cluster_mca6==7),1:5])
summary(dataf[which(dataf$cluster_mca6==8),1:5])
















ggplot(dataf) + geom_histogram(aes(x=lmin,color=as.factor(feeding),fill=migration_type,alpha=as.factor(spawning)))+ facet_wrap(~cluster_mca)


ggplot(data) + geom_histogram(aes(x=length_min,fill=as.factor(lmin)))
ggplot(data) + geom_histogram(aes(x=length_max,fill=as.factor(lmax)))

# write.csv(traits,file='species_by_cluster_2015-06-03.csv')



# 3axes
ggplot(mcf_sp) + geom_point(aes(x=mcf_sp[,1],y=mcf_sp[,2],colour=as.factor(mcf_sp[,3]>0)))+#,shape=traits$order2),size=5)+
  xlim(-0.04,0.06)+
  geom_text(aes(x=mcf_sp[,1],y=mcf_sp[,2],label=rownames(mcf_sp),colour=as.factor(mcf_sp[,3]>0)),hjust=-0.2, vjust=-0.2,size=3)+
  geom_text(data=bli,aes(x=bli$X1,y=bli$X2,label=rownames(bli),colour=as.factor(bli$X3>0)),size=8)+
  xlab("MCA1 (16.76%)") +
  ylab("MCA2 (15.25%)") +labs(title="MCA on copepod traits (MCA3 : 13.52%)")          



# 2 axes, by order
ggplot(mcf_sp) + geom_point(aes(x=dataf.mca$rs[,1],y=dataf.mca$rs[,2],color=traits$order2,shape=traits$order2),size=5)+
  xlim(-0.04,0.035)+
  geom_text(aes(x=dataf.mca$rs[,1],y=dataf.mca$rs[,2],label=rownames(dataf.mca$rs),color=as.factor(traits$order2)),hjust=-0.2, vjust=-0.2,size=3)+
  geom_text(data=bli,aes(x=bli$X1,y=bli$X2,label=rownames(bli)),size=8)+
  xlab("MCA1 (24.53%)") +
  ylab("MCA2 (22.21%)") +labs(title="MCA on copepod traits")          

# ggplot(mcf_sp) + geom_point(aes(x=dataf.mca$rs[,1],y=dataf.mca$rs[,2],color=as.factor(traits$cluster),shape=traits$order2),size=5)+
#   xlim(-0.04,0.035)+
#   geom_text(aes(x=dataf.mca$rs[,1],y=dataf.mca$rs[,2],label=rownames(dataf.mca$rs),color=as.factor(traits$cluster)),hjust=-0.2, vjust=-0.2,size=3)+
#   geom_text(data=bli,aes(x=bli$X1,y=bli$X2,label=rownames(bli)),size=8)+
#   xlab("MCA1 (24.53%)") +
#   ylab("MCA2 (22.21%)") +labs(title="MCA on copepod traits")          


# MCA with FactoMinR
#----------------------------------
?MCA

datafc <- cbind(dataf,traits$cluster)
dataf.MCA <- MCA(datafc, ncp = 3, quanti.sup = 5, 
    graph = TRUE, axes = c(1,2),  method="Burt",#"Indicator"
    na.method="Average")#NA")




#-------------------------------
# MCA with  3 size classes
#-------------------------------



# Convert size into 3 quantiles
# so same number of classes in each traits

lminq3 <- quantile(datab$length_min,seq(0,1,.333))
lmaxq3 <- quantile(datab$length_max,seq(0,1,.333))

lmin3 <- cut(datab$length_min,unique(lminq3),include.lowest=TRUE)
lmax3 <- cut(datab$length_max,unique(lmaxq3),include.lowest=TRUE)

dataf3<-cbind(datab,lmin3,lmax3)
dataf3<-dataf3[,-1]
dataf3<-dataf3[,-1]
dataf3

levels(dataf3$feeding) <- list("f_C"="C", "f_F"="FF", "f_O"="O")
levels(dataf3$migration_type) <-   list("m_n"="No", "m_r"="Reverse", "m_s"="Strong", "m_w"="Weak")
levels(dataf3$lmin3) <- c("s_1","s_2","s_3")
levels(dataf3$lmax3)<- c("S_1","S_2","S_3")


dataf3.mca <- mca(dataf3, abbrev=TRUE)   # To use the abbrevations of the factors
dataf3.mca
# Multiple correspondence analysis of 187 cases of 4 factors
# Correlations 0.719 0.625  cumulative % explained 23.95 44.79

# un peu moins bien

plot(dataf3.mca)

mcf3_sp <- as.data.frame(dataf3.mca$rs)
mcf3_var <- as.data.frame(dataf3.mca$cs)
colnames(mcf3_var) <- c("X1","X2")

bli3 <- data.frame(mcf3_var)
rownames(mcf3_var)


# cluster sur MCA coordinates: 

mca_distance23 <- dist(mcf3_sp, method = "euclidean") # distance matrix
fit_mca3 <- hclust(mca_distance23, method= "ward.D") 
# clustering on the distance matrix computed from the MCA, aggregation method = Ward's
plot(fit_mca3, hang= -1, labels = traits$species) # display dendogram
mca3_groups <- cutree(fit_mca3, k = 3) # cut tree into 3 clusters
traits$cluster_mca3 <- mca3_groups

ggplot(mcf3_sp) + geom_point(aes(x=dataf3.mca$rs[,1],y=dataf3.mca$rs[,2],color=as.factor(traits$cluster_mca3),shape=traits$order2),size=5)+
  xlim(-0.04,0.035)+
  geom_text(aes(x=dataf3.mca$rs[,1],y=dataf3.mca$rs[,2],label=rownames(dataf3.mca$rs),color=as.factor(traits$cluster_mca3)),hjust=-0.2, vjust=-0.2,size=3)+
  geom_text(data=bli3,aes(x=bli3$X1,y=bli3$X2,label=rownames(bli3)),size=8)+
  xlab("MCA1 (23.95%)") +
  ylab("MCA2 (20.84%)") +labs(title="MCA on copepod traits (3 size classes)")          



# Gower distance and trait dendrogram : ONLY 4 Traits -> length_min      length_max       feeding   migration_type
#-----------------------------------------------------------------------------------------------


data_gower <- gowdis(datab)
attributes(data_gower)

# Daniel Müllner, fastcluster: Fast Hierarchical, Agglomerative Clustering Routines for R and Python,
# Journal of Statistical Software 53 (2013), no. 9, 1–18, URL http://www.jstatsoft.org/v53/i09/
  
clust_gower <- fastcluster::hclust(data_gower,method="average")
clust_gower <- stats::hclust(data_gower,method="average")
summary(clust_gower)
par(cex=0.8);plot(clust_gower)
rect.hclust(clust_gower, 4)



## Colored dendrogram
#----------------------------------

length(unique(traits$superfamily))
length(unique(traits$order))

# vectors of colors 
color_sf <-  c(brewer.pal(12,"Set3"), brewer.pal(9,"Set1") )
labelColors <- c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0")


dend <- as.dendrogram(clust_gower)
# let's add some color and sort them based on their order in dend:
#colors_to_use <- as.numeric(traits$order)

# by MCA cluster:
colors_to_use <- as.numeric(traits$cluster_mca)

colors_to_use <- colors_to_use[order.dendrogram(dend)]
# Now we can use them
labels_colors(dend) <- labelColors[colors_to_use]
plot(dend, main = "A color for each MCA cluster")

colors_to_use_sf <- as.numeric(traits$superfamily)
colors_to_use_sf <- colors_to_use[order.dendrogram(dend)]
# Now we can use them
dendsf <- dend
labels_colors(dendsf) <- color_sf[colors_to_use_sf]
plot(dendsf, main = "A color for each Super Family")





#----------------------------------
# Quick MDS on species traits
#----------------------------------

# Classical MDS

fit <- cmdscale(data_gower,eig=TRUE, k=2) # k is the number of dim
fit # view results

res_mds <- data.frame(fit$points)

ggplot(res_mds) + geom_point(aes(x=X1,y=X2,color=traits$order,shape=traits$order))

ggplot(res_mds) + geom_point(aes(x=X1,y=X2,color=traits$superfamily,shape=traits$order))+
  xlim(-0.65,0.5)+
  geom_text(aes(x = res_mds$X1, y = res_mds$X2,label=rownames(res_mds),color=as.factor(traits$superfamily)),hjust=-0.2, vjust=-0.2,size=4)



ggplot(res_mds) + geom_point(aes(x=X1,y=X2,color=traits$superfamily2,shape=traits$order2,size=10))+
  xlim(-0.65,0.5)+
  geom_text(aes(x = res_mds$X1, y = res_mds$X2,label=rownames(res_mds),color=as.factor(traits$superfamily2)),hjust=-0.2, vjust=-0.2,size=3)


library(data.table)
dt <- data.table(res_mds)
level <- traits$order2
hulls <- dt[, .SD[chull(X1, X2)], by = level]
levelsf <- traits$superfamily2
hullssf <- dt[, .SD[chull(X1, X2)], by = levelsf]

ggplot(dt,aes(x=X1,y=X2,color=level)) +
  geom_point(aes(shape=level)) +
  geom_polygon(data = hulls,aes(fill=level),alpha=0.01)

ggplot(dt,aes(x=X1,y=X2,color=levelsf)) +
  geom_point(aes(shape=traits$order2)) +
  geom_polygon(data = hullssf,aes(x=X1,y=X2,alpha=1),fill=NA)

ggplot(res_mds) + geom_point(aes(x=X1,y=X2,color=traits$superfamily2,shape=traits$order2))+
  xlim(-0.65,0.5)+
  #geom_text(aes(x = res_mds$X1, y = res_mds$X2,label=rownames(res_mds),color=as.factor(traits$superfamily2)),hjust=-0.2, vjust=-0.2,size=3)+
  #geom_polygon(data = hullssf,aes(x=X1,y=X2,colour=levelsf,alpha = 0),fill=NA)+
  geom_polygon(data = hulls,aes(x=X1,y=X2,alpha = 0,colour=level),fill=NA,colour="black")




ggplot(res_mds) + geom_point(aes(x=X1,y=X2,color=as.factor(traits$cluster),shape=traits$order2))+
  xlim(-0.65,0.5)



















