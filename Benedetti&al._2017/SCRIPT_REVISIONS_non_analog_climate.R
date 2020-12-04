
# ---------------------------------------------------------------------------------------------------------------------

library("raster")
library("rgdal")
library("rgeos")
library("ggplot2")
library("stringr")
library("foreign")
library("stringr")
library("plyr")
library("reshape2")
library("parallel")

WD <- getwd()

### Mediterranean coastline :
cl <- read.csv("gshhg_medit.csv")
names(cl) <- c("lon", "lat")
coast <- list(
 	 	# the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  	  	geom_polygon(aes(x=lon, y=lat), data= cl, fill = "grey45"),
  	  	geom_path(aes(x=lon, y=lat), data= cl, colour = "black", linetype = 1),
  	  	# appropriate projection
  	  	coord_map(),
  	  	# remove extra space around the coast
  	  	scale_x_continuous(name = "Longitude", breaks = c(0,10,20,30), labels= c("0°","10°E","20°E","30°E"), 
                    	 expand=c(0,0)), 
  		scale_y_continuous(name = "Latitude", breaks = c(30,35,40,45), labels= c("30°N","35°N","40°N","45°N")
    					, expand=c(0,0)),
  		# dark gray background for the panel and legend
 	   	theme(
	   	 	panel.background=element_rect(fill="white"),  # background
	    		legend.key=element_rect(fill="black"),
	    		panel.grid.major=element_line(colour="white")
  		)
)

# ----------------------------------------------------------------------------------------------------------------------

# List of metadata: scenarios and prevalence levels:
prevalences <- c("n1","n1.5","n2","n10","n50")
scenarios <- c("A2","A2F","A2RF","A2ARF","A1BARF","B1ARF")

### Get global env predictors
varT0_glob <- brick("clims_woa13_baseline_04.grd")
varT0_glob <- dropLayer(varT0_glob, c(1,2,6))
varT0_glob <- as.data.frame(varT0_glob, xy = TRUE)
varT0_glob$id <- paste(varT0_glob$x, varT0_glob$y, sep= "_")

### Species names list:
setwd("./niche_modelling")
second.wd <- getwd()
setwd(paste(second.wd,"/","pabs_files","/", sep=""))
sp.names <- dir() ; length(sp.names)
setwd(second.wd)

# Fix period and pseudo abs run
p <- "2050"
r <- "1"

library("dismo")
#?mess

### Ok, let's test this for 1 species and 1 projection from NEMOMED8...
#sp <- "Pleuromamma_gracilis"
#preval <- "n1"
#s <- "A2"


### Define a FUn for mclapply:
non_analog <- function(sp = sp.names) {

		for(s in scenarios) {	
			
			# Load proper climatology
			setwd(WD)	
			futenv <- read.table(paste("clims",s,p,"final_Levin.txt", sep= "_"), h=T, sep="\t", dec=".")
			futenv <- futenv[,c(1:5)]
			futenv$id <- paste(futenv$x, futenv$y, sep= "_")	
			
			for(preval in prevalences) {
					
					message(paste("Doing MESS for ",sp," ",s," ",preval, sep = ""))
					# Got to proper dir
					setwd( paste(second.wd,"/","pabs_files","/",sp,"/", sep="") )
					# Get presence data
					ddf <- read.csv( paste("pabs_",sp,"_",preval,"_run",r,".txt", sep=""), header = T, sep="\t", dec=".")
					ddf$id <- paste(ddf$x, ddf$y, sep= "_")
					# Ok, now, for each presence point, retrieve corresponding values of average SST, average SSS and ∆SST
					matching <- varT0_glob[which(varT0_glob$id %in% ddf$id),]
					ddf <- cbind(ddf, matching[,c(3:6)]) ; rm(matching)
					# dim(ddf); head(ddf)
					ddf$sp.name <- sp
					
					### Prepare data for MESS
					library("raster")
					# set up an 'empty' raster, here via an extent object derived from your data
					e <- extent(futenv[,1:2])
					#r <- raster(e, ncol=10, nrow=2)
					R <- raster(xmn= -5.91045, xmx= 36.1829, ymn= 30.29393, ymx= 45.71447)
					# you need to provide a function 'fun' for when there are multiple points per cell
					ras <- rasterize(futenv[,1:2], R, futenv[,3:5], fun= mean)
					pr <- projectRaster(ras, R, method='bilinear')
					# pr
					
					### Perform MESS
					mess <- dismo::mess(x = pr, v = ddf[,c(5,6,7)], full = T)
					m <- as.data.frame(mess$rmess, xy = T)
					# dim(m) ; head(m)
					m$id <- paste(m$x, m$y, sep = "_")
					# length(unique(m$id))
					### 64800 rows because of land points, I guess
					# dim(m[which(m$rmess != 'Inf'),])
					M <- m[which(m$rmess != 'Inf'),]
					# head(M) ; summary(M)

					map2 <- ggplot() + geom_point(aes(x = x, y = y), data = M[which(M$rmess < 0),], colour = "red", shape = "+") + coast + 
							scale_colour_gradient2("Multivariate dissimilarity", guide = "colorbar") + 
							coord_quickmap() + theme_classic()

					setwd( paste(WD,"/","maps_non_analogs","/","MESS_T1","/", sep ="") )
					
					if( preval == "n1" && s == "A2" ) {
						dir.create(sp)
						setwd( paste(WD,"/","maps_non_analogs","/","MESS_T1","/",sp,"/", sep ="") )
					} else {
						setwd( paste(WD,"/","maps_non_analogs","/","MESS_T1","/",sp,"/", sep ="") )
					} # eo if else
					
					ggsave( paste("map","_",sp,"_",preval,"_",s,".pdf", sep = ""), plot = map2, width = 8, height = 6)					
					save(M , file= paste(sp,"_",s,"_",preval,".Rdata", sep = "") )
					
					setwd(WD)
					rm(M, map2, m, mess, e, ras, pr, ddf, R)
					gc()
	
			} # eo for preval in prevalences
			
			
		} # eo for s in scenarios


} # eo FUN


library("parallel")
mclapply(X = sp.names, FUN = non_analog, mc.cores = 25)


