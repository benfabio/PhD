
##### 03/11/2015 - ETHZ - Fabio Benedetti & Pieter Provoost
##### Script for : 
#					- Retrieve species name list for the 'Neocopepoda' group
#					- Restict the list to planktonic copepods (manual screening with Razouls et al., 2005-2015)
#					- With a mclapply, download all species occurrences from OBIS, and save the datasets, as .Rdata, of those who have more than 
#					  100 occurrences between 0 and 300m depth


### Last update : 03/11/2015


# ---------------------------------------------------------------------------------------------------------------------------------

library("rgdal")
library("raster")
library("Hmisc")
library("ggplot2")
library("stringr")
library("parallel")
library("RPostgreSQL")
library("devtools")

WD <- getwd()

# ---------------------------------------------------------------------------------------------------------------------------------

### 1) Connect to server with the ID given by Pieter
require("RPostgreSQL")

host <- "obisdb-stage.vliz.be"
db <- "obis"
user <- "obisreader"
password <- "0815r3@d3r"
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname=db, host=host, user=user, password=password)

# Get all accepted neocopepod species
all_spp <- dbGetQuery(con, "with neo as (select valid_id, storedpath || valid_id || 'x%' as path from obis.tnames where tname = 'Scyphozoa') select tnames.* from obis.tnames inner join neo on tnames.storedpath like neo.path where tnames.rank_id = 220 and tnames.id = tnames.valid_id")
colnames(all_spp)

head(all_spp)
unique(all_spp$tname)

# Load list of id and sp.names of the chosen spp.: 
list <- read.csv("species_list_neocop_OBIS.csv", h=T, sep=";")
dim(list)
head(list)
summary(list)
list[which(list$flag == "yes"),]

# Restrict all_spp object to the species which have a 'yes' flag in 'spp'
# Identify the IDs that match
i <- match( list[which(list$flag == "yes"), c("id")], all_spp[,c("id")] )
ddf <- all_spp[i,]
dim(ddf)
head(ddf)

### OK ! 


### 2) Define a function that will, for each species id, download the occurrences from OBIS. Print it if it has more than 75 occurrences between 0 and 300m depth. 
###    More selection to come (time, depth...)

# Test: 
#i <- ddf[which(ddf$tname == "Neocalanus gracilis"),"id"]

CleverDownload <- function(i = ids) {

					# Useless message:
					message(paste("Downloading ", ddf[which(ddf$id == i),"tname"], sep=""))
					require(RPostgreSQL)
					host <- "obisdb-stage.vliz.be"
					db <- "obis"
					user <- "obisreader"
					password <- "0815r3@d3r"
					drv <- dbDriver("PostgreSQL")
					con <- dbConnect(drv, dbname=db, host=host, user=user, password=password)
					
					# Choose species according to id: 
					ID <- i
					
					# Get path to file: 
					path <- paste0("'", ddf$storedpath[i], ddf$id[i], "x%'")
					
					# Download data:
					data <- dbGetQuery(con, paste0("select p.* from portal.points_ex p left join obis.tnames t on p.valid_id = t.id where t.id = ", ID, " or t.storedpath like ", path))
					#dim(data)
					
					# Get number of occurrences between surface and 300m depth
					n <- nrow(data[which(data$maximumdepth <= 300),])
					
					# If n >= 75, go to dir and save data as .Rdata: 
					if(n >= 75) { 
						
						setwd(paste(WD,"/","OBIS_files","/", sep=""))
						# Add an underscore to the sp.name
						sp.name <- stringr::str_replace_all(ddf[which(ddf$id == i),"tname"]," ","_")
						# Save as .Rdata
						save(data, file = paste(sp.name,".Rdata", sep=""))	
						
					} else {
						rm(n, path, data)
						gc()
						setwd(WD)
					} # eo if else loop

} #eo FUN 

ids <- unique(ddf$id)
mclapply(X= ids, FUN = CleverDownload, mc.cores = 12)


### 3) 346 copepod species have been retrieved this way. From your computer, load each Rdata within a lapply, rbind and map occurrences. 

# Load global coastline
setwd("~/Desktop/PhD/WOA_2009/WD/WOA_2009_NODC")
cl <- read.csv("world_coast.csv", h=T)

coast <- list(
  # the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  geom_polygon(aes(x=lon, y=lat), data= cl, fill="grey40"),
  geom_path(aes(x=lon, y=lat), data= cl, colour = "black", linetype = 1),
  # appropriate projection
  coord_quickmap(),
  # remove extra space around the coast
  scale_x_continuous(name = "Longitude", 
                     breaks = c(-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180), 
                     labels= c("-180°W","-150°W","-120°W","-90°W","-60°W","-30°W","0°","30°E","60°E","90°E","120°E","150°E","180°E"), 
                     expand=c(0,0)), 
  
  scale_y_continuous(name = "Latitude", 
                     breaks = c(-90,-60,-30,0,30,60,90), 
                     labels= c("-90°S","-60°S","-30°S","0°","30°N","60°N","90°N"),
                     expand=c(0,0)),
  # dark gray background for the panel and legend
  theme(
    panel.background=element_rect(fill="white"),  # background
    legend.key=element_rect(fill="grey50"),
    panel.grid.major=element_line(colour="white")
  )
  
) # eo coast


setwd(WD)


# Load all .Rdata, rbind: 
files <- dir()

res <- lapply(files, function(f) {
				# Load/get
				ddf <- load(f)
				ddf <- get(ddf)
				# Return spp ddf
				return(ddf)
		} # eo fun
	
) # eo lapply

table <- do.call(rbind, res) ; rm(res)

dim(table)
str(table)
colnames(table)
head(table)

# Restrict according to some criteria: 
t <- table[which(table$maximumdepth <= 300 & table$yearcollected >= 1955),]

ggplot() + geom_point(aes(x=longitude, y=latitude), colour= "red", alpha = 0.65, data = t) + coast + coord_quickmap() + theme_linedraw()

### Check the samples' characteristics
summary(table$maximumdepth)
summary(as.numeric(table$yearcollected))
length(unique(table$tname))


res2 <- lapply(files, function(f) {

				# Load/get
				ddf <- load(f)
				ddf <- get(ddf)
				
				ndepth_300 <- nrow(data[which(data$maximumdepth <= 300),])
				
				nmonth <- nrow( data[!is.na(data$monthcollected),] )
				
				nyear <- nrow( data[!is.na(data$yearcollected),] )
				
				meanmaxdepth <- mean(data$maximumdepth, na.rm = TRUE)
				
				meanmonth <- mean(as.numeric(data$monthcollected), na.rm = TRUE)
				
				table <- data.frame(sp.name= f, ndepth_300= ndepth_300, nmonth= nmonth, nyear= nyear, meanmaxdepth= meanmaxdepth, meanmonth= meanmonth)
				
				# Return spp ddf
				return(table)
				
		} # eo fun
	
) # eo lapply

table2 <- do.call(rbind, res2) ; rm(res2)

dim(table2)
head(table2)
summary(table2)

nrow( table2[which(table2$nyear <= 100),] )


### 04/11/2015 : Plot occurrences for each 346 spp.
### Load all .Rdata, rbind: 
files <- dir()[grep(".Rdata",dir())] 
length(files)

for(f in files) {
	
		# Useless message
		message(paste("Doin' ",f, sep=""))

		# Load/get
		ddf <- load(f)
		ddf <- get(ddf)
		
		t <- ddf[which(ddf$maximumdepth <= 300 & ddf$monthcollected %in% c(1:12)),]
		
		plot <- ggplot() + geom_point(aes(x=longitude, y=latitude), colour= "red", alpha = 0.65, data = t) + coast + coord_quickmap() + theme_linedraw()
		
		setwd( paste(WD,"/","_maps_occ","/", sep="") )
		
		sp.name <- stringr::str_replace_all(f,".Rdata","")
		
		ggsave(plot, file= paste(sp.name,".pdf",sep=""), dpi = 300)
		
		setwd(WD)

}


ddf <- read.csv("salpida_obis.csv", h = T, dec = ".")
dim(ddf)
colnames(ddf)
head(ddf)
