### Pair UKCP09 data with site indices ###

site <- read.csv("Data/ukbms_site.csv",header=T)
site <- subset(site, !is.na(LATITUDE) & !is.na(LONGITUDE) & !duplicated(site))
  # one site was entered twice

# Convert lat-long to BNG -------------------------------------------------

require(rgdal)
site_rg <- site
coordinates(site_rg) <- c("LONGITUDE","LATITUDE")
wgs84 <- '+proj=longlat +datum=WGS84'
bng <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
site_rg@proj4string = CRS(wgs84)
site_bng = spTransform(site_rg, CRS(bng)) 
  # conversion to British National Grid (OSGB36)

site$EAST <- site_bng@coords[,1]
site$NORTH <- site_bng@coords[,2]

# http://blogs.casa.ucl.ac.uk/2013/12/05/british-national-grid-transformation-and-reprojection-in-r/

# Convert UKBMS site 1km coordinates to 5km references --------------------

for(i in 1:length(site$EAST)){
	ifelse(site$EAST[i]%%10000<5000,y<-2500,y<-7500)
	site$E5[i]<-(site$EAST[i]%/%10000)*10000+y
}
for(i in 1:length(site$NORTH)){
	ifelse(site$NORTH[i]%%10000<5000,y<-2500,y<-7500)
	site$N5[i]<-(site$NORTH[i]%/%10000)*10000+y
}	
  # coordinates at centre points of grid

# Climate extraction function ---------------------------------------------

climextract <- function(climname,ymin=1970){

	require(data.table)
	require(maptools)

	dir <- paste0("Data\\UKCP09\\",climname,"\\")
	
	fnf <- list.files(dir,pattern=".txt")

	fl <- nchar(fnf[1])
	yearfull <- as.numeric(substr(fnf,fl-9,fl-6))
	monthfull <- as.numeric(substr(fnf,fl-5,fl-4))

	year <- yearfull[yearfull>=ymin]
	month <- monthfull[yearfull>=ymin]
	fn <- fnf[yearfull>=ymin]
	
	nf <- length(fn) # number of files
	ytl <- as.list(rep(NA,nf))

	site <- as.data.table(site)
	site$E5 <- as.factor(site$E5)
	site$N5 <- as.factor(site$N5)

	nf <- length(ytl)
	
	for(i in 1:nf){
	    
	  cd <- data.frame(readAsciiGrid(paste0(dir,fn[i])))
	  names(cd) <- c(climname,"E5","N5")
	  cd$year <- year[i]
	  cd$month <- month[i]
	    
	  # data.table needs coords as factor:
	  cd$E5 <- as.factor(cd$E5)
	  cd$N5 <- as.factor(cd$N5)
	  cd <- as.data.table(cd)
	    
	  ytl[[i]] <- merge(site,cd,by=c("E5","N5"))
	    
	  cat(paste(year[i],month[i],collapse="-"),sep="\n")
	  
	} # close i loop
	
	fd <- rbindlist(ytl)
	return(fd)

	}

# Extract climate data ----------------------------------------------------

# Temperature
system.time(
	temp <- climextract("mean-temperature")
	)

# Rain
system.time(
  rain <- climextract("rainfall")
	)

# Combined Temperature and Rain
all <- merge(temp,rain)
names(all) <- tolower(names(all))
all <- subset(all, select=!names(all) %in% c("e5","n5"))
library(plyr)
all <- rename(all, c(
  "latitude"="lat",
  "longitude"="long",
  "mean-temperature"="temp",
  "rainfall"="rain"
  ))

write.table(all,
            file=paste0("Output/UKCP09_",format(Sys.Date(),"%d%b%Y"),".txt"),
            row.names=F
            )

# Checks for consistency with 2013 data -----------------------------------

site0 <- read.csv("Data/2013/ukbms_site.csv",header=T)
names(site0)[names(site0)=="SITENO"] <- "SITE"
sitem <- merge(site,site0,by="SITE",all.x=F,all.y=F)
summary(sitem$NORTH.y-sitem$NORTH.x)
summary(sitem$EAST.y-sitem$EAST.x)
plot(sitem$EAST.y,sitem$NORTH.y,col="red",pch="+")
points(sitem$EAST.x,sitem$NORTH.x,pch="+")
  # lat and long are equivalent, so must be difference in transformation

all0 <- read.table("Data/2013/all.sites.CIP.mean.tempANDrainfall.1971_2012.txt",header=T)

allm <- merge(all,all0,by=c("site","year","month"))
summary(allm$temp-allm$mean.temp)
summary(allm$rain-allm$rainfall)

### Read in ready-made files
#final.temps <- as.data.table(read.table("all.site.CIP.mean.temp.1971_2012.txt",header=T)) 
#final.rain <- as.data.table(read.table("all.site.CIP.rain.all.1971_2012.txt",header=T)) 

###################
# TEMPERATURE MAP #
###################

#  AVERAGING - TO BE COMPLETED
#file.names <- list.files("Met office climate data\\Meantemp_2007_2013",pattern=".txt")
#library(maptools)
#f <- file.names[1]
#cip.array <- array(dim=c(52200,1,length(file.names)))
#for (f in file.names[]){
#	cip.array[,,f] <- readAsciiGrid(
#		paste("Met office climate data\\Meantemp_2007_2013\\",f,sep="")
#		)
#	}

cip.tab <- readAsciiGrid("Data\\UKCP09\\mean-temperature\\ukcp09_gridded-land-obs-monthly_5km_mean-temperature_201612.txt")

tranblack <- rgb(0,0,0,alpha=0.5)
require("fields")

pdf("bms_site_map_aug2012.pdf",width=12,height=12)
image.plot(cip.tab)
with(site, points(EAST,NORTH,pch=16,col=tranblack))
dev.off()

