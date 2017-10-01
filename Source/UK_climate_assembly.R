# Pair UKCP09 data with site indices

site <- read.csv("Data/ukbms_site.csv",header=T)
site <- subset(site, !is.na(LATITUDE) & !is.na(LATITUDE))

### Convert UKBMS lat / long to UTM

require(rgdal)
mrc = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'
coordinates(site) <- c("LONGITUDE","LATITUDE")
wgs84 = '+proj=longlat +datum=WGS84'
bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
mrc = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'
site@proj4string   # slot will be empty
site@proj4string = CRS(wgs84)

site_merc = spTransform(site, CRS(mrc)) # reproject data to Mercator coordinate system
bbox = site@bbox # download Mercator tiles covering the data area
site_bng = spTransform(site, CRS(bng)) # conversion to British National Grid (OSGB36)

site$EAST <- site_bng@coords[,1]
site$NORTH <- site_bng@coords[,2]

# supposed to be 543000	110700

# http://blogs.casa.ucl.ac.uk/2013/12/05/british-national-grid-transformation-and-reprojection-in-r/

### Convert UKBMS site 1km coordinates to 5km references (rounding)

for(i in 1:length(site$EAST)){
	ifelse(site$EAST[i]%%10000<5000,y<-2500,y<-7500)
	site$E5[i]<-(site$EAST[i]%/%10000)*10000+y
	}
for(i in 1:length(site$NORTH)){
	ifelse(site$NORTH[i]%%10000<5000,y<-2500,y<-7500)
	site$N5[i]<-(site$NORTH[i]%/%10000)*10000+y
	}	

### Climate extraction function

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

### Create and write files

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

write.table(final.all,"all.site.CIP.mean.tempANDrainfall.1971_2012.txt",
	sep="\t",row.names=F)

remove(final.temps)
remove(final.temps2)
remove(final.rain)
remove(final.rain2)
remove(final.all)
gc()

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

