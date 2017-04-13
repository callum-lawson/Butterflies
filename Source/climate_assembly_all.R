setwd("C:/Users/Callum/Documents/Dropbox/PhD/Analyses/BMS/Data")

### Read in sites file

sites.full <- read.csv("ukbms_site.csv",header=T)
site.keepnames <- c("NORTH","EAST","SITENO")
sites <- sites.full[,names(sites.full) %in% site.keepnames]

### Convert UKBMS site 1km coordinates to 5km references (rounding)

for (i in 1:length(sites$EAST)){
	ifelse(sites$EAST[i]%%10000<5000,y<-2500,y<-7500)
	sites$east.5k[i]<-(sites$EAST[i]%/%10000)*10000+y
	}
for (i in 1:length(sites$NORTH)){
	ifelse(sites$NORTH[i]%%10000<5000,y<-2500,y<-7500)
	sites$north.5k[i]<-(sites$NORTH[i]%/%10000)*10000+y
	}	

### Climate extraction function

climextract <- function(weathertype){

	require(data.table)
	require(maptools)

	directory <- paste("Met office climate data\\",weathertype,"\\",sep="")

	file.names <- list.files(directory,pattern=".txt")
	final.table <- NULL
	yeartablist <- list()
	sites <- as.data.table(sites)
	# data.table needs coords as factor:
	sites$east.5k <- as.factor(sites$east.5k)
	sites$north.5k <- as.factor(sites$north.5k)

	for (i in 1:length(file.names)){
	
		f <- file.names[i]
		temp.table <- NULL
		cip.tab <- data.frame(readAsciiGrid(paste(directory,f,sep="")))
		if(weathertype=="MeanTemp") names(cip.tab)[1] <- "mean.temp"
		if(weathertype=="Rainfall") names(cip.tab)[1] <- "rainfall"
		cip.tab$year <- substr(f,10,13)
		cip.tab$month <- substr(f,15,16)

		# Convert to data.table (substantially faster):
		names(cip.tab)[names(cip.tab)=="s1"] <- "east.5k"
		names(cip.tab)[names(cip.tab)=="s2"] <- "north.5k"
		# data.table needs coords as factor:
		cip.tab$east.5k <- as.factor(cip.tab$east.5k)
		cip.tab$north.5k <- as.factor(cip.tab$north.5k)
		cip.tab <- as.data.table(cip.tab)

		yeartablist[[i]] <- merge(sites,cip.tab,by=c("east.5k","north.5k"))
	
		cat(f,sep="\n")
	
		} # close i loop
	
	final.table <- rbindlist(yeartablist)
	return(final.table)

	}

### Create and write files

# Temperature
system.time(
	final.temps <- climextract("MeanTemp")
	)	# 9.5 mins

final.temps2 <- final.temps[,
	names(final.temps) %in% c("SITENO","year","month","mean.temp"),
	with=F]

# Rain
system.time(
	final.rain <- climextract("Rainfall")
	)	# 9.3 mins

final.rain2 <- final.rain[,
	names(final.rain) %in% c("SITENO","year","month","rainfall"),
	with=F]

# Combined Temperature and Rain

final.all <- merge(final.temps2,final.rain2,by=c("SITENO","year","month"))
setnames(final.all,"SITENO","site")

write.table(final.all,"all.sites.CIP.mean.tempANDrainfall.1971_2012.txt",
	sep="\t",row.names=F)

remove(final.temps)
remove(final.temps2)
remove(final.rain)
remove(final.rain2)
remove(final.all)
gc()

### Read in ready-made files
#final.temps <- as.data.table(read.table("all.sites.CIP.mean.temp.1971_2012.txt",header=T)) 
#final.rain <- as.data.table(read.table("all.sites.CIP.rain.all.1971_2012.txt",header=T)) 

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

cip.tab <- readAsciiGrid("Met office climate data\\Meantemp_2007_2013\\MeanTemp_2012-08_ACTUAL.txt")

tranblack <- rgb(0,0,0,alpha=0.5)
require("fields")

pdf("bms_site_map_aug2012.pdf",width=12,height=12)
image.plot(cip.tab)
with(sites, points(EAST,NORTH,pch=16,col=tranblack))
dev.off()

