### Working directory and functions

#setwd("C:/Users/Callum/Documents/Dropbox/PhD/Analyses/BMS")
setwd("C:/Users/Callum/Documents/My Dropbox/PhD/Analyses/BMS")
source("Code/data_processing_functions_28Apr2013.R")

### Packages

require(reshape2)
require(data.table)

### Read in data

colform <- c("factor","factor","integer","integer","numeric","factor","integer","integer","integer")
bms <- read.table("Data/Callum_Lawson_weeklycounts.txt",sep=",",header=T,	colClasses=colform)
names(bms) <- c("site","gref","east","north","tlen","species","year","week","count")

sum(is.na(bms$tlen))/length(bms$tlen) # 6.6% of transect lengths missing

### Prepare summed count dataframe

SDT2 <- sumframe(bms)

### Create data list ordered by species 

bsplist <- levels(bms$species)
bspl <- list()
for(s in 1:length(bsplist)){
	bspd2 <- dataprep(bms,bsplist[[s]],SDT2) # bms goes straight in
	# bspd3 <- subcounts(bspd2,taucounts=2,tauyears=15) # exclude low counts
		# Remove "extinctions" and "colonisations"
	bspl[[s]] <-  bspd2
	}
names(bspl) <- bsplist

### Merge with site data

# WARNING: 
# MAY NEED TO SORT OUT X AND Y COORDS TO MATCH CLIMATE DATA
#

site <- read.csv("ukbms_site.csv",header=T)
site2 <- site[,names(site) %in% c("SITENO","EAST","NORTH","LENGTH")]
names(site2) <- c("site","east","north","tlen")

site2$site <- as.factor(site2$site)
site3 <- as.data.table(site2)

fspl <- list()
for(i in 1:length(bsplist)){
	bspdat <- as.data.table(bspl[[i]])
	fspl[[i]] <- as.data.frame(merge(bspdat,site2,by="site"))
	}
names(fspl) <- bsplist

#####################
# DISPLAY LOCATIONS #
#####################

gridplot <- function(curdat,...){
	with(curdat,plot(east,north,type="n",...))
	gridseq <- (0:1000)*5000
	for(i in 1:length(gridseq)){
		abline(h=gridseq[i],v=gridseq[i],col="gray")
		}
	with(curdat,points(east,north,pch="+"))
	}

pdf("gridplots.pdf",width=10,height=10)
for(i in 1:length(fspl)){
	curdat <- fspl[[i]]
	onepersite <- curdat[match(levels(curdat$site),curdat$site),]
	gridplot(onepersite,main=bsplist[[i]])
	}
dev.off()

#############
# OCCUPANCY #
#############

curdat <- bspl$"Silver-spotted Skipper"

occ.calc <- function(curdat){
	mbyw <- melt(curdat,id=c("site","year","week"))
	is.occ <- function(x){ 
		if(!(F %in% is.na(x))){
			NA
			}
		else{
			sum(x,na.rm=T) > 0 
			}
		}
	occ <- dcast(mbyw, site + year ~ variable, is.occ)
	return(occ)
	}

ospl <- lapply(bspl,occ.calc)

yearsums <- lapply(ospl, function(occ){
	with(occ, tapply(count,year, function(x) mean(x,na.rm=T)))
	})

pdf("occplot.pdf",width=7*3.5,height=7*3.5)

par(mfrow=c(8,8))
for(i in 1:length(yearsums)){
	yeardata <- yearsums[[i]]
	plot(yeardata~names(yeardata),
		type="l",
		las=1,
		xlab="year",
		ylab="Proportion of sites occupied",
		main=bsplist[[i]])
	}

dev.off()

pdf("occimage.pdf",width=10,height=10)
for(i in 1:length(ospl)){
	occdata <- ospl[[i]]
	occtab <- with(occdata, 
		tapply(count,list(site,year),function(x) mean(x,na.rm=T))
		)
	image(x=as.numeric(colnames(occtab)),
		y=1:nrow(occtab),
		z=t(occtab),col=c("red","blue"),
		xlab="year",
		ylab="site",
		main=bsplist[[i]],
		las=1)
	legend("bottomleft",legend=c(0,1),fill=c("red","blue"),bty="n")
	}
dev.off()

#############
# PHENOLOGY #
#############

### All species - means by week

mean.nona <- function(x) mean(x,na.rm=T)

pdf("mean_phenology.pdf",width=6,height=5)
for(i in 1:length(bspl)){
	curdat <- bspl[[i]]
	if(nrow(curdat)<1) next
	weekmeans <- with(curdat,tapply(count,week,mean.nona))
	plot(weekmeans ~ names(weekmeans),
		type="l",
		xlim=c(1,26),
		main=paste(bsplist[[i]]),
		las=1,
		ylab="Mean count",
		xlab="Week number"
		)
	}
dev.off()

### SSSk by year

curdat <- bspl$"Silver-spotted Skipper"
weekmeans <- with(curdat,tapply(count,list(year,week),mean.nona))
weeknames <- colnames(weekmeans)
phencols <- rainbow(length(weeknames))

pdf("phenology_byweek.pdf",width=6,height=5)
plot(1,1,
	xlim=c(1,26),
	ylim=c(0,max(weekmeans[!is.nan(weekmeans)])),
	type="l",
	ylab="Mean count",
	xlab="Week number",)
for(i in 1:nrow(weekmeans)){
	lines(weekmeans[i,]~weeknames,
		lty=1,col=phencols[i])
	} # close i loop
dev.off()

### Calculating flight periods

# Function to calculate month based on BMS week:
monthcalc <- function(flyweek){
	basetime <- as.POSIXlt("2012 14 01", format = "%Y %W %w")
	baseweek <- format(basetime,"%W")
	sweek <- as.numeric(flyweek) + as.numeric(baseweek)
	stime <-  as.POSIXlt(paste("2012",sweek,"01"), format = "%Y %W %w")
	smonth <- as.numeric(format(stime,"%m"))
	return(smonth)
	}

monthnums <- function(curdat){

	weekmeans <- with(curdat,tapply(count,week,mean.nona))
	weekdens <- weekmeans/sum(weekmeans,na.rm=T)
	weeknames <- names(weekmeans)
	flyweeks <- names(which(weekdens > 0.01)) 
		# less than 1% of total densities

	firstmonth <- flyweeks[1]
	lastmonth <- flyweeks[length(flyweeks)]
		# only works for single flight period at present

	firstnum <- monthcalc(firstmonth)
	lastnum <- monthcalc(lastmonth)

	return(c(first=firstnum,last=lastnum))

	}

monlist <- t(sapply(bspl,monthnums))

write.csv("monlist")

############
# OLD CODE #
############

### separate out sites

bms2 <- bms[,names(bms) %in% c("site","year","week","count")]
sites <- bms[,names(bms) %in% c("site","year","week","count")]
# sort out data so that have first week, last week, duration, etc.

mbms <- melt(bms, measure=c("count"))
cbms <- cast(mbms, site + gref + east + north + tlen + year + week  ~ variable , sum)
sitedat <- cast(mbms, site + gref + east + north  ~ variable, length)
names(sitedat)[names(sitedat)=="count"] <- "nrecords"
# NB: transects always the same length

### Data add function

bms2 <- subset(bms,site %in% c("1018","1019") &
	species %in% c("Silver-spotted Skipper","Adonis Blue")
	)
	
spname <- "Silver-spotted Skipper"

fulldat <- function(spname){
	spdat <- subset(bms,species==spname)
	spdat <- spdat[,names(spdat) %in% c("site","year","week","count")]
	DT <- as.data.table(spdat)
	setkey(DT,site,year,week)
	DT2 <- DT[CJ( 
		unique(site),
		seq(min(year),max(year)),
		seq(min(week),max(week))
		)]
	return(DT2)
	}

sss <- fulldat("Silver-spotted Skipper")

	# merge(DT2,) # work out how to merge info?
	# add in zeroes for surveyed sites (year,week)

### Adding zeroes for surveyed sites

based <- subset(sss,year=="2007" & site=="991")
based <- sss

sites <- levels(as.factor(as.character(based$site)))
years <- unique(based$year)
weeks <- unique(based$week)

nsites <- length(sites)
nyears <- length(years)
nweeks <- length(weeks)

for(i in 1:nsites){
	for(j in 1:nyears){
		curseq <- which(based$site==sites[i] & based$year==years[j])
		basedweek <- based$week[curseq]
		cbmsweek <- cbms$week[cbms$site==sites[i] & cbms$year==years[j]]
		based$count[curseq][basedweek %in% cbmsweek & is.na(based$count[curseq])] <- 0
		}
	}
	