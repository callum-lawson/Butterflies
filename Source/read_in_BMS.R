
### Packages

require(reshape)
require(data.table)

### Working directory and functions

setwd("C:/Users/Callum/Documents/Dropbox/PhD/Analyses/BMS")
source("data_processing_functions_10Mar2013.R")

### Read in data

colform <- c("factor","factor","integer","integer","numeric","factor","integer","integer","integer")
bms <- read.table("Callum_Lawson_weeklycounts.txt",sep=",",header=T,	colClasses=colform)
names(bms) <- c("site","gref","east","north","tlen","species","year","week","count")

sum(is.na(bms$tlen))/length(bms$tlen) # 6.6% of transect lengths missing

### Create data list ordered by species 

bsplist <- levels(bms$species)
bspl <- list()
for(s in 1:length(bsplist)){
	bspd2 <- dataprep(bms,bsplist[[s]]) # bms goes straight in
	bspd3 <- subcounts(bspd2,taucounts=2,tauyears=15) # exclude low counts
		# Remove "extinctions" and "colonisations"
	bspl[[s]] <-  bspd3
	}
names(bspl) <- bsplist

with(curdat,sum(na.omit(count[week<24])))

#############
# PHENOLOGY #
#############

mean.nona <- function(x) mean(na.omit(x))

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

pdf("phenology_byweek.pdf",width=6,height=5)

curdat <- bspl$"Silver-spotted Skipper"
mean.nona <- function(x) mean(na.omit(x))
weekmeans <- with(curdat,tapply(count,list(year,week),mean.nona))
weeknames <- colnames(weekmeans)
phencols <- rainbow(length(weeknames))
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
	