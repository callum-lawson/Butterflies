#################################################################
# BMS MONTLY ANALYSIS
# This code analyses population dynamics at a site level, using
# monthly weather variables 
#
#################################################################

##### EXCLUDE SITES ABOVE 10000 ###########
##### BEFORE 1976 - PILOT YEARS ###########

#############################
# READ IN DATA AND PACKAGES #
#############################

### Packages
require(reshape2)

### Site-level data
source("Source/data_processing_functions.R") # Change to most recent version
colform <- c("NULL","factor","factor","factor","integer","integer",rep("NULL",8))
fullind <- read.csv("Data/ukbms_index.csv",header=T,colClasses=colform)
names(fullind) <- c("species","brood","site","year","count")

ind <- subset(fullind,count>=0) # gets rid of other count codes

########################
# ADD TRANSECT LENGTHS #
########################

site <- read.csv("Data/ukbms_site.csv",header=T)
site2 <- site[,names(site) %in% c("SITENO","EAST","NORTH","LENGTH")]
names(site2) <- c("site","east","north","tlen")
for (i in 1:length(site$EAST)){
	ifelse(site$EAST[i]%%10000<5000,y<-2500,y<-7500)
	site2$east[i]<-(site$EAST[i]%/%10000)*10000+y
	}
for (i in 1:length(site$NORTH)){
	ifelse(site$NORTH[i]%%10000<5000,y<-2500,y<-7500)
	site2$north[i]<-(site$NORTH[i]%/%10000)*10000+y
	}
site2$site <- as.factor(site2$site)
site2$gridno <- gridcat(site2,5000)
ind2 <- merge(ind,site2,by="site")

############################
# ADD MONTHLY CLIMATE DATA #
############################

mc <- read.table("Data/all.sites.CIP.mean.tempANDrainfall.1971_2012.txt",header=T)
mc$site <- as.factor(mc$site)
# climnames <- c("temp","runtemp","rain","runrain") # running means not present
# names(mc)[4:7] <- climnames 
climnames <- c("temp","rain")
names(mc)[4:5] <- climnames 
mmc <- melt(mc,id=c("site","year","month"))
cmc <- dcast(mmc, site + year ~ ...)
monthcols <- paste(rep(climnames,12),rep(1:12,each=length(climnames)),sep="")
names(cmc)[-(1:2)] <- monthcols 
indm <- merge(ind2,cmc,by=c("site","year"))

# table(levels(ind2$site) %in% levels(cmc$site))
# 	sites are missing from climate data, but seem to be filled in incorrectly

##########################
# SUBSET DATA BY SPECIES #
##########################

### Create visit dataset

sitevis <- sumframe(indm)

### Read in species traits

sps <- read.csv("Data/speciestraits208.csv",header=T) # Roger Dennis species data
spi <- read.csv("Data/species_list_with_codes_08Mar2013.csv",header=T) # BMS species index

### Create list of species subsets

	# using spi because more species
splist.full <- levels(spi$english.name)
enoughcounts <- with(indm,tapply(count,species,function(x) sum(x,na.rm=T)))>100
	 # at least 100 counts to include
keepspecies <- levels(indm$species)[!is.na(enoughcounts) & enoughcounts==T]
	# BMS codes of species to keep
matchspeciesnames <- spi$english.name[match(keepspecies,spi$BMScode)]
keepspeciesnames <- na.omit(matchspeciesnames)

### Split by brood

broodlist <- with(indm,tapply(brood,species,table))
broodmat <- matrix(unlist(broodlist),nc=3,byrow=T)#
twobroods <- apply(broodmat,1,function(x) sum(x>0)>1)
names(twobroods) <- names(broodlist)[!sapply(broodlist,is.null)]

# In general, second brood larger, so using that
# but could use brood with the largest sample size

### Create species datasets

splist <- splist.full[splist.full %in% keepspeciesnames]
tblist <- twobroods[names(twobroods) %in% keepspecies[!is.na(matchspeciesnames)]]
#usebrood <- splist[splist.full %in% keepspeciesnames]

spl <- list()
for(s in 1:length(splist)){
	spno <- spi$BMScode[spi$english.name==splist[s]] # species number
	if(!spno %in% indm$species) next
	spd <- droplevels(subset(indm,species==spno))
	spd2 <- dataprep(spd,spno,sitevis,brood=tblist[s])
	#spd2 <- subcounts(spd2,taucounts=2,tauyears=15) # exclude low counts
	spd3 <- subset(spd2,is.finite(growth)) # exclude extinctions?
		# Remove "colonisations" but not "extinctions"
	spl[[s]] <-  droplevels(spd3)
	}
names(spl) <- splist

### select only datasets with a sufficient number of records:

datapresent <- lapply(spl,function(x) nrow(x)>100) 
spl <- spl[unlist(datapresent)]

#splitbrood <- function(mydata){ split(mydata,mydata$brood) }
#splb <- lapply(spl,splitbrood)
#curdat <- splb$"Adonis Blue"
#sapply(splb,length)
#sapply(splb,function(x)lapply(x,nrow))

### Match up to Jeremy Thomas data
thom <- read.csv("Data/roy_thomas_table_namecompatible.csv",header=T,
	stringsAsFactors=F)
thom <- sapply(thom,function(x)replace(x,x=="",NA))
	# blank spaces = NA

# Bivoltine, low number of sites 

###########################################################
# LOG COUNTS BINNED BY DENSITY AND CATEGORISED BY CLIMATE #
###########################################################

# Have changed averaging method to MEDIAN

sem <- function(x) 1.253*sd(x)/sqrt(length(x))
# standard error of the MEDIAN (assumption of normal distribution)
# http://davidmlane.com/hyperstat/A106993.html

curdat <- spl[[20]]
curname <- names(spl)[20]
climname <- "pretemp8"
nclim=3
ndens=10

colvec <- c("blue","orange","red")
pchvec <- rep(16,3)

splotbyclim <- function(curdat,curname,climname,nclim=2,ndens=10,incldat=F,thomas=T,axlabs=F,...){
  
  require("plyr") # for ddply
  require("Hmisc") # for cut2
  require("gplots") # for plotCI
  
  curdat2 <- subset(curdat,!is.na(lpredens) & count>0) # change?  
  curdat2$climcat <- cut2(curdat2[,climname],g=nclim,levels.mean=T)
  # split clim into n=nclim quantiles
  
  curdat2$dencat <- cut2(curdat2$lpredens,g=ndens,levels.mean=T)
    # chopping up densities in the same way for all clim categories

  # sordat <- ddply(curdat2,.(climcat),function(d){
  #   data.frame(d,dencat=cut2(d$lpredens,g=ndens,levels.mean=T))
  #   })
  
  sordat2 <- ddply(curdat2,.(climcat,dencat),function(d){
    c(gmeans=median(log(d$growth)),gses=sem(log(d$growth)),climcat=d$climcat[1])
    })
  # calculate means, sem, and density for each density/temperature class
  
  with(sordat2, 
    plotCI(as.numeric(as.character(dencat)),gmeans,
      ui=gmeans+1.96*gses,
      li=gmeans-1.96*gses,
      pch=pchvec[as.numeric(climcat)],
      col=colvec[as.numeric(climcat)],
      gap=0,
      ylab=ifelse(axlabs==T,"r",""),
      xlab=ifelse(axlabs==T,paste0("ln N (",climname,")"),""),
      bty="l",
      las=1,
      ...
      )
    )
  abline(h=0,lty=3)
  # plot results
  
  if(thomas==T){
    thom.col <- match(climname,colnames(thom))
    thom.row <- match(curname,thom[,colnames(thom)=="english.name"])
    thomtext <- thom[thom.row,thom.col]
    if(!is.na(thomtext)) mtext(thomtext,side=3,line=0,col="blue")
    # find and plot result from Roy & Thomas 2003
    }
  
  if(incldat==T) return(sordat2)
  
  }

splotbyclim(spl[[13]],which(names(spl)=="Brimstone"),"pretemp8",nclim=3,ndens=10)
splotbyclim(spl[[46]],which(names(spl)=="Holly Blue"),"pretemp8",nclim=2,ndens=10)

### Inputs
nclim <- 3
ndens <- 10 # was 20
climnames <- paste(rep(c("temp","rain"),each=12),rep(1:12,2),sep="")
fclimnames <- paste(rep(c("","pre"),times=24),rep(climnames,each=2),sep="")

denbins <- function(curdat,curname){
  sapply(fclimnames,function(x){
    splotbyclim(curdat=curdat,curname=curname,climname=x,nclim=nclim,ndens=ndens,axlabs=T)
    })
  mtext(curname,side=3,cex=2,outer=T)
  }

pdf(paste0("Plots/binnedgrow_density_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=13.5,height=18)
par(mfrow=c(8,6),mar=c(4,4,2,2)+0.1,oma=c(0,0,4,0))
mapply(denbins,spl,names(spl))
dev.off()

pdf(paste0("Plots/binnedgrow_density_pretemp8_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=18,height=18)
par(mfrow=c(8,8),mar=c(4,4,2,2)+0.1,oma=c(0,0,0,0))
for(i in 1:length(spl)){
  splotbyclim(spl[[i]],names(spl)[i],nclim=3,ndens=10,climname="pretemp8",main=names(spl)[i])
  }
dev.off()

pdf(paste0("Outputs/binnedgrow_density_pretemp8_fixedaxes_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=18,height=18)
par(mfrow=c(8,8),mar=c(4,4,2,2)+0.1,oma=c(0,0,0,0))
for(i in 1:length(spl)){
  splotbyclim(spl[[i]],names(spl)[i],climname="pretemp8",main=names(spl)[i],
    xlim=c(-8,-2),ylim=c(-1,1))
  }
dev.off()


### MARIE-CURIE FIGURES

xrange <- c(-8.5,-1.5)
yrange <- c(-1,1)

addxlab <- function(xname){
  mtext(xname,side=1,line=3,las=1,xpd=T,cex=0.9)
  }

addylab <- function(yname){
  mtext(yname,side=2,line=3,las=3,xpd=T,cex=0.9)
  }

ndens <- 20

tiff(filename=paste0("Outputs/mariecurie_data_illustration_",format(Sys.Date(),"%d%b%Y"),".tiff"),
  width=11, height=3.5, units="in", res=600)

par(mfrow=c(1,4),mar=c(2,2,2,2)+0.1,oma=c(3,3,3,3))
# bluedat <- splotbyclim(spl[[20]],names(spl)[20],climname="temp6",incldat=T,thomas=F,ndens=ndens,
#   axlabs=F)
# mtext(expression(bold((a))~italic(Celastrina~argiolus)),side=3,line=0.5,xpd=T,adj=0.05,cex=0.9)

gatdat <- splotbyclim(spl[[13]],names(spl)[13],climname="temp4",incldat=T,thomas=F,ndens=ndens,
  axlabs=F)
mtext(expression(bold((a))~italic(Pyronia~tithonus)),side=3,line=0.5,xpd=T,adj=0.05,cex=0.9)

addylab(expression(log~population~growth~rate~(yr^-1)))
addxlab(expression(log~population~density~(m^2)))

# bluedat$dencon <- as.numeric(as.character(bluedat$dencat))
# bluemod <- glm(gmeans ~ dencon + climcat,data=bluedat)
# lines(bluedat$dencon[1:10],predict(bluemod)[1:10],col="blue")
# lines(bluedat$dencon[11:20],predict(bluemod)[11:20],col="red")

tortdat <- splotbyclim(spl[[45]],names(spl)[45],climname="temp1",incldat=T,thomas=F,ndens=ndens,
  axlabs=F)
mtext(expression(bold((b))~italic(Aglais~urticae)),side=3,line=0.5,xpd=T,adj=0.05,cex=0.9)

addxlab(expression(log~population~density~(m^2)))

paintdat <- splotbyclim(spl[[29]],names(spl)[29],climname="temp2",incldat=T,thomas=F,ndens=ndens,
  axlabs=F)
mtext(expression(bold((c))~italic(Vanessa~cardui)),side=3,line=0.5,xpd=T,adj=0.05,cex=0.9)

addxlab(expression(log~population~density~(m^2)))

plot(1,1,type="n",xlab="",ylab="n",xaxt="n",yaxt="n",bty="n",ann=F)
legend("left",pch=c(16,21),col=c("blue","red"),
  legend=c("cold","warm"),
  title="Temperature",bty="n",cex=1.15)

# brimdat$dencon <- as.numeric(as.character(brimdat$dencat))
# bluemod <- glm(gmeans ~ dencon * climcat,data=brimdat)
# lines(brimdat$dencon[1:10],predict(bluemod)[1:10],col="blue")
# lines(brimdat$dencon[11:20],predict(bluemod)[11:20],col="red")

dev.off()

spi$sci.name[spi$english.name==names(spl)[13]] # temp 4
spi$sci.name[spi$english.name==names(spl)[45]] # temp 1
spi$sci.name[spi$english.name==names(spl)[29]] # temp 2

nrow(spl[13][[1]]) # 7750
nrow(spl[45][[1]]) # 7706
nrow(spl[29][[1]]) # 5033

############################# PREVIOUS CODE ###########################################

#################
# CHECKING DATA #
#################

curdat <- spl$"Silver-spotted Skipper"
olddat <- subset(ind, species==58) # 58 = SSSk
scounts1 <- with(curdat,tapply(count,list(site,year),function(x) mean(x,na.rm=T)))
scounts2a <- with(olddat,tapply(count,list(site,year),function(x) mean(x,na.rm=T)))
scounts2b <- scounts2a[rownames(scounts2a) %in% rownames(scounts1),
	colnames(scounts2a) %in% colnames(scounts1)]
sum(scounts1-scounts2b,na.rm=T) # counts in spl match counts in ind (original data)

yearmeans <- with(cmc, tapply(temp8,year,mean))
yearmeans2 <- with(curdat, tapply(temp8,year,mean))
plot(yearmeans~names(yearmeans),type="l",xlim=c(1980,2010))
lines(yearmeans2~names(yearmeans2),col="red") # temporal climate data seems fine

spatmeans <- with(indm,tapply(log(temp8),list(east,north),mean))
image(x=as.numeric(rownames(spatmeans)),y=as.numeric(colnames(spatmeans)),z=spatmeans)
spatmeans2 <- with(curdat,tapply(log(temp8),list(east,north),mean))
image(x=as.numeric(rownames(spatmeans2)),y=as.numeric(colnames(spatmeans2)),z=spatmeans2)
	# spatial climate also looks OK (but big variation for SSSk)

pdf("Outputs/dynamcheck.pdf",width=5,height=5)
lapply(split(curdat,curdat$site),function(d){
	matplot(d$year,d[,names(d)%in%c("count","precount")],type="b",pch=c(16,17))
	})
dev.off()

tomclim <- read.table("Data/all.sites.CIP.mean.tempANDrainfall.1971_2006_running.3month.means.txt",
	header=T)
ttemps1 <- with(curdat,tapply(count,list(site,year),function(x) mean(x,na.rm=T)))
ttemps2a <- with(olddat,tapply(count,list(site,year),function(x) mean(x,na.rm=T)))
ttemps2b <- ttemps2a[rownames(ttemps2a) %in% rownames(ttemps1),
	colnames(ttemps2a) %in% colnames(ttemps1)]
sum(ttemps1-ttemps2b,na.rm=T) # Temperatures match Tom's temps

########################
# CLIMATE CORRELATIONS #
########################

mcby <- dcast(mmc,year ~ month + variable,mean)

tempby <- mcby[,grep("temp",names(mcby))]
rainby <- mcby[,grep("rain",names(mcby))]
pdf("Outputs/climcors.pdf",width=20,height=20)
pairs.panels(tempby)
pairs.panels(rainby)
dev.off()

#########
# PLOTS #
#########

# Set species number 
i <- 34
cursp <- spl[[i]]
tred <- rgb(1,0,0,alpha=0.5)

# Some general functions
zerline <- function() abline(h=0,col="blue",lty=2)

# Population dynamics
pdf("popdynamics.pdf",width=6,height=5)
for(j in 1:length(spl)){
	popdata <- spl[[j]]
	sitelist <- levels(popdata$site)
	nsites <- length(sitelist)
	plot(log(count)~year,data=popdata,
		type="n",main=paste(names(spl)[j]))
	mycols <- rainbow(nsites)
	for(i in 1:nsites){
		with(popdata,
			lines(log(count[site==sitelist[i]])~year[site==sitelist[i]],
				col=mycols[i])
			)
		} # close i loop
	} # close j loop
dev.off()

# Density-dependence
dplot <- function(curspd,curname){
	plot(log(growth)~log(precount),
		data=curspd,las=1,
		main=curname,
		ylab=expression(log(lambda)),
		xlab=expression(log(N[t-1])),
		pch=16,col=tred)
	zerline()
	}
pdf("densityplots.pdf",width=5,height=5)
mapply(dplot,spl,names(spl))
dev.off()

# Temperature
gtplot <- function(curspd,curtempvar,prelogged=F,...){

	if(prelogged==F) growvar <- log(curspd$growth)
	if(prelogged==T) growvar <- curspd$loggrowth

	plot(growvar~curtempvar,
		data=curspd,
		ylab=expression(log(lambda)),
		xlab="mean temperature (?C)",
		las=1,
		pch=16,
		...)		

	if(nrow(curspd)>3){
		#mod <- glm(growvar~curtempvar,data=curspd)
		#abline(mod,col="red")
		with(curspd,lines(supsmu(curtempvar,growvar,bass=5),
			col="black",lwd=2))
		}

	zerline()

	}

alltemps <- function(curspd,curname){

	alltempcols <- grep("temp",names(curspd))
	runtempcols <- grep("runtemp",names(curspd))
	tempcols <- alltempcols[!alltempcols %in% runtempcols]
	tempnames <- names(curspd)[tempcols]
	months <- format(ISOdatetime(2000,1:12,1,0,0,0),"%b")
	monthnames <- paste(rep(c("cur","pre"),each=12),months)

	thom.col <- match(tempnames,colnames(thom))
	thom.row <- match(curname,thom[,colnames(thom)=="english.name"])
	thomtext <- thom[thom.row,thom.col]

	par(mfrow=c(4,6),oma=c(0,0,5,0))

	for(i in 1:length(tempcols)){
		curtempvar <- curspd[,tempcols[i]]
		gtplot(curspd, curtempvar, 
			col=as.numeric(curspd$site),main=monthnames[i])
		if(!is.na(thomtext[i])) mtext(thomtext[i],side=3,line=0)
		}

	mtext(curname,side=3,line=1,outer=T,cex=1.75)

	}

pdf("Outputs/alltemps.pdf",width=16,height=15)
mapply(alltemps,spl,names(spl))
dev.off()

### Space Vs Time

spacetime <- function(curspd,curname){

	alltempcols <- grep("temp",names(curspd))
	runtempcols <- grep("runtemp",names(curspd))
	tempcols <- alltempcols[!alltempcols %in% runtempcols]
	tempnames <- names(curspd)[tempcols]
	months <- format(ISOdatetime(2000,1:12,1,0,0,0),"%b")
	monthnames <- paste(rep(c("cur","pre"),each=12),months)
	
	curspd$loggrowth <- log(curspd$growth)
	
	mcursp <- melt(curspd,measure=c(tempnames,"loggrowth"))
	yeardat <- dcast(mcursp,year ~ variable,mean)
	sitedat <- dcast(mcursp,site ~ variable,mean)
	allgrowth <- c(yeardat$loggrowth,sitedat$loggrowth)
	ylim.both <- c(min(allgrowth),max(allgrowth))
	newtempcols <- which(names(yeardat) %in% tempnames)

	thom.col <- match(tempnames,colnames(thom))
	thom.row <- match(curname,thom[,colnames(thom)=="english.name"])
	thomtext <- thom[thom.row,thom.col]

	aveplot <- function(curspd,curtempvar,...){
		gtplot(curspd,curtempvar,prelogged=T,
			ylim=ylim.both,col="gray",...)
		if(!is.na(thomtext[i])) mtext(thomtext[i],side=3,line=0)
		}

	for(i in 1:length(tempcols)){
		aveplot(yeardat, yeardat[,newtempcols[i]], main=paste("time",monthnames[i]))
		aveplot(sitedat, sitedat[,newtempcols[i]], main=paste("space",monthnames[i]))
		}
	
	mtext(curname,side=3,line=1,outer=T,cex=1.75)

	}

pdf("space_vs_time.pdf",width=6,height=70)
par(mfrow=c(24,2),oma=c(0,0,5,0))
mapply(spacetime,spl,names(spl))
dev.off()

###########################
# SPACE VS TIME - DENSITY #
###########################

# Temperature
dplot <- function(curspd,curtempvar,...){

	plot(log(dens)~curtempvar,
		data=curspd,
		ylab=expression(log(bar(Nm^-1))),
		xlab="mean temperature (?C)",
		las=1,
		pch=16,
		...)		

	if(nrow(curspd)>3){
		mod <- glm(log(dens)~curtempvar,data=curspd)
		abline(mod,col="red")
		with(curspd,lines(supsmu(curtempvar,log(dens),bass=5),
			col="black",lwd=2))
		}

	}

dplot.spacetime <- function(curspd,curname){

	alltempcols <- grep("temp",names(curspd))
	runtempcols <- grep("runtemp",names(curspd))
	tempcols <- alltempcols[!alltempcols %in% runtempcols]
	tempnames <- names(curspd)[tempcols]
	months <- format(ISOdatetime(2000,1:12,1,0,0,0),"%b")
	monthnames <- paste(rep(c("cur","pre"),each=12),months)
		
	mcursp <- melt(curspd,measure=c(tempnames,"dens"))
	yeardat <- dcast(mcursp,year ~ variable,mean)
	sitedat <- dcast(mcursp,site ~ variable,mean)
	alldens <- c(log(yeardat$dens),log(sitedat$dens))
	ylim.both <- c(min(alldens[is.finite(alldens)],na.rm=T),
		max(alldens[is.finite(alldens)],na.rm=T))
	newtempcols <- which(names(yeardat) %in% tempnames)

	thom.col <- match(tempnames,colnames(thom))
	thom.row <- match(curname,thom[,colnames(thom)=="english.name"])
	thomtext <- thom[thom.row,thom.col]

	aveplot <- function(curspd,curtempvar,...){
		dplot(curspd,curtempvar,
			ylim=ylim.both,col="gray",...)
		if(!is.na(thomtext[i])) mtext(thomtext[i],side=3,line=0)
		}

	for(i in 1:length(tempcols)){
		aveplot(yeardat, yeardat[,newtempcols[i]], main=paste("time",monthnames[i]))
		aveplot(sitedat, sitedat[,newtempcols[i]], main=paste("space",monthnames[i]))
		}
	
	mtext(curname,side=3,line=1,outer=T,cex=1.75)

	}

pdf("space_vs_time_dens.pdf",width=6,height=70)
par(mfrow=c(24,2),oma=c(0,0,5,0),mar=c(5,5,4,2)+0.1)
mapply(dplot.spacetime,spl,names(spl))
dev.off()

### FlexParamCurve
#nlsList(log(count/precount)~SSposnegRichards(temp8),data=spd3)

###########################################################
# LOG COUNTS BINNED BY CLIMATE AND CATEGORISED BY DENSITY #
###########################################################

sem <- function(x) sd(x)/sqrt(length(x))
	# standard error of the mean

splotbyden <- function(curdat,curname,climname,nde=2,ncl=10){

	require("plyr") # for ddply
	require("Hmisc") # for cut2
	require("gplots") # for plotCI

	curdat2 <- subset(curdat,!is.na(lpredens) & count>0) # change?  
	curdat2$dencat <- cut2(curdat2$lpredens,g=nde,levels.mean=T)
		# split density into n=nde quantiles

	sordat <- ddply(curdat2,"dencat",function(d){
		data.frame(d,climcat=cut2(d[,names(d)==climname],g=ncl,levels.mean=T))
		})

	sordat2 <- ddply(sordat,c("dencat","climcat"),function(d){
		c(gmeans=mean(log(d$growth)),gses=sem(log(d$growth)),dencat=d$dencat[1])
		})
		# calculate means, sem, and density for each density/temperature class

	with(sordat2, 
		plotCI(as.numeric(as.character(climcat)),gmeans,
			ui=gmeans+1.96*gses,
			li=gmeans-1.96*gses,
			pch=16,
			col=dencat,
			gap=0,
			ylab=expression(log(lambda)),
			xlab=climname,
			bty="n",
			las=1,
			)
		)
	# plot results

	thom.col <- match(climname,colnames(thom))
	thom.row <- match(curname,thom[,colnames(thom)=="english.name"])
	thomtext <- thom[thom.row,thom.col]
	if(!is.na(thomtext)) mtext(thomtext,side=3,line=0,col="blue")
	# find and plot result from Roy & Thomas 2003

	}

### Inputs
nde <- 2
ncl <- 10
climnames <- paste(rep(c("temp","rain"),each=12),rep(1:12,2),sep="")
fclimnames <- paste(rep(c("","pre"),times=24),rep(climnames,each=2),sep="")

speciesbins <- function(curdat,curname){
	sapply(fclimnames,function(x){
		splotbyden(curdat=curdat,curname=curname,climname=x)
		})
	mtext(curname,side=3,cex=2,outer=T)
	}

pdf("Outputs/binnedgrow.pdf",width=13.5,height=18)
par(mfrow=c(8,6),mar=c(4,4,2,2)+0.1,oma=c(0,0,4,0))
mapply(speciesbins,spl,names(spl))
dev.off()

########################
# SINGLE FLIGHT PERIOD #
########################

### Checking brood numbers by species
broodno <- with(fullind, tapply(brood,species,table))
broodpres <- lapply(broodno,function(x) x>0)
broodmat <- matrix(unlist(broodpres),
	ncol=3, byrow=T,
	dimnames=list(levels(ind$species),0:2))
rownames(broodmat) <- spi$english.name[match(as.numeric(rownames(broodmat)),spi$BMScode)]
broodmat2 <- data.frame(broodmat[!is.na(rownames(broodmat)),])
names(broodmat2) <- c("none","one","two")
sapply(broodmat2,sum) # approx 15 species have two broods

onebroodnames <- sort(rownames(subset(broodmat2,one==F & two==F)))
# but:
# Comma has overwintering
# large white may only have one peak
# Orange tipe looks like several peaks
# painted lady has two peaks
# small blue might have one
# small heath might have two
# small white may only have one
# speckled wood might have several

###############################
# VARIANCE COMPONENT ANALYSES #
###############################

partvar <- function(model){
	varcomps <- c(unlist(VarCorr(model)), resid=attributes(VarCorr(model))$sc^2)
	sapply(varcomps, function(x) round(x / sum(varcomps), 2))
	}

varcomps <- function(curdat){
	if ( length(unique(curdat$year))*length(levels(curdat$site)) > nrow(curdat)
		& length(unique(curdat$year))>5 & length(levels(curdat$site))>5){
		cmod <- lmer(count ~ offset(log(tlen)) + (1|site) + (1|year),
			family="gaussian",data=curdat)
		partvar(cmod)
		}
	}

spvarcomps <- sapply(spl, varcomps)
spvarmat <- matrix(unlist(spvarcomps),nc=3,byrow=T)
spvartab <- as.table(spvarmat)
colnames(spvartab) <- names(spvarcomps[[2]])
rownames(spvartab) <- names(spvarcomps)[!sapply(spvarcomps,is.null)]

pdf("sizevars_withoffset.pdf",width=15,height=8)
par(mar=c(12,4,4,2)+0.1)
barplot(t(spvartab),beside=T,las=3)
dev.off()

pdf("sizevars_withoffset_grouped.pdf",width=12,height=6)
par(mar=c(5,5,4,2)+0.1)
barplot(spvartab,beside=T,las=1,
	cex.axis=1.5,cex.names=1.5,
	cex.lab=1.5,ylab="variance explained (%)")
dev.off()

tempvars <- sapply(split(mc,mc$month), function(d){
	tmod <- lmer(temp ~ 1 + (1|site) + (1|year), data=d)
	partvar(tmod)
	})

rainvars <- sapply(split(mc,mc$month), function(d){
	tmod <- lmer(rain ~ 1 + (1|site) + (1|year), data=d)
	partvar(tmod)
	})

pdf("tempvars.pdf",width=8,height=7)
barplot(tempvars,beside=T,las=1,ylim=c(0,1),
	ylab="proportion of variance explained",xlab="month",
	main="mean temperature",
	legend.text=rownames(tempvars),
	args.legend=list("topright",bty="n")
	)
dev.off()

pdf("rainvars.pdf",width=8,height=7)
barplot(rainvars,beside=T,las=1,ylim=c(0,1),
	ylab="proportion of variance explained",xlab="month",
	main="total rainfall",
	legend.text=rownames(tempvars),
	args.legend=list("topright",bty="n")
	)
dev.off()

####################
# SPATIAL PATTERNS #
####################

# See code for weekly if want to do separately

corplot <- function(curdat,curname){

	nsites <- length(levels(curdat$site))
	cormat <- matrix(nc=nsites,nr=nsites)
	usite <- curdat[match(levels(curdat$site),curdat$site),
		names(curdat) %in% c("east","north")]
	rownames(usite) <- levels(curdat$site)

	distmat <- dist(usite)
	growbysite <- with(curdat, tapply(growth,list(year,site),print))

	if(is.numeric(distmat) & is.numeric(growbysite)){
		cormat <- as.dist(cor(growbysite,use="pairwise.complete.obs")) 
		plot(distmat,cormat,main=curname)
		lines(supsmu(distmat,cormat),col="red")
		}

	}
	
pdf("Outputs/corplots.pdf",width=7,height=7)
mapply(corplot,spl,names(spl))
dev.off()

############
# ANALYSES #
############

### SSSk models

curdat <- spl$"Silver-spotted Skipper"
require(lme4)
require(mgcv)

curdat$tpretemp8 <- with(curdat, pretemp8-mean(pretemp8))
curdat$tlpredens <- with(curdat, lpredens - mean(lpredens,na.rm=T))

mod1a <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens) 
	+ (1|site),
	family="poisson",data=curdat)

mod1b <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens) 
	+ (1|year),
	family="poisson",data=curdat)

mod1c <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens) 
	+ (1|site) + (1|year),
	family="poisson",data=curdat)

mod1d <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens)
	+ lpredens 
	+ (1|site) + (1|year),
	family="poisson",data=curdat)

mod1e <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens)
	+ lpredens 
	+ (1+lpredens|site) + (1|year),
	family="poisson",data=curdat)

mod1f <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens)
	+ lpredens 
	+ (1|site) + (1+lpredens|year),
	family="poisson",data=curdat)

mod1g <- glmer(count ~ offset(log(tlen))
	+ offset(lpredens) 
	+ lpredens 
	+ (1+lpredens|site) + (1+lpredens|year),
	family="poisson",data=curdat)

aicm <- function(l) sapply(l,function(x) AIC(logLik(x)) ) 

aicm(list(siteint=mod1a,yearint=mod1b,siteandyearint=mod1c,
	dd=mod1d, ddsite=mod1e,ddyear=mod1f,ddyearandsite=mod1g))
	# clearly density-dependence going on
	# density-dependence varies among sites and years

### With climate

mod2a <- glmer(count ~ offset(log(tlen))
 	+ offset(lpredens) 
	+ lpredens 
	+ temp8
	+ (1+lpredens|site),
	family="poisson",data=curdat)

	# random variation by year and site causes singularities here

mod2b <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens) 
	+ lpredens 
	+ temp8
	+ lpredens:temp8
	+ (1+lpredens|site),
	family="poisson",data=curdat)

mod2c <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens) 
	+ lpredens 
	+ pretemp8
	+ (1+lpredens|site),
	family="poisson",data=curdat)

mod2d <- glmer(count ~ offset(log(tlen)) 
	+ offset(tlpredens) 
	+ tlpredens 
	+ tpretemp8
	+ tlpredens:tpretemp8
	+ (1+tlpredens|site),
	family="poisson",data=curdat)

mod2e <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens) 
	+ lpredens 
	+ temp8
	+ pretemp8
	+ lpredens:pretemp8
	+ (1+lpredens|site),
	family="poisson",data=curdat)
	
	# best model, but weird interaction

mod2g <- glmer(count ~ offset(log(tlen)) 
	+ offset(lpredens) 
	+ lpredens 
	+ temp8
	+ pretemp8
	+ (1+lpredens|site),
	family="poisson",data=curdat)

aicm(list(curad=mod2a,curint=mod2b,prevad=mod2c,
	prevint=mod2d,adprevint=mod2e,adcurint=mod2f,
	alladd=mod2g))

### Quick GAMM

curdat$obslevel <- factor(1:nrow(curdat))
intgam <- gam(count ~ 
	te(tpretemp8,tlpredens)
	# + s(east,north,k=25)
	# + s(gridno,bs="re"),
	+ s(site,bs="re"),
	offset=lprecount,
	data=curdat,
	family=poisson,
	scale=-1	 
	)

sepgam <- gam(count ~ 
	s(tpretemp8)
	+ tlpredens
	+ s(site,bs="re"),
	offset=lprecount,
	data=curdat,
	family=poisson,
	scale=-1	 
	)

# "scale" forces the scale parameter of the Poisson to be treated 
# as unknown, and smoothing parameters to be estimated by
# GCV, rather than UBRE, which is the Poisson default: 
# hence the model is employing an ?overdispersed Poisson? structure

gam.check(intgam)
vis.gam(intgam,theta=135,phi=30)
vis.gam(intgam,theta=-45,phi=30,view=c("east","north"))

vis.gam(sepgam,theta=135,phi=30)

############################
# DENSITY-DEPENDENCE PLOTS #
############################

curdat <- spl$"Silver-spotted Skipper"

### Same year

n.seq <- 100
temp8seq <- with(curdat, seq(min(temp8),max(temp8),length.out=n.seq))
lpdseq <- with(curdat, seq(min(lpredens[is.finite(lpredens)]),max(lpredens[is.finite(lpredens)]),length.out=n.seq))
temp8seq2 <- rep(temp8seq,each=n.seq)
lpdseq2 <- rep(lpdseq,times=n.seq)

params <- fixef(mod2b)

mm1A <- model.matrix(~lpdseq*rep(mean(temp8seq),n.seq))
plambda1A <- matrix(mm1A %*% params,nc=1)
mm1B <- model.matrix(~rep(mean(lpdseq),n.seq)*temp8seq)
plambda1B <- matrix(mm1B %*% params,nc=1)

par(mfrow=c(1,2))
plot(plambda1A ~ exp(lpdseq),type="l")
with(curdat, points(log(growth)~exp(lpredens)))
plot(plambda1B ~ temp8seq,type="l")
with(curdat, points(log(growth)~temp8))

mm2 <- model.matrix(~lpdseq2*temp8seq2)
plambda2 <- matrix(mm2 %*% params,nc=n.seq)
persp(exp(lpdseq),temp8seq,plambda2,
	theta=45,phi=15,
	xlab=expression(Density (nm^-2)),
	ylab="Mean Temperature",
	zlab="log(lambda)",
	col="red",ticktype="detailed")

require("rgl")
plot3d(temp8seq2,exp(lpdseq2),plambda2)

### Previous year

curdat_mod <- subset(curdat,is.finite(tlpredens))

n.seq <- 100
tpretemp8seq <- with(curdat_mod, seq(min(tpretemp8),max(tpretemp8),length.out=n.seq))
tlpdseq <- with(curdat_mod, seq(min(tlpredens),max(tlpredens),length.out=n.seq))
tpretemp8seq2 <- rep(tpretemp8seq,each=n.seq)
tlpdseq2 <- rep(tlpdseq,times=n.seq)

params <- fixef(mod2d)

mm1A <- model.matrix(~tlpdseq*rep(mean(tpretemp8seq),n.seq))
plambda1A <- matrix(mm1A %*% params,nc=1)
mm1B <- model.matrix(~rep(mean(tlpdseq),n.seq)*tpretemp8seq)
plambda1B <- matrix(mm1B %*% params,nc=1)

par(mfrow=c(1,2))
plot(plambda1A ~ exp(tlpdseq),type="l")
with(curdat, points(log(growth)~exp(tlpredens)))
plot(plambda1B ~ tpretemp8seq,type="l")
with(curdat, points(log(growth)~tpretemp8))

mm2 <- model.matrix(~tlpdseq2*tpretemp8seq2)
plambda2 <- matrix(mm2 %*% params,nc=n.seq)
par(mfrow=c(1,1))
persp(tlpdseq,tpretemp8seq,plambda2,
	theta=45,phi=15,
	xlab=expression(Density (nm^-2)),
	ylab="Mean Temperature",
	zlab="log(lambda)",
	col="red",ticktype="detailed")

require("rgl")
plot3d(tpretemp8seq2,exp(tlpdseq2),plambda2)

### Both (three vars)

n.seq <- 100
pretemp8seq <- with(curdat, seq(min(pretemp8),max(pretemp8),length.out=n.seq))
lpdseq <- with(curdat, seq(min(lpredens),max(lpredens),length.out=n.seq))
pretemp8seq2 <- rep(pretemp8seq,each=n.seq)
lpdseq2 <- rep(lpdseq,times=n.seq)
temp8seq2 <- rep(mean(curdat$temp8),n.seq^2)

params <- fixef(mod2e)

mm3 <- model.matrix(~lpdseq2+temp8seq2+pretemp8seq2+lpdseq2:pretemp8seq2)
plambda3 <- matrix(mm3 %*% params,nc=n.seq)

pdf("sss_persp_pretemp8.pdf",width=7,height=7)
persp(lpdseq,pretemp8seq,plambda3,
	theta=45,phi=0,
	xlab="log(N/m)",
	ylab="Mean Temperature (?C)",
	zlab="log(lambda)",
	col="red",ticktype="detailed")
dev.off()

nlines <- 5
sel <- seq(0,n.seq,length.out=nlines)

pdf("sss_dd_pretemp8_lines.pdf",width=3.5,height=3.5)
layout(matrix(rep(c(1,2),c(20,5)),nc=5,nr=5,byrow=T))
par(mar=c(4,5,4,2)+0.1)
matplot(lpdseq,plambda3[,sel],
	type="l",
	ylab=expression(log(lambda)),
	xlab=expression(log(Nm^-1)),
	lty=rep(1,nlines),
	las=1,
	bty="n"
	)
legend("topright",legend=round(pretemp8seq[sel],1),
	col=1:nlines,lty=rep(1,nlines),
	bty="n",title="Temperature (?C)")
abline(h=0,lty=3,col="orange")
par(mar=c(0,5,0,2)+0.1,bty="n")
boxplot(curdat$lpredens,horizontal=T,range=0,col="grey",xaxt="n")
dev.off()

pdf("sss_pretemp8_byd_lines.pdf",width=3.5,height=3.5)
layout(matrix(rep(c(1,2),c(20,5)),nc=5,nr=5,byrow=T))
par(mar=c(4,5,4,2)+0.1)
matplot(pretemp8seq,t(plambda3[sel,]),
	type="l",
	ylab=expression(log(lambda)),
	xlab="Temp",
	lty=rep(1,nlines),
	las=1,
	bty="n"
	)
legend("topleft",legend=round(lpdseq[sel],1),
	col=1:nlines,lty=rep(1,nlines),
	bty="n",title=expression(log(Nm^-1)) )
abline(h=0,lty=3,col="orange")
par(mar=c(0,5,0,2)+0.1,bty="n")
boxplot(curdat$pretemp8,horizontal=T,range=0,col="grey",xaxt="n")
dev.off() # NB: hasn't been checked carefully

########
# SIMS #
########

curdat_sim <- subset(curdat,is.finite(tlpredens))

pfun <- function(m,x,z,s){
	predict(m, 
		newdata=data.frame(tlpredens=x,tpretemp8=z,tlen=10^3,site=s),
		type="link",
		re.form=NULL
		)
	}

gfun <- function(m,x,z,s){
	predict.gam(m, 
		newdata=data.frame(tlpredens=x,tpretemp8=z,tlen=10^3,site=s),
		type="link",
		re.form=NULL
		)
	}
	# unfinished

n.seq <- 100
tpretemp8sim0 <- rnorm(n.seq,mean=0,sd=sd(curdat_sim$tpretemp8))
tpretemp8sim1 <- tpretemp8sim0 + sd(curdat_sim$tpretemp8)
tpretemp8sim2 <- tpretemp8sim0 - sd(curdat_sim$tpretemp8)

tlpd1a <- tlpd2a <- tlpd1b <- tlpd2b <- vector(mode="numeric",length=n.seq)
tlpd1a[1] <- tlpd2a[1] <- tlpd1b[1] <- tlpd2b[1] <- 3.75 	# start at density on log scale

for(t in 2:n.seq){	
	tlpd1a[t] <- pfun(m=mod2d,x=tlpd1[t-1],z=tpretemp8sim1[t],s=1309)
	tlpd2a[t] <- pfun(m=mod2d,x=tlpd2[t-1],z=tpretemp8sim2[t],s=1309)
	tlpd1b[t] <- pfun(m=mod2d,x=tlpd1[t-1],z=tpretemp8sim1[t],s=1312)
	tlpd2b[t] <- pfun(m=mod2d,x=tlpd2[t-1],z=tpretemp8sim2[t],s=1312)
	}

csim <- cbind(tlpd1a,tlpd2a,tlpd1b,tlpd2b)
matplot(csim,type="l",lty=1)
mean(tlpd1a) - mean(tlpd2a)
mean(tlpd1b) - mean(tlpd2b)
sd(tlpd1a) - sd(tlpd2a)
sd(tlpd1b) - sd(tlpd2b)

plot(density(tlpd1a),xlim=range(csim),col="black",lty=1)
lines(density(tlpd2a),col="black",lty=2)
lines(density(tlpd1b),col="red",lty=1)
lines(density(tlpd2b),col="red",lty=2)

###########################
# AUTOMATED MODEL-FITTING #
###########################

monlist <- read.csv("monlist.csv",header=T)
names(monlist)[names(monlist)=="X"] <- "species"

# Set up data
curdat <- splb$"Silver-spotted Skipper"
curname <- "Silver-spotted Skipper"
curb <- curdat$"0"

# Set x variable names to keep
fl <- subset(monlist,species==curname)
weathname <- c("temp","rain") # ignores running means
preweathname <- c("pretemp","prerain")

pastenosep <- function(x,y){ paste(x,y,sep="") }
prekeep <- as.vector(outer(preweathname,(fl$first):12,pastenosep))
	# exclude vars before flight period in prev 
postkeep <- as.vector(outer(weathname,1:(fl$first-1),pastenosep))
	# exclude vars at or during flight period in current
keepnames <- c(prekeep,postkeep)

# Fit models 
require(lme4)
fullmod <- lmer(count~log(precount)*pretemp8+(1+log(precount)|gridno),
	data=curbroo,family="poisson")

modtab <- dredge(fullmod)

#################
# GAMM ANALYSES #
#################

###################################################################
# A note on offsets: 
#
# "gam" requires that offset must be specified with "offset=", 
# otherwise predict.gam gives incorrect values.
#
# "gamm" requires offset to specified as model term:
# + offset(variable)
# AND the variable must be log-transformed prior to inclusion
###################################################################

# fit model additively (without interaction)
# do for all pairs and rank by AIC
# additional issues: overdispersion, synchrony

# NB: offset must be reparameterised in combined form! 
# (or tlen is missed)

# NB: log(predens)+log(tlen)=log(predens*tlen)=log(precount)
# So created new variable, "lprecount", to ease fitting

require(mgcv)
require(MASS)
curdat <- spl$"Gatekeeper"
curdat$tpretemp8 <- with(curdat, pretemp8-mean(pretemp8))
curdat$tlpredens <- with(curdat, lpredens - mean(lpredens,na.rm=T))

par(mfrow=c(2,2))
plot(log(growth)~lpredens,data=curdat)
plot(log(growth)~temp8,data=curdat)
with(curdat, lines(supsmu(temp8,log(growth)),col="red"))
plot(log(growth)~pretemp8,data=curdat)
with(curdat, lines(supsmu(pretemp8,log(growth)),col="red"))

### Estimate theta for overdispersion
nbm <- glm.nb(count ~ offset(lprecount) 
	+ pretemp8*lpredens
	+ gridno,
	link=log,
	data=curdat)
thetahat <- nbm$theta

### Model
intglm <- gam(count ~ 
	pretemp8*lpredens
	+ s(east,north),
	#+ s(gridno,bs="re"),
	offset=lprecount,
	family=negbin(theta=thetahat),
	data=curdat)

system.time(
intgam <- gam(count ~ 
	te(pretemp8,lpredens)
	+ s(east,north,k=25)
	+ s(gridno,bs="re"),
	offset=lprecount,
	#gamma=1.4,
	family=negbin(theta=thetahat),
	data=curdat)
	)
# 40 minutes with negative binom and full spatial grid
# 25 minutes with 10 x 10 grid (100 knots)
# 40 minutes with 5 x 5 grid (25 knots) and gamma=1.4
# 44 minutes with 5 x 5 grid (25 knots) and gamma=1 (!)

addgamm <- gamm(count ~ 
	s(pretemp8) + s(lpredens)
	+ s(east,north)
	+ offset(lprecount),
	#random=list(gridno=~1),
	#correlation=corExp(form=~year|gridno), # does year work? 
	family=negbin(theta=thetahat), 
	data=curdat)

fixgam <- gam(count ~ 
	s(tpretemp8)
	+ tlpredens
	+ s(site,bs="re"),
	offset=lprecount,
	data=curdat,
	family=poisson,
	scale=-1	 
	)

system.time(
intgamm <- gamm(count ~ 
	te(pretemp8,lpredens)	#k=...
	+ s(temp8)			#k=...
	+ s(east,north,k=100)
	+ offset(lprecount),
	random=list(gridno=~1), 
	correlation=corExp(form=~year), # SORT OUT AUTOCORRELATION
	family=negbin(theta=thetahat), 
	gamma=1.4,
	niterPQL=100,
	data=curdat)
	)	

model1$modelStruct$corStruct
plot(ACF(model1),alpha=0.05)
plot(ACF(model1,resType="normalized"),alpha=0.05)
	# taking correlation into account
	
	# convergence failure even after 40 iterations when unrestricted
	# convergence after 9 iterations with knots=10,4, but v. slow 
	# including exp autocorrelation, convergence after 1 hour. 
	# may need to use gc()
	
# Accounting for overdispersion greatly changes effect of temperature

vis.gam(intglm,view=c("pretemp8","lpredens"),theta=45,ticktype="detailed")
vis.gam(intglm,view=c("east","north"),theta=45,ticktype="detailed")
gam.check(intglm)

vis.gam(intgam,view=c("pretemp8","lpredens"),theta=45,ticktype="detailed")
vis.gam(intgam,view=c("east","north"),theta=45,ticktype="detailed")
gam.check(intgam)

vis.gam(intgamm$gam,view=c("pretemp8","lpredens"),theta=45)#se=1.96
vis.gam(intgamm$gam,view=c("east","north"),theta=-45)
vis.gam(intgamm$gam,view=c("lprecount","temp8"),theta=-45)
gam.check(intgamm$gam)

vis.gam(addgamm$gam,view=c("pretemp8","lpredens"),theta=30)#se=1.96
vis.gam(addgamm$gam,view=c("east","north"),theta=-45)
gam.check(addgamm$gam)

### Frankfurt work

bigsites <- with(curdat, names(which(tapply(lprecount,site,length)>1)))

scurdat <- subset(curdat,site %in% bigsites)
scurdat$obslevel <- factor(1:nrow(scurdat)) # can't figure out how to implement this
scurdat$gridno <- as.factor(as.character(scurdat$gridno))
scurdat$site <- as.factor(as.character(scurdat$site))

system.time(
  intgam2 <- gam(count ~ 
                 te(pretemp8,lpredens)
                 + s(site,bs="re"),
                 offset=lprecount,
                 family=quasipoisson(link=log),
                 data=scurdat
  )
)

vis.gam(intgam2,view=c("pretemp8","lpredens"),theta=135,ticktype="detailed")
vis.gam(intgam2,view=c("pretemp8","lpredens"),theta=270,ticktype="detailed")
summary(intgam2)
gam.check(intgam2)

### Predictions

# Choose inputs
npred <- 100
premod <- intgam2 # intgam
climvar <- curdat$pretemp8
lpdvar <- curdat$lpredens
climname <- "pretemp8"

# Make predictions
climseq <- seq(min(climvar),max(climvar),length.out=npred)
lpdseq <- seq(min(lpdvar,na.rm=T),max(lpdvar,na.rm=T),length.out=npred)
pframe <- expand.grid(lpredens=lpdseq,climvar=climseq)
# tlenseq <- rep(1,npred^2)
pframe$site = names(sort(table(scurdat$site), decreasing=T)[1])
	# east=mean(curdat$east),north=mean(curdat$north),
	# tlen=tlenseq)
names(pframe)[names(pframe)=="climvar"] <- climname

preds <- predict.gam(premod,pframe,type="link",se.fit=T)
	# using untransformed predictions because want log(lambda)
predmat <- matrix(preds$fit,nc=npred,byrow=F)
rownames(predmat) <- signif(lpdseq,3)
colnames(predmat) <- signif(climseq,3)

# 2D plot (could make this into function)
require(fields)
nlines <- 10
sel <- round(seq(1,npred,length.out=nlines),0)
layout(matrix(rep(c(1,2),c(20,5)),nc=5,nr=5,byrow=T))
par(mar=c(4,5,4,2)+0.1)
matplot(lpdseq,predmat[,sel],
	type="l",
	ylab=expression(log(lambda)),
	xlab=expression(log(Nm^-1)),
	col=tim.colors(nlines),
	lty=rep(1,nlines),
	las=1,
	bty="n"
	)
legend("topright",legend=round(climseq[sel],1),
	col=tim.colors(nlines),lty=rep(1,nlines),
	bty="n",title="Temperature (?C)")
abline(h=0,lty=3,col="orange")
par(mar=c(0,5,0,2)+0.1,bty="n")
boxplot(curdat$lpredens,horizontal=T,range=0,col="grey",xaxt="n")

layout(matrix(rep(c(1,2),c(20,5)),nc=5,nr=5,byrow=T))
par(mar=c(4,5,4,2)+0.1)
matplot(climseq,t(predmat[sel,]),
        type="l",
        ylab=expression(log(lambda)),
        xlab="Temperature (C)",
        col=tim.colors(nlines),
        lty=rep(1,nlines),
        las=1,
        bty="n"
)
legend("topright",legend=round(lpdseq[sel],1),
       col=tim.colors(nlines),lty=rep(1,nlines),
       bty="n",title=expression(log(Nm^-1)))
abline(h=0,lty=3,col="orange")
par(mar=c(0,5,0,2)+0.1,bty="n")
boxplot(curdat$pretemp8,horizontal=T,range=0,col="grey",xaxt="n")


### Spatial patterns

# Choose inputs
evar <- curdat$east
nvar <- curdat$north

# Make predictions
eseq <- seq(min(climvar),max(climvar),length.out=npred)
nseq <- seq(min(lpdvar,na.rm=T),max(lpdvar,na.rm=T),length.out=npred)
eseq2 <- rep(climseq,each=npred)
nseq2 <- rep(lpdseq,times=npred)

pframe <- data.frame(climvar=rep(mean(climvar),npred^2),
	lpredens=rep(mean(curdat$lpredens,na.rm=T),npred^2),
	east=eseq2,north=nseq2,
	tlen=tlenseq)
names(pframe)[names(pframe)=="climvar"] <- climname

preds <- predict.gam(premod,pframe,type="link",se.fit=T)
	# using untransformed predictions because want log(lambda)
predmat <- matrix(preds$fit,nc=npred,byrow=F)

# 3D plot
persp(eseq,nseq,predmat,theta=45,phi=30) # double-check that axis vars are right

#########################
# DETECTABILITY EXAMPLE #
#########################

pdf("detection_example.pdf",width=5,height=5)
N <- 20
t <- 20
tseq <- 1:t
Nseq <- rep(N,t)
set.seed(1000)
Cseq <- rbinom(t,N,0.1)
plot(Cseq ~ tseq,type="b",col="red",ylim=c(0,max(c(Cseq,Nseq))),
	pch=16,ylab=expression(N[t]),xlab="year")
lines(Nseq ~ tseq,type="b",col="blue",pch=16)
legend("right",legend=c("true","observed"),
	col=c("blue","red"),pch=rep(16,2),lty=rep(1,2),
	bty="n")
dev.off()

#######################
# FUNCTIONAL RESPONSE #
#######################

richards <- function(x,k1,k2,k3,k4){
	k1/(1+(k1/k2-1)*exp(-k3*k4*x))^(1/k4)
	}

curve(richards((x-15),10,15,1,1),xlim=c(0,30),
	xlab="Temperature (?C)",
	ylab=expression(log(lambda)),
	las=1)

####################
# CLIMATE ANALYSES #
####################

pdf("yearbysite_clims_sssk.pdf",width=10,height=10)
with(curdat, plot(pretemp8,temp8,type="n"))
with(curdat,text(pretemp8,temp8,labels=paste(year),col=as.numeric(site)))
dev.off()

with(curdat, plot(
	tapply(temp8,site,mean),
	tapply(temp8,site,function(x) sd(x))
	))

# make map
with(site2,plot(east,north))

