#################################
# BMS DATA PROCESSING FUNCTIONS #
#################################

#########################################################
# addnas												
#
# adds na values to missing site/year/week combinations
# NB: note that brood is removed
#########################################################

#data=bms
#spref="Silver-spotted Skipper"
#data=ind
#spref=54

# Sumframe = summing counts for all species 

sumframe <- function(mydata){
	
	require(data.table)

	mainnames <- c("species","site","year","count")

	if("week" %in% names(mydata)){
		bms <- mydata
		mbms <- melt(bms[,names(bms) %in% c(mainnames,"week")],measure="count")
		sumbms <- dcast(mbms,site+year+week~variable,sum)
		names(sumbms)[names(sumbms)=="count"] <- "sumcount"

		SDT <- as.data.table(sumbms)
		setkey(SDT,site,year,week)
		SDT2 <- SDT[CJ( 
			unique(site),
			seq(min(year),max(year)),
			seq(min(week),max(week))
			)]
		}

	else{
		bms <- mydata
		mbms <- melt(bms[,names(bms) %in% mainnames],measure="count")
		sumbms <- dcast(mbms,site+year~variable,sum)
		names(sumbms)[names(sumbms)=="count"] <- "sumcount"

		SDT <- as.data.table(sumbms)
		setkey(SDT,site,year)
		SDT2 <- SDT[CJ( 
			unique(site),
			seq(min(year),max(year))
			)]
		}

	return(SDT2)

	}

# SDT2 is a summed data.table object (use "sumframe" to create)
# SDT2 = sitevis 
# brood = T
# data=droplevels(subset(indm,species==54))

dataprep <- function(data,spref,SDT2,brood=0){

	require(data.table)
	spdat <- droplevels(subset(data,species==spref))
	exclude <- c("species","gref") # removed tlen
	
	matchnames <- ! names(spdat) %in% c(exclude) # removed week
	spdat <- spdat[,matchnames]

	# Week included:
	if("week" %in% names(spdat)){

		DT <- as.data.table(spdat)
		setkey(DT,site,year,week)
		DT2 <- DT[CJ( 
			unique(site),
			seq(min(year),max(year)),
			seq(min(week),max(week))
			)]

		CDT <- merge(DT2,SDT2)

		}

	# Week not included:
	if(! "week" %in% names(spdat)){	

		if(length(brood)==1){

			DT <- droplevels(subset(spdat,brood==paste(brood)))
			DTb <- as.data.table(DT[,names(DT)!="brood"])
			setkey(DTb,site,year)
			DT2 <- DTb[
					CJ(unique(site),
						seq(min(year),max(year))
						),
					allow.cartesian=T]

			}

		if(length(brood)==2){

			DT <- droplevels(subset(spdat,brood %in% c("1","2")))
				# for some species, brood=NULL; may cause probs if included
			DTb <- as.data.table(DT)
			setkey(DTb,brood,site,year)
			DT2 <- DTb[
				CJ(unique(brood),
					unique(site),
					seq(min(year),max(year))
					),	
				allow.cartesian=T]

			}
		
		CDT <- merge(DT2,SDT2)

		}

	CDT$count[is.na(CDT$count)] <- 0
	CDT$count[is.na(CDT$sumcount)] <- NA

	DT3 <- droplevels(as.data.frame(CDT))
	DT3 <- DT3[,names(DT3)!="sumcount"]

	if("week" %in% names(DT3)) DT7 <- DT3
	
	if(! "week" %in% names(DT3)){

		excluderun <- names(DT3)[grep("run",names(DT3))]
		excludenames <- c(c("year","week","gridno","brood","tlen","north","east"),excluderun)
			# exclude names and remove last row:
		
		dlist <- lapply(split(DT3,brood),function(d){

			# do once for each brood:

			preDT3 <- d[-nrow(d),! names(d) %in% excludenames]
				# add in extra row at top: 
			preDT4 <- rbind(rep(NA,ncol(preDT3)),preDT3)
				# this shifts entire dataframe down one row, and because data.table
				# has ordered, means each row in DT3 will be next to previous
			siteseq <- as.numeric(levels(preDT4$site))
			realsites <- as.numeric(as.character(preDT4$site))
			matchfirstyear <- match(siteseq,realsites)
			preDT5 <- as.matrix(preDT4[,names(preDT4)!="site"])
				# have to use matrix because want to add rows of missing values
			preDT5[matchfirstyear,] <- NA

			return(preDT5)
	
			})

		if(length(dlist)==1){
			preDT5b <- dlist[[1]] # only one brood, so no need to rbind
			}
		if(length(dlist)==2){
			preDT5b <- rbind(dlist[[1]],dlist[[2]]) # rbind separate brood dataframes
			}

		preDT6 <- as.data.frame(preDT5b) 
		names(preDT6) <- paste("pre",names(DT3)[!names(DT3) %in% c(excludenames,"site")],sep="")
		# names(preDT4[names(preDT4)!="site"])

		DT7 <- data.frame(DT2,preDT6)

		DT7$lprecount <- with(DT7, log(precount))
		DT7$dens <- with(DT7, count/tlen)
		DT7$lpredens <- with(DT7, log(precount/tlen))
		DT7$growth <- with(DT7, count/precount)
		
		}
		
	return(DT7)

	# head(cbind(DT7$year,DT7$site,DT7$count,DT7$precount,DT7$temp1,DT7$pretemp1),1000)

	}

########################################################
#subcounts 
#
#subsets data according to:
#	taucounts = threshold number of counts to be included as a yea
#	tauyears = threshold number of years with counts equal to or above that count
########################################################

subcounts <- function(data,taucounts,tauyears){

	abovecount <- function(x) sum(na.omit(x)>=taucounts)
	presyears <- with(data, tapply(count,list(site,year),abovecount))
	npresyear <- apply(presyears,1,function(x){
		sum(na.omit(x))
		})
	keepsites <- names(which(npresyear>tauyears))
	data2 <- subset(data,site %in% keepsites)
	data2 <- droplevels(data2)
	
	return(data2)

	}

#########################################################
#gridcat
#
#categorises each patch by a grid cell of a certain size 
#
############################################################

#indata <- site2
#gridsize <- 5000

gridcat <- function(indata, gridsize){
	coordcut <- function(x){
		cut(x, 
			seq(signif(min(x,na.rm=T)-gridsize, 2), 
			signif(max(x,na.rm=T)+gridsize, 2), gridsize)
			)
		}
	xcut <- coordcut(indata$east)
	ycut <- coordcut(indata$north)
	coordtab <- table(xcut, ycut)
	gridmat <- matrix(1:(length(levels(xcut))*length(levels(ycut))), ncol=length(colnames(coordtab)))
	xpos <- match(xcut, rownames(coordtab)) 
	ypos <- match(ycut, colnames(coordtab))
	gridno <- vector()
	for(i in 1:length(indata$east)){
		gridno[i] <- gridmat[xpos[i], ypos[i]]
		}
	#gridno[is.na(gridno)] <- round(
	#	runif(sum(is.na(gridno)),1,max(na.omit(gridno))),
	#	0)	# randomly assigns - ignored here
	gridno <- as.factor(gridno)
	gridno
	}

############
# GRAPHICS #
############

### PAIRS CORRELATION

panel.cor.scale <- function(x, y, digits=2, prefix="", cex.cor)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r = (cor(x, y,use="pairwise"))
txt <- format(c(r, 0.123456789), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
text(0.5, 0.5, txt, cex = cex * abs(r))
}


panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r = (cor(x, y,use="pairwise"))
txt <- format(c(r, 0.123456789), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
text(0.5, 0.5, txt, cex = cex )
}


panel.hist <- function(x, ...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )
h <- hist(x, plot = FALSE)
breaks <- h$breaks; nB <- length(breaks)
y <- h$counts; y <- y/max(y)
rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE)
{if (smooth ){
if (scale) {
pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=panel.smooth)
}
else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
} #else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
}
else #smooth is not true
{ if (scale) {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale)
} else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor) }
} #end of else (smooth)
} #end of function






