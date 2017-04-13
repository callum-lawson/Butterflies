setwd("C:/Users/Callum/Documents/Dropbox/PhD/Analyses/BMS")

###
###        # convert UKBMS site 1km coordinates to 5km references
library(maptools)
sites<-read.csv("ukbms_site.csv",header=T)
head(sites)
for (i in 1:length(sites$EAST)){
ifelse(sites$EAST[i]%%10000<5000,y<-2500,y<-7500)
sites$east.5k[i]<-(sites$EAST[i]%/%10000)*10000+y
}
for (i in 1:length(sites$NORTH)){
ifelse(sites$NORTH[i]%%10000<5000,y<-2500,y<-7500)
sites$north.5k[i]<-(sites$NORTH[i]%/%10000)*10000+y
}

nrow(sites)      # 1465 sites #DOESN'T MATCH TOM'S - CHECK


#  MEAN TEMPERATURE
file.names<-list.files("Met office climate data\\Meantemp_2007_2013",pattern=".txt")
library(maptools)
final.table<-NULL
for (f in file.names[]){
temp.table<-NULL
cip.tab<- data.frame(readAsciiGrid(paste("Met office climate data\\Meantemp_2007_2013\\",f,sep="")))
names(cip.tab)[1]<-"mean.temp"
cip.tab$year<-substr(f,10,13)
cip.tab$month<-substr(f,15,16)
# cip.tab[cip.tab$s1==587500&cip.tab$s2==347500,]
#str(cip.tab)
#with(cip.tab,plot(s1,s2))   # s1 = east s2 = north
# points(587500,347500,col=2)
#head(cip.tab)
#nrow(cip.tab)
temp.table<-merge(sites,cip.tab,by.x=c("east.5k","north.5k"),by.y=c("s1","s2"))    # 1435 rows
final.table<-rbind(final.table,temp.table)
}
nrow(final.table)
write.table(final.table,"all.sites.CIP.mean.temp.2007_2013.txt",sep="\t",col.names=NA)

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


#####
######
# RAINFALL
file.names<-list.files("C:\\TOLIVER\\BMS_DATABASES\\UKCIP_monthly_1971_2006\\Rainfall",pattern=".txt")

final.table<-NULL
for (f in file.names[]){
temp.table<-NULL
cip.tab<- data.frame(readAsciiGrid(paste("C:\\TOLIVER\\BMS_DATABASES\\UKCIP_monthly_1971_2006\\Rainfall\\",f,sep="")))
names(cip.tab)[1]<-"rainfall"
cip.tab$year<-substr(f,10,13)
cip.tab$month<-substr(f,15,16)
#str(cip.tab)
#with(cip.tab,plot(s1,s2))
#head(cip.tab)
#nrow(cip.tab)
temp.table<-merge(sites,cip.tab,by.x=c("east.5k","north.5k"),by.y=c("s1","s2"))    # 1435 rows
final.table<-rbind(final.table,temp.table)
}
nrow(final.table)
head(final.table)

 # write.table(final.table,"C:\\TOLIVER\\R\\OUTPUT\\all.sites.CIP.rainfall.1971_2006.txt",sep="\t",col.names=NA)

#####

