#################
### LOAD DATA ###
#################

cdc <- read.csv("Data/CDC/regional_averages_tm_02.csv",header=T,skip=1)

martemps <- cdc$Sachsen
year <- cdc$Jahr

pcex <- 0.65

col1 <- rgb(0,102,153,maxColorValue=255)
col2 <- rgb(226,240,214,maxColorValue=255)
col3 <- "white"

climplot <- function(points,trend){
	par(mar=c(5,5,2,2)+0.1,bg=col1,fg=col3,col.axis=col3)
	plot(martemps~year,
	type=ifelse(points==T,"l","n"),
	xlab="", # "Year",
	ylab="", # "Mean February Temperature (°C)",
	las=1,
	pch=16,
	bty="l",
	col=col2
	)
	if(points==T) points(martemps~year,pch=16,cex=pcex)
	if(trend==T) abline(lm(martemps~year),lty=2,col=col3)
	}

pdf(paste0("Plots/Feb_Meantemps_CDC_",format(Sys.Date(),"%d%b%Y"),".pdf"),
	width=6,height=4)

climplot(points=F,trend=F)
climplot(points=T,trend=F)
climplot(points=T,trend=T)

dev.off()


