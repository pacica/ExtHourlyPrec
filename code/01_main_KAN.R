rm(list=ls())

setwd("~/Downloads/ExtHourlyPrec-main")

source("functions/tv.th.mh.R")

# load precipitation data for Kansas City
load("input/KAN.RData")

dt <- KAN

# adding covariate NAO
nao<-read.csv("input/nao.txt",sep="")
nao$YEAR<-as.numeric(rownames(nao))
NAO <- reshape(nao, idvar = "YEAR", varying=list(1:12), direction = "long")
NAO$time<-as.numeric(NAO$time)
names(NAO)<-c("year","month","NAO")

months <- unique(dt$month)
dt_merge <- merge(dt,NAO,by=c("year","month"),all.x=T)

dt_merge_1 <- dt_merge[with(dt_merge, order(as.Date(DATE))), ]

# tv threshold estimation
mode_th <- "hour"
quant_level <- c(0.95,0.96,0.97)
nq <- length(quant_level)
ylag <- 3
yburn <- 2005:2007

tvth <- vector(mode="list")

for (q in 1:nq){
  tvth <- append(tvth, list(tv.th.mh(dt_merge_1,quant_level[q],ylag,yburn,mode_th,T)))
}  

# save output
save.image("output/01_KAN_GPD_static.RData")

