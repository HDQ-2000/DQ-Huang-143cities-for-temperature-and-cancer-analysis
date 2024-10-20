###################################
####Load required packages#########
#################################
library(data.table);library(dplyr);library(dlnm);library(splines);library(tidyr)
library(tsModel);library(forestplot);library(reshape2);library(mvmeta);library(metafor);

####preparation stage
t <-fread("D:/xxx/cancer_sample.csv") #exposure, cancer death, and city characteristics
regions <- unique(t$City) 
dlist <- lapply(regions, function(x)t[t$City == x,])
names(dlist) <- regions

ranges <- t(sapply(dlist, function(x)range(x$Tem_avg, na.rm = T)))
bound <- colMeans(ranges)

#####two stage analysis for total cancer disease as an example#########
#####city-specific analysis
yall <- matrix(NA, length(dlist), 5, dimnames = list(regions, paste("b", seq(5), sep = "")))
sall <-  vector("list", length(dlist))

for(i in seq(dlist)){
  print(paste(i, regions[i], sep="-")) 
  sub<-as.data.frame(dlist[[i]])
  sub$time <- c(1:length(sub$Tem_avg))
  cb1 <- crossbasis(sub$Tem_avg,lag=lag,argvar=list(fun="bs", degree=2, 
                                                    knots=quantile(sub$Tem_avg,perc/100,na.rm=T)), arglag = list(knots=logknots(lag,3)))
  mfirst <- glm(total ~ cb1 + as.factor(Dow) + ns(time, df=6*9) + ns(Rhu_avg,4), 
                family=quasipoisson(), data=sub, na.action="na.exclude")
  cen <- mean(sub$Tem_avg, na.rm=T)
  pred <- crosspred(cb1, mfirst, cen=cen)
  crall <- crossreduce(cb1, mfirst, cen=cen)
  yall[i,] <- coef(crall)
  sall[[i]] <- vcov(crall)
}

#####meta-analysis
method <- "reml"
control = list(showiter = T)
mvall <- mvmeta(yall ~ lat + lon + region + climate + avgtmean + rangetmean + gdp, 
                S = sall, method = method)

blup <- blup(mvall, vcov=T)

sall2 <- vector("list", length(dlist))
for (i in seq(length(blup))) {
  sall2[[i]] <- blup[[i]]$vcov}
mvall2 <- mvmeta(blup ~ 1, S = sall2, method = method)

xvar <- seq(bound[1], bound[2], by = 0.1)  
bvar <- do.call("onebasis", c(list(x = xvar), attr(cb1, "argvar")))  
xlag <- 0:lag
blag <- do.call("onebasis", c(list(x = xlag), attr(cb1, "arglag")))  

cpall.cb <- crosspred(bvar, coef = coef(mvall2), vcov = vcov(mvall2),
                      model.link = "log", by=0.1, from=bound[1], to=bound[2])
MMT <- round(as.numeric(names(which.min(cpall.cb$allRRfit))),1)
cpall.cb <- crosspred(bvar, coef = coef(mvall2), vcov = vcov(mvall2),
                      model.link = "log", by=0.1, from=bound[1], to=bound[2], cen=MMT)

####plot
xlab <- expression(paste("Temperature(",degree,"C)"))
ylab <- 'Relative Risk'
plot(cpall.cb, type="n", xlab = xlab, ylab = "Relative Risk", ylim=c(0.9,1.4))
title(main = "Total")
lines(cpall.cb, col=2, lwd=2)
abline(h=1, lwd=1)  