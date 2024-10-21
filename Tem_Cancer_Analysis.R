###################################
####load required packages
library(data.table);library(dplyr);library(dlnm);library(splines);library(tidyr)
library(tsModel);library(forestplot);library(reshape2);library(mvmeta);library(metafor)

###################################
####two stage analysis for total cancer deaths as an exampple
#load the data
t <-fread("D:/xxx/cancer_sample.csv") #exposure, cancer death, and city characteristics

#preparation stage
regions <- unique(t$City) 
dlist <- lapply(regions, function(x)t[t$City == x,])
names(dlist) <- regions

ranges <- t(sapply(dlist, function(x)range(x$Tem_avg, na.rm = T)))
bound <- colMeans(ranges)

#city-specific analysis
yall1 <- yall2 <- yall3 <- matrix(NA, length(dlist), 5, dimnames = list(regions, paste("b", seq(5), sep = "")))
sall1 <- sall2 <- sall3 <- vector("list", length(dlist))

for(i in seq(dlist)){
  print(paste(i, regions[i], sep="-")) 
  sub<-as.data.frame(dlist[[i]])
  sub$time <- c(1:length(sub$Tem_avg))
  cb1 <- crossbasis(sub$Tem_avg,lag=lag,argvar=list(fun="bs", degree=2, 
                                                    knots=quantile(sub$Tem_avg,perc/100,na.rm=T)), arglag = list(knots=logknots(lag,3)))
  mfirst <- glm(total ~ cb1 + as.factor(Dow) + ns(time, df=6*9) + ns(Rhu_avg,4), 
                family=quasipoisson(), data=sub, na.action="na.exclude")
  
  cen <- mean(sub$Tem_avg, na.rm=T)
  
  crall <- crossreduce(cb1, mfirst, cen=cen)
  yall1[i,] <- coef(crall)
  sall1[[i]] <- vcov(crall)
  
  crcold <- crossreduce(cb1, mfirst, type="var", value=-0.9, cen=cen)  #extreme cold temperature
  yall2[i,] <- coef(crcold)
  sall2[[i]] <- vcov(crcold)
  
  crheat <- crossreduce(cb1, mfirst, type="var", value=29.6, cen=cen)  #extreme heat temperature
  yall3[i,] <- coef(crheat)
  sall3[[i]] <- vcov(crheat)
}

#multivariate meta-analysis
method <- "reml"
control = list(showiter = T)
mvall <- mvmeta(yall1, sall1, method=method)
mvcold <- mvmeta(yall2, sall2, method=method)
mvheat <- mvmeta(yall3, sall3, method=method)

#multivariate meta-regression models
mvreall <- mvmeta(yall1 ~ lat + ... + ..., 
                S = sall1, method = method)
mvrecold <- mvmeta(yall2 ~ lat + ... + ..., 
                S = sall2, method = method)
mvreheat <- mvmeta(yall3 ~ lat + ... + ..., 
                S = sall3, method = method)

#prediction for simple meta-analysis with no predictors
xvar <- seq(bound[1], bound[2], by = 0.1) 
bvar <- do.call("onebasis", c(list(x = xvar), attr(cb1, "argvar"))) 
xlag <- 0:lag
blag <- do.call("onebasis", c(list(x = xlag), attr(cb1, "arglag")))  

cpall.cb <- crosspred(bvar, coef = coef(mvall), vcov = vcov(mvall),
                      model.link = "log", by=0.1, from=bound[1], to=bound[2])
MMT <- round(as.numeric(names(which.min(cpall.cb$allRRfit))),1)
cpall.cb <- crosspred(bvar, coef = coef(mvall), vcov = vcov(mvall),
                      model.link = "log", by=0.1, from=bound[1], to=bound[2], cen=MMT)
cpcold.cb <- crosspred(blag, coef = coef(mvcold), vcov = vcov(mvcold),
                         model.link = "log", at = 0:140/10, cen=MMT)
cpheat.cb <- crosspred(blag, coef = coef(mvhot), vcov = vcov(mvhot),
                         model.link = "log", at = 0:140/10, cen=MMT)


#prediction for multivariate meta-regression models: as an example for predictor of latitude
predrellat <- predict(mvreall,data.frame(lat=percmod[,"lat"]),vcov=T)
cprellat75 <- crosspred(bvar,coef=predrellat[[1]]$fit,
                        vcov=predrellat[[1]]$vcov,model.link="log",by=0.5)
cprellat25 <- crosspred(bvar,coef=predrellat[[2]]$fit,
                        vcov=predrellat[[2]]$vcov,model.link="log",by=0.5)

#wald test
fwald <- function(mfirst, var){
  ind <- grep(var, names(coef(mfirst)))
  coef <- coef(mfirst)[ind]
  vcov <- vcov(mfirst)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}

#heterogeneity test
qall <- qtest(mvall)  
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)  


####plot
xlab <- expression(paste("Temperature(",degree,"C)"))
ylab <- 'Relative Risk'
plot(cpall.cb, type="n", xlab = xlab, ylab = "Relative Risk", ylim=c(0.9,1.4))
title(main = "Total")
lines(cpall.cb, col=2, lwd=2)
abline(h=1, lwd=1)  