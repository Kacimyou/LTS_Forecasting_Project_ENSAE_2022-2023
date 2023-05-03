########################################
###             Packages             ###
########################################


require(forecast)
require(zoo)
require(tseries)
require(fUnitRoots) 

########################################
### Import the data and set the path ###
########################################
path <- "C:\\Users\\youns\\Documents\\GitHub\\LTS_Forecasting_Project_ENSAE_2022-2023"
setwd(path) #definit l'espace de travail (working directory ou "wd")
getwd() #affiche le wd
list.files() #liste les elements du wd

datafile <- "valeurs_mensuelles.csv" #definit le fichier de donnees
data <- read.csv(datafile,sep=";") #importe un fichier .csv dans un objet de classe data.frame



# Define time and interest variable #

dates_char <- as.character(data[[1]])
T <- length(dates_char)
dates_char <- dates_char[4:T]
dates_char[1]
tail(dates_char,1)
dates <- as.yearmon(seq(from=2008+10/12, to=2023+3/12, by=1/12))
y <-zoo(data[[2]][4:T], order.by=dates)
y_num = as.numeric(y)
y_diff = diff(y_num, differences = 1)#first difference
y_diff = zoo(y_diff)

y_second_diff = diff(y_diff, differences = 1)#second difference
y_second_diff = zoo(y_second_diff)


########################################
### Plot the data ###
########################################

#plot(y, ylim=c(60,160))  

plot(y, ylim =c(60,160))
plot(y_diff)
plot(y_second_diff)

acf(y_num)

summary(lm( formula = y ~ dates))





summary(lm( formula = y_diff ~ dates[-1]))
adf <- adfTest_valid(y_diff,24, type="nc")
adf


########################################
### ARMA modelisation ###
########################################
#Since our data is now stationary, we can try to model it with an ARMA(p,q) model.

### Identification of p and q ###
par(mfrow=c(1,1)) #puts the graphs into 1 column and 2 lines
acf(y_diff)
pacf(y_diff)

#The ACF and PACF diagram suggests to choose p=2 and q=2 for the ARMA model. 

### Nullity test of the ARMA coefficients ###


arima(y_wo_trend,c(2,0,2)) # p= 2 q= 2 
#The model is properly adjusted since the t statistics of each coefficient is superior to 1.96

arima202 <- arima(y_wo_trend,c(2,0,2)) 

### Ljung Box test to check if residuals are not correlated ####
Box.test(arima202$residuals, lag=5, type="Ljung-Box", fitdf=4) #
#Ljung Box reject the absence of autocorrelation of residuals at order 8 because p_value = 0.015
# The model is thus not valid.


Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
#We do LB test for two periodicity (24 tests since data are monthly)
round(Qtests(arima202$residuals,24,fitdf=4),3) # fitdf = p+q
# In this test, Ho: residuals are not correlated VS H1: residuals are correlated
# Thus we would like to accept H0 in order to have a valid model

adfTest_valid <-
  function(series,kmax,type){ #ADF tests until no more autocorrelated residuals
    k <- 0
    noautocorr <- 0
    while (noautocorr==0){
      cat(paste0("ADF with ",k, " lags: residuals OK? "))
      adf <- adfTest(series,lags=k,type=type)
      pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
      if (sum(pvals<0.05,na.rm=T) == 0) {
        noautocorr <- 1; cat("OK \n")}
      else cat("nope \n")
      k <- k + 1
    }
    return(adf)
  }




### Model selection  ###

# We know that our model is at most an ARMA(5,2) but given the test above this model is not
# well adjusted and not valid. We want to find simpler model that are valid and well adjusted 

#Function that tests the significance of coefficients 
#Allow us to test if the model is well adjusted 

signif <- function(estim){ 
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

signif(arima202) #

#Function that apply "signif" function and tests if the residuals are correlated and 
#Allow us to test if the model is valid and well adjusted. 

arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("Nullity test of the coefficients :\n")
  print(adjust)
  cat("\n Test of absence of residuals autocorrelation : \n")
  print(pvals)
}
# We apply this function on every simpler model that respect p < p_star = 5 and q < q_star = 2


estim <- arima(y_wo_trend,c(1,0,0)); arimafit(estim) # Nope:The model is not valid

estim <- arima(y_wo_trend,c(2,0,0)); arimafit(estim) # OK:The model is well adjusted and valid
ar2 <- estim

estim <- arima(y_wo_trend,c(0,0,1)); arimafit(estim) # Nope:The model is not valid

estim <- arima(y_wo_trend,c(0,0,2)); arimafit(estim) # OK:The model is well adjusted and valid
ma2 <- estim

estim <- arima(y_wo_trend,c(1,0,1)); arimafit(estim) # OK:The model is well adjusted and valid
ar1ma1 <- estim

estim <- arima(y_wo_trend,c(1,0,2)); arimafit(estim) # Nope:The model is not properly adjusted

estim <- arima(y_wo_trend,c(2,0,1)); arimafit(estim) # Nope:The model is not properly adjusted

estim <- arima(y_wo_trend,c(2,0,2)); arimafit(estim) # OK:The model is well adjusted and valid
ar2ma2 <- estim


### Compute the AIC BIC matrix of valid models #####
models <- c("ar2","ma2",'ar1ma1',"ar2ma2"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))

#

#### Q7 ####

##
models <-  c("ar2","ma2",'ar1ma1',"ar2ma2")
preds <- zoo(matrix(NA,ncol=3,nrow=4),order.by=tail(index(y.source),4))
colnames(preds) <- models
desaisonp <- preds #
xmp <- preds #

##
for (m in models){
  pred1 <- mean(desaison) + zoo(predict(get(m),4)$pred, order.by=tail(index(xm.source),4))
  pred2 <- as.numeric(tail(xm,12))[1:4] + pred1
  desaisonp[,m] <- pred1
  xmp[,m] <- pred2
}

obs <- tail(xm.source,4) #
cbind(obs,xmp) #
apply(xmp,2, function(x) sqrt(sum((x-obs)^2)/4)/sd(xm.source)) #

#

#### Q8 ####
datafile <- "Donnees2.csv" #definit le fichier de donnees

data <- read.csv(datafile,sep=";") #importe un fichier .csv dans un objet de classe data.frame
xm.source <- zoo(data[[1]]) #convertit le premier element de data en serie temporelle de type "zoo"
T <- length(xm.source)
xm <- xm.source[1:(T-4)] #supprime les 4 dernieres valeurs
dev.off() #reinitialise les parametre de graphique
plot(xm)
### 

trend <- 1:length(xm)
lt <- lm(xm ~ trend) #
summary(lt) #
r <- lt$residuals #
par(mfrow=c(1,2))
plot(r)
acf(r)
### 

pp.test(xm) 
### 

acf(r,24);pacf(r,24) 
### 
### 
pmax=4; qmax=21

### 



## fonction pour estimer un arima et en verifier l'ajustement et la validite
modelchoice <- function(p,q,data=r, k=24){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

## fonction pour estimer et verifier tous les arima(p,q) avec p<=pmax et q<=max
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}

armamodels <- armamodelchoice(pmax,qmax) #estime tous les arima (patienter...)


selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modeles bien ajustes et valides
selec
### On a ? modeles bien ajustes et valides

pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) #cree une liste des ordres p et q des modeles candidats
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")") #renomme les elements de la liste
models <- lapply(pqs, function(pq) arima(r,c(pq[["p"]],0,pq[["q"]]))) #cree une liste des modeles candidats estimes
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) #calcule les AIC et BIC des modeles candidats
### L'ARMA(?,?) minimise les criteres d'information.

rps <- lapply(models, function(m) as.zoo(predict(m,4)$pred)) #previsions de r
xmps <- lapply(rps, function(rp) rp+cbind(1,c((T-3):T))%*%lt$coefficients) #previsions de xm
rmse <- vapply(xmps, FUN.VALUE=numeric(1), function(xmp) sqrt(sum((as.zoo(xmp)-tail(xm.source,4))^2))) #calcule les rmse out-of-sample
rmse
### L'ARMA(?,?) fait aussi la meilleure prevision
