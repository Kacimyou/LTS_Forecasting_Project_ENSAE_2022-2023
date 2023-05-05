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
# path <- "C:\\Users\\youns\\Documents\\GitHub\\LTS_Forecasting_Project_ENSAE_2022-2023"
path <- "C:\\Users\\mira_\\Documents\\GitHub\\LTS_Forecasting_Project_ENSAE_2022-2023"
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

#Q2)
 
########################################
### Plot the data ###
########################################
 
plot(y, ylim =c(60,160))

#The series in level seems to have a increasing linear trend. 

plot(y_diff)

#The first difference series seems relatively stable around a null constant and could be stationary. 

#We guess that the series is probably I(1)


########################################
### Stationary ###
########################################
# Before modeling our series with ARMA model we need to check that it is stationary
# If it is not, we need to correct it by differentiating it or deseasonalizing it.


#Before performing the unit root tests to check stationarity, we need to check if there is an intercept and / or a non null linear trend.
#The graph representation of spread showed that the trend is probably linear and increasing.

# Let’s regress y on its dates to check :
summary(lm( formula = y ~ dates))


#The coefficient associated with the linear trend (dates) is indeed positive, thus we
#need to study the case of unit root tests with intercept and possibly non zero trends

adf <- adfTest(y_num, lag=0, type="ct") # ct here take into account the fact that y has an intercept and non zero trend.
#Before interpreting the test, let’s check that the model’s residuals are not autocorrelated, otherwise the test
#would not be valid.

Qtests(adf@test$lm$residuals,24,length(adf@test$lm$coefficients))

#We reject the absence of residual autocorrelation for every lag, thus invalidating the ADF
#test without lags. Let’s add lags of ∆Xt until the residuals are no longer autocorrelated.

adf <- adfTest_valid(y_num,24, type="ct") # ct here take into account the fact that y has an intercept and non zero trend.
#We have had to consider 6 lags on the ADF test to erase residual autocorrelation.

adf
#The unit root is not rejected at the 95% - level for the series in levels, the series is thus at least I(1).


#Let’s now test the unit root for the first differenciated series. The previous graph representation seems to
#show the absence of a constant and non zero trend. Let’s check with a regression :

summary(lm( formula = y_diff ~ dates[-1]))

#There isn’t any constant or significant trend. Let’s perform the ADF test in the no-constant and no-trend case,
#and control for the absence of residual autocorrelation.

adf <- adfTest_valid(y_diff,24, type="nc") # nc here take into account the fact that y has no intercept and zero trend.

#It was necessary to include 1 lags in the ADF test
adf

#The test rejects the unit root hypothesis (p-value<0.05), we will thus say that the differenciated series is
#”stationary”. y is therefore I(1) as guessed by the plot.


#Q3)

plot(y_diff)

### Useful function ###
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

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

#Q4)


########################################
### ARMA modelization ###
########################################
#Since our data is now stationary, we can try to model it with an ARMA(p,q) model.


### Identification of p and q ###
par(mfrow=c(1,2)) #puts the graphs into 1 column and 2 lines
acf(y_diff)
pacf(y_diff)

#Since the series is stationary, it is integrated of order d = 0.
#The complete autocorrelation functions are statistically significant (i.e. bigger than the bounds ±1, 96/√n of
#the confidence interval of a null test of the autocorrelation at the 95% level) until q∗ = 2 and the partial
#autocorrelation until p∗ = 2. If y follows an ARIMA(p,d,q), it follows at most an ARIMA(p∗ =2, d∗ = 0,q∗ = 2), which we can estimate.

### Model selection  ###

# We know that our model is at most an ARIMA(2,0,2)
#The potential models are all the ARIMA(p,0,q) for spread where p ≤ 2 and q ≤ 2.
#We are looking for a model that is :
# — well adjusted : the estimated coefficients (notably the coefficients of the higher AR and MA orders) are
# statistically significant.
# — valid : the residuals are not correlated.


#Function that tests the significance of coefficients 
#Allow us to test if the model is well adjusted 

signif <- function(estim){ 
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}
#Function that apply "signif" function and tests if the residuals are autocorrelated
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


# First we want to apply this function on every simpler model that respect p ≤ 2 and q ≤ 2
# to test if there exist some simpler model that are both valid and well adjusted.

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

# To choose between all the models that are well adjusted and valid we compute AIC and BIC of each model
# We will probably select both models with minimum AIC and BIC.

### Compute the AIC BIC matrix of valid models #####
models <- c("ar2","ma2",'ar1ma1',"ar2ma2"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))

#We can see that the model AR(2) has both minimum AIC and BIC, that's why we select it
# Thus our final model for the first difference series is an ARIMA(2,0,0)



#Q6)


models <-  c("ar2","ma2",'ar1ma1',"ar2ma2")
preds <- zoo(matrix(NA,ncol=4,nrow=4),order.by=tail(index(y.source),4))
colnames(preds) <- models
y_diff_pred <- preds #
y_pred <- preds #



for (m in models){
  pred1 <- mean(y_diff) + zoo(predict(get(m),4)$pred, order.by=tail(index(y.source),4))
  pred2 <- as.numeric(tail(y,12))[1:4] + pred1
  y_diff_pred[,m] <- pred1
  y_pred[,m] <- pred2
}

obs <- tail(y.source,4) #
cbind(obs,y_pred) #
apply(y_pred,2, function(x) sqrt(sum((x-obs)^2)/4)/sd(y.source)) #

#
