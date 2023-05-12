arima(y_wo_trend,c(2,0,2)) # regresses the ARIMA(2,0,2) for the y series
#The model is properly adjusted since the t statistics of each coefficient is superior to 1.96

arima202 <- arima(y_wo_trend,c(2,0,2)) 

### Ljung Box test to check if residuals are not correlated ####
Box.test(arima202$residuals, lag=5, type="Ljung-Box", fitdf=4) #
#Ljung Box reject the absence of autocorrelation of residuals at order 8 because p_value = 0.015
# The model is thus not valid.



#We do LB test for two periodicity (24 tests since data are monthly)
round(Qtests(arima202$residuals,24,fitdf=4),3) # fitdf = p+q
# In this test, Ho: residuals are not correlated VS H1: residuals are correlated
# Thus we would like to accept H0 in order to have a valid model


ar2i1ma0

#Q6)
y.source <- y

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
plot(cbind(obs,y_pred)) #
apply(y_pred,2, function(x) sqrt(sum((x-obs)^2)/4)/sd(y.source)) #

#

