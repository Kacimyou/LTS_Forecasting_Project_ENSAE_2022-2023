require(zoo)
require(tseries)
library("lubridate")
library("magrittr")

#nettoyage l'environnement de travail
rm(list=objects()) 
graphics.off()

#importation des données
list
path <- "C:\\Users\\mira_\\Documents\\ENSAE\\2A\\S2\\séries temp\\projet\\serie_010537368_09042023"
setwd(path)
getwd()
list.files()

datafile <- "valeurs_mensuelles-eau.csv"
data <- read.csv(datafile,sep=";") 

#liste des valeurs 
ch.source  <- zoo(data[[2]]) 
T = length(ch.source)
ch <- ch.source[4:T]

#axe des temps 
tps.source <- zoo(data[[1]])
tps <- tps.source[4:T]#on supprime les 3 premieres valeurs
tps <- ym(tps) #libradate, permet de convertir les strings en dates


plot.zoo(tps, ch, xaxt="n")
axis.Date(1, at = tps)
?axis.Date

plot.zoo(tps, ch)

plot(ch, ylim=c(60,140))             

df_tap = as.numeric(ch)
sup_trend_num = diff(df_tap, differences = 1)
sup_trend = zoo(sup_trend_num)
plot(sup_trend) 

pp.test(sup_trend)
