library(cluster)
library(TSclust)
library(tseries)
library(factoextra)
library(tidyverse)
##Hirarki Clustering
#Import data
data<- read.table("E:/BAHAN SKRIPSI/Time Series Clustering/Nee/Dasarian nomiss.txt", header = T)
par(mfrow=c(4,5))
View(data)
datats<- ts(data[,2:35])
ts.plot(data=data$METEO.BENGKULU, ylab="CH", main = "Meteo Bengkulu")
ts.plot(data=data$PADANG.HARAPAN, ylab="CH", main = "Padang Harapan")
ts.plot(data=data$STASIUN.KLIMATOLOGI.PULAU.BAAI.BENGKULU, ylab= "CH", main = "Statklim Pulau Baai")
ts.plot(data=data$AGRISINAL, ylab="CH", main = "Agrisinal")
ts.plot(data=data$BATU.BANDUNG, ylab="CH", main = "Batu Bandung")
ts.plot(data=data$BUNGA.MAS, ylab="CH", main = "Bunga Mas")
ts.plot(data=data$MANNA, ylab="CH", main = "Manna")
ts.plot(data=data$SEGINIM, ylab="CH", main = "Seginim")
ts.plot(data=data$SELALI, ylab="CH", main = "Selali")
ts.plot(data=data$SULAU, ylab= "CH", main ="Sulau")
ts.plot(data=data$TALANG.PAUH, ylab="CH", main = "Talang Pauh")
ts.plot(data=data$ARGAMAKMUR, ylab="CH", main = "Argamakmur")
ts.plot(data=data$KARANG.PULAU, ylab="CH", main = "Karang Pulau")
ts.plot(data=data$KURO.TIDUR, ylab="CH", main = "Kuro Tidur")
ts.plot(data=data$PT.Julang.Oca.Permana..JOP..Sebayur, ylab="CH", main = "PT. JOP Sebayur")
ts.plot(data=data$KANPEL.LINAU, ylab="CH", main = "Kanpel Linau")
ts.plot(data=data$KAUR, ylab="CH", main = "Kaur")
ts.plot(data=data$KABAWETAN, ylab="CH", main = "Kabawetan")
ts.plot(data=data$KELOBAK, ylab="CH", main = "Kelobak")
ts.plot(data=data$STASIUN.GEOFISIKA.KEPAHYANG, ylab="CH", main = "Geofisika Kepahyang")
ts.plot(data=data$TEBAT.KARAI, ylab="CH", main = "Tebat Karai")
ts.plot(data=data$TES, ylab="CH", main = "TES")
ts.plot(data=data$AIR.MANJUNTO, ylab="CH", main = "Air Manjunto")
ts.plot(data=data$PENARIK, ylab="CH", main = "Penarik")
ts.plot(data=data$AIR.BENING, ylab="CH", main = "Air Bening")
ts.plot(data=data$BUKIT.KABA, ylab="CH", main = "Bukit Kaba")
ts.plot(data=data$CURUP.DIPERTA, ylab="CH", main = "Curup Diperta")
ts.plot(data=data$PAL.8, ylab="CH", main = "Pal.8")
ts.plot(data=data$PT.AGRO.TEH.BUKIT.DAUN, ylab="CH", main = "PT. Agro Teh Bukit Daun")
ts.plot(data=data$KEMBANG.MUMPO, ylab="CH", main = "Kembang Mumpo")
ts.plot(data=data$MASMAMBANG, ylab="CH", main = "Masmambang")
ts.plot(data=data$RIMBO.KEDUI, ylab="CH", main = "Rimbo Kedui")
ts.plot(data=data$SUKARAJA, ylab="CH", main = "Sukaraja")
ts.plot(data=data$TALANG.DANTUK, ylab="CH", main = "Talang Dantuk")
View(datats)
#Hitung jarak/distance dengan euclidean dan DTW
disdtw<-diss(datats, METHOD = "DTW")
disdtw
matriksdtw<-as.matrix(disdtw)
View(matriksdtw)
diseuc<-diss(datats, METHOD = "EUCL")
diseuc
matrikseuc<-as.matrix(diseuc)
View(matrikseuc)
#pengklasteran hirarki dengan jarak DTW
cdtw1<-hclust(disdtw, method = "average")
cdtw2<-hclust(disdtw, method = "single")
cdtw3<-hclust(disdtw, method = "complete")
plot(cdtw1, main = "Dendogram Average Jarak DTW", cex=0.8)
plot(cdtw2, main = "Dendogram Single Jarak DTW", cex=0.8)
plot(cdtw3, main = "Dendogram Complete Jarak DTW",cex=0.8)
#pengklasteran hirarki dengan jarak euclidean
ceuc1<-hclust(diseuc, method = "average")
ceuc2<-hclust(diseuc, method = "single")
ceuc3<-hclust(diseuc, method = "complete")
plot(ceuc1, main = "Dendogram Average Jarak Euclidean", cex=0.8)
plot(ceuc2, main = "Dendogram Single Jarak Euclidean", cex=0.8)
plot(ceuc3, main = "Dendogram Complete Jarak Euclidean", cex=0.8)
#perbandingan nilai korelasi cophenetic antar metode hirarki
cave<-cophenetic(cdtw1)
cdtwave = cor(disdtw,cave)
cdtwave
csin<-cophenetic(cdtw2)
cdtwsin = cor(disdtw,csin)
cdtwsin
ccom<-cophenetic(cdtw3)
cdtwcom = cor(disdtw,ccom)
cdtwcom
eave<-cophenetic(ceuc1)
ceucave = cor(diseuc,eave)
ceucave
esin<-cophenetic(ceuc2)
ceucsin = cor(diseuc,esin)
ceucsin
ecom<-cophenetic(ceuc3)
ceuccom = cor(diseuc,ecom)
ceuccom
#K optimum
ck<-fviz_nbclust(datats, kmeans, method = "silhouette")
ck
ck$data
#KMEANS
set.seed(3434)
kmecl<-kmeans(disdtw, 2)
kmecl
fviz_cluster(kmecl, data = disdtw)+ggtitle("K=2")
finaltabel=data.frame(kmecl$cluster)
view(finaltabel)

#ARIMA
library(tseries)
library(forecast)
library(lmtest)
library(FitAR)
library(stats)
library(EnvStats)
library(TSA)
library(DescTools)


#Cluster 1
cluster1<- read.table("E:/BAHAN SKRIPSI/Time Series Clustering/Nee/C1K2.txt", header = T)
tsc1<-ts(cluster1)
str(tsc1)
View(tsc1)
datrain<-tsc1[1:180,]
datest<-tsc1[181:183,]
#Identifikasi Stasioner dalam varians
ts.plot(datrain)
cx<-BoxCox.ar(datrain)
nlambda1<-as.matrix(c(cx$lambda))
nloglik1<-as.matrix(c(cx$loglike))
nbind1<-cbind(nlambda1,nloglik1)
lambda1<-cx$mle
#lambda<-BoxCox.lambda(datrain, method = "loglik", lower = -2, upper = 2)
#plot.default(bcx$lambda, bcx$objective, type = "l")
#bcx<-boxcox(datrain, objective.name = "Log-Likelihood")
#plot(bcx$lambda, bcx$objective, type = "l", xlab = "Lambda", ylab = "Log-Likelihood")
#plot(bcx, ylim = c(-700,-300))
#BoxCox.ts(datrain)

#transformasi data jika tidak stasioner dalam varians
tdata<-datrain^lambda1
View(as.matrix(tdata)) #nilai hasil transformasi
tlambda1<-BoxCox.ar(tdata)
tlambda1<-tlambda1$mle
#Identifikasi stasioner dalam rataan
adf.test(tdata)#pvalue lebih kecil dari alpha maka stasioner
#melakukan differencing data tidak stasioner dalam rataan
ddata<-diff(tdata, differences = 1)
View(as.matrix(ddata))
adf.test(ddata)
#Menentukan orde ARIMA(p,d,q)
#MA melihat plot ACF terpotong pada lag ke-q
#AR melihat plot PACF terpotong pada lag ke-p
par(mfrow=c(1,1))
acf(ddata)
pacf(ddata)
#AR(1),AR(2), MA(1),MA(2),MA(3), MA(4), MA(5)

#estimasi parameter
#model ARIMA (0,1,1)
model1=Arima(ddata, order=c(0,0,1))
coeftest(model1)
summary(model1)
#model ARIMA (1,1,0)
model2=Arima(ddata, order=c(1,0,0))
coeftest(model2)
summary(model2)
#model ARIMA (2,1,0)
model3=Arima(ddata, order=c(2,0,0))
coeftest(model3)
summary(model3)
#model ARIMA (1,1,1)
model4=Arima(ddata, order=c(1,0,1))
coeftest(model4)
summary(model4)
#model ARIMA (2,1,1)
model5=Arima(ddata, order=c(2,0,1))
coeftest(model5)
summary(model5)
#parameter dikatakan signifikan jika pvalue < alpha

##Uji Diagnosa Model
#plot residual
ts.plot(model1$residuals)
#uji Normalitas visual residuals
library(nortest)
library(normtest)
qqnorm(model1$residuals)
qqline(model1$residuals)
jb.norm.test(model1$residuals)
#autokorelasi residuals
acf(model1$residuals, lag=10) #jika tidak terdapat plot ACF terpotong pada lag, maka tidak ada autokorelasi di residuals
#Ljung-Box Test
Box.test(model1$residual,lag=3, type="Ljung")
#model dikatakan baik jika nilai p valu Ljung Box lebih besar dari alpha

##Forecasting
nilaiprediksi<-forecast(model1, h=3, level=c(99.5))
nilaiprediksi
#Melakukan re-transformasi data hasil prediksi
pred<-as.data.frame(nilaiprediksi)
nilaiprediksi<-pred$`Point Forecast`+tsc1[180:182]
hasilprediksi<-data.frame(nilaiprediksi)
hasilprediksi
#Menggabungkan data hasil prediksi kedalam data training
datrain<-data.frame(Data_Training=datrain)
nilaipredik<-data.frame(Data_Training=nilaiprediksi)
hasilprediksi<-rbind(datrain,nilaipredik)
hasilprediksi<-ts(hasilprediksi)
hasilprediksi
#plot data hasil prediksi
plot(hasilprediksi, type = "l", col = "red")
lines(datrain, type = "l", col = "blue")
plot(datest, type = "l", col = "blue", xlab = "periode", ylab = "data", ylim = c(80,160))
lines(nilaiprediksi, type = "l", col = "red")
legend("topright", legend = c("Nilai Prediksi", "Nilai Aktual"), col = c("red", "blue"), lty = 1:1, cex = 0.6, title = "Lines Type", bg = 'transparent')
#Menghitung nilai MAPE
library(MLmetrics)
nilaimape<-MAPE(nilaiprediksi,datest)*100
nilaimape

#Cluster 2
cluster2<- read.table("E:/BAHAN SKRIPSI/Time Series Clustering/Nee/C2K2.txt", header = T)
tsc2<-ts(cluster2)
str(tsc2)
View(tsc2)
datrain2<-tsc2[1:180,]
datest2<-tsc2[181:183,]
#Identifikasi Stasioner dalam varians
ts.plot(datrain2)
cx2<-BoxCox.ar(datrain2)
nlambda2<-as.matrix(c(cx2$lambda))
nloglik2<-as.matrix(c(cx2$loglike))
nbind2<-cbind(nlambda2,nloglik2)
lambda2<-cx2$mle
#BoxCox.ar(datrain2)
#boxcox(datrain2)
#transformasi data jika tidak stasioner dalam varians
tdata2<-datrain2^lambda2
View(as.matrix(tdata2)) #nilai hasil transformasi
tlambda2<-BoxCox.ar(tdata2)
tlambda2<-tlambda2$mle
#boxcox(tdata2)
#Identifikasi stasioner dalam rataan
adf.test(tdata2)#pvalue lebih kecil dari alpha maka stasioner
#melakukan differencing data tidak stasioner dalam rataan
#ddata2<-diff(tdata2, differences = 1)
#View(as.matrix(ddata2))
#adf.test(ddata2)
#Menentukan orde ARIMA(p,d,q)
#MA melihat plot ACF terpotong pada lag ke-q
#AR melihat plot PACF terpotong pada lag ke-p
par(mfrow=c(2,1))
acf(tdata2)
pacf(tdata2)
#MA(1),MA(2), AR(1)

#estimasi parameter
#model ARIMA (1,0,0)
model1=Arima(tdata2, order=c(1,0,0))
coeftest(model1)
summary(model1)
#model ARIMA (2,0,0)
model2=Arima(tdata2, order=c(2,0,0))
coeftest(model2)
summary(model2)
#model ARIMA (3,0,0)
model3=Arima(tdata2, order=c(3,0,0))
coeftest(model3)
summary(model3)
#model ARIMA (0,0,1)
model4=Arima(tdata2, order=c(0,0,1))
coeftest(model4)
summary(model4)
#model ARIMA (0,0,2)
model5=Arima(tdata2, order=c(0,0,2))
coeftest(model5)
summary(model5)
#model ARIMA (0,0,3)
model6=Arima(tdata2, order=c(0,0,3))
coeftest(model6)
summary(model6)
#model ARIMA (1,0,1)
model7=Arima(tdata2, order=c(1,0,1))
coeftest(model7)
summary(model7)
#model ARIMA (1,0,2)
model8=Arima(tdata2, order=c(1,0,2))
coeftest(model8)
summary(model8)
#model ARIMA (1,0,3)
model9=Arima(tdata2, order=c(1,0,3))
coeftest(model9)
summary(model9)
#model ARIMA (2,0,1)
model10=Arima(tdata2, order=c(2,0,1))
coeftest(model10)
summary(model10)
#model ARIMA (2,0,2)
model11=Arima(tdata2, order=c(2,0,2))
coeftest(model11)
summary(model11)
#model ARIMA (2,0,3)
model12=Arima(tdata2, order=c(2,0,3))
coeftest(model12)
summary(model12)
#model ARIMA (3,0,1)
model13=Arima(tdata2, order=c(3,0,1))
coeftest(model13)
summary(model13)
#model ARIMA (3,0,2)
model14=Arima(tdata2, order=c(3,0,2))
coeftest(model14)
summary(model14)
#model ARIMA (3,0,3)
model15=Arima(tdata2, order=c(3,0,3))
coeftest(model15)
summary(model15)
#parameter dikatakan signifikan jika pvalue < alpha
par(mfrow=c(1,1))
##Uji Diagnosa Model
#plot residual
ts.plot(model13$residuals)
#uji Normalitas visual residuals
library(nortest)
library(normtest)
qqnorm(model13$residuals)
qqline(model13$residuals)
jb.norm.test(model13$residuals)
#autokorelasi residuals
acf(model13$residuals, lag=10)
#jika tidak terdapat plot ACF terpotong pada lag, maka tidak ada autokorelasi di residuals
#Ljung-Box Test
Box.test(model13$residual,lag=3, type="Ljung")
#model dikatakan baik jika nilai p valu Ljung Box lebih besar dari alpha

##Forecasting
nilaiprediksi2<-forecast(model13, h=3, level=c(99.5))
nilaiprediksi2
#Melakukan re-transformasi data hasil prediksi
pred2<-as.data.frame(nilaiprediksi2)
nilaiprediksi2<-pred2$`Point Forecast`^(1/0.5)
hasilprediksi2<-data.frame(nilaiprediksi2)
hasilprediksi2
#Menggabungkan data hasil prediksi kedalam data training
datrain2<-data.frame(Data_Training=datrain2)
nilaipredik2<-data.frame(Data_Training=nilaiprediksi2)
hasilprediksi2<-rbind(datrain2,nilaipredik2)
hasilprediksi2<-ts(hasilprediksi2)
hasilprediksi2
#plot data hasil prediksi
plot(hasilprediksi2, type = "l", col = "red")
lines(tsc2, type = "l", col = "blue")
plot(datest2, type = "l", col = "blue", xlab = "periode", ylab = "data")
lines(hasilprediksi2[181 :183,], type = "l", col = "red")
legend("topleft", legend = c("Nilai Prediksi", "Nilai Aktual"), col = c("red", "blue"), lty = 1:1, cex = 0.6, title = "Lines Type", bg = 'transparent')
#Menghitung nilai MAPE
MAPE(nilaiprediksi2,datest2)
nilaimape2<-MAPE(nilaiprediksi2,datest2)*100
nilaimape2
