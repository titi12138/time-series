##导入R包

library(tseries)
library(rugarch)
library(aTSA)
library(parallel)
library(broom)
library(fracdiff)
library(ltsa)
library(arfima)
library(forecast)
library(zoo)
library(xts)
library(psych)
library(fBasics)
library(timeDate)
library(timeSeries)
library(ggplot2)
library(dplyr)
library(pracma)
library(wesanderson)
library(MCS)


##有关于自定义函数

#定义高频数据
RVjisuan<-function(x,x1,x2,x3,fre) {
  #x为每fre分钟的股票收盘价
  #x1为每fre分钟的股票开盘价
  #x2为每日收盘价
  
  RVt_=c()#未转化的已实现波动率
  RVt =c()#转化后的已实现波动率
  InRv=c()#对数已实现波动率
  Rt  =c()#已实现波动率
  Rsd =c()#每天对数收益率的标准差
  
  Rt<-100*(diff(log(x2)))#得出日收益率，标识为相邻两个交易日收盘价的自然对数之差
  #Rt的长度是交易日长度减1
  
  Rt_2<-c(Rt^2)#得出日收益率的平方
  Rt_2<-round(Rt_2,digits = 9)
  
  len=length(x)#有多少个交易日
  day = len/fre#得出有多少个交易日
  
  for(i in 1:day) {
    y<-x[(((i-1)*fre)+1):(((i-1)*fre)+fre)]#取出第i天的所有高频数据
    y1<-x1[((i-1)*fre+1)]#取出第i天的第一个时间段的开盘价
    Itd<-c(y1,y)#将第一个开盘价和fre个收盘价组合在一起,长度为fre+1
    Rtd<-(diff(log(Itd))*100) #第i天的高频收益率
    Rsd[i] <- sd(Rtd)*10
    RVt_[i]<-sum(Rtd^2)#得到未转化的已实现波动率
  }
  
  Rt2_sum<-sum(Rt^2)
  RVt_sum<-sum(RVt_)
  gama<-(Rt2_sum/day)/(RVt_sum/day)#尺度函数对已实现波动率进行转化
  RVt <- gama*RVt_#得到已实现波动率
  InRv<- log(RVt)
  RVSD<- RVt/Rsd
  
  par(mfrow=c(1,2))
  plot(RVt,pch=15)#输出RVt的图像
  plot(InRv,pch=15)#输出InRV的图像
  
  RV1<-list(Rt,Rt_2,RVt,InRv,RVSD,gama)
  return(RV1)
}

```

#高频数据描述函数
mydescribe<-function(x){
  result<-data.frame()
  for (i in 1:5){
    a<-unlist(x[[i]])
    b<-describe(a)[-c(1,2,5,6,7,8,9,10,13)]
    JB <- normalTest(a,method='jb')@test$statistic#JB统计量
    JB_p<-normalTest(a,method='jb')@test$p.value#JB检验的P量
    Q_10 <-   Box.test(a,lag = 10,type = 'Ljung-Box')$statistic#LB统计量
    Q_10_p <- Box.test(a,lag = 10,type = 'Ljung-Box')$p.value#LB P量
    Q_20 <-   Box.test(a,lag = 20,type = 'Ljung-Box')$statistic#LB统计量
    Q_20_p <- Box.test(a,lag = 20,type = 'Ljung-Box')$p.value#LB P量
    Q_30 <-   Box.test(a,lag = 30,type = 'Ljung-Box')$statistic#LB统计量
    Q_30_p <- Box.test(a,lag = 30,type = 'Ljung-Box')$p.value#LB P量
    hurst <- hurstexp(a,display = FALSE)[[1]]
    #hurst指数，当0.5<h<1时，表示是长期记忆的时间序列
    b <- cbind(b,JB,JB_p,Q_10,Q_10_p,Q_20,Q_20_p,Q_30,Q_30_p,hurst)
    result<-rbind(result,b)
  }
  result<-t(result)
  return(result)
}


##数据导入

#数据导入
shuj1 <- read.csv("C:\\Users\\titi\\Desktop\\论文\\终数据\\shuj1.csv")
shuj5 <- read.csv("C:\\Users\\titi\\Desktop\\论文\\终数据\\shuj5.csv")
shuj15 <- read.csv("C:\\Users\\titi\\Desktop\\论文\\终数据\\shuj15.csv")
shuj30 <- read.csv("C:\\Users\\titi\\Desktop\\论文\\终数据\\shuj30.csv")
shuj60 <- read.csv("C:\\Users\\titi\\Desktop\\论文\\终数据\\shuj60.csv")
shujday<- read.csv("C:\\Users\\titi\\Desktop\\论文\\终数据\\shujday.csv")

x<-shuj5$close
plot(x,type = 'l')
plot(diff(x),type = 'l')


##序列输出

Rt <- RVjisuan(shuj5$close,shuj5$open,shujday$close,fre=48);Rt
Rt1 <- Rt[[1]]
tRt1 <- ts(Rt1)
Rt2 <- Rt[[2]]
tRt2 <- ts(Rt2)
RV <- Rt[[3]]
InRV <- Rt[[4]]
RVSD <- Rt[[5]]
Rt[[6]]
for (i in 1:5) print(length(Rt[[i]]))

##序列描述性统计量

des <- mydescribe(Rt)
colnames(des) <- c('Rt','Rt_2','RVt','InRv','RVSD')
des
#write.csv(des,file = 'C:\\Users\\titi\\Desktop\\论文\\数据\\describe1.csv')
```

##对RV-ARFIMA模型建立模型,根据AIC和BIC信息准则判断拟合模型精度

plot(RV,type='l')
tRV <- ts(RV)
par(mfrow=c(2,2))
mod.fit <- arfima::arfima(z=tRV,order = c(0,0,0))

summary(mod.fit)
acf(x= residuals(mod.fit)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit)[[1]],lag.max = 100)

mod.fit_1 <- arfima::arfima(z=tRV,order = c(1,0,0))
summary(mod.fit_1)
acf(x= residuals(mod.fit_1)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_1)[[1]],lag.max = 100)

mod.fit_2 <- arfima::arfima(z=tRV,order = c(0,0,1))
summary(mod.fit_2)
acf(x= residuals(mod.fit_2)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_2)[[1]],lag.max = 100)


mod.fit_3 <- arfima::arfima(z=tRV,order = c(1,0,1))
summary(mod.fit_3)
acf(x= residuals(mod.fit_3)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_3)[[1]],lag.max = 100)

AIC_RV <- t(c(AIC(mod.fit),AIC(mod.fit_1),AIC(mod.fit_2),AIC(mod.fit_3)));AIC_RV
BIC_RV <- t(c(BIC(mod.fit),BIC(mod.fit_1),BIC(mod.fit_2),BIC(mod.fit_3)));BIC_RV
IC <- rbind(AIC_RV,BIC_RV);IC


#综合以上结果 可以看出RV-AFRIMA(1,d,1)是其中拟合的最好的模型

##对InRV-ARFIMA模型建立模型,根据AIC和BIC信息准则判断拟合模型精度

plot(InRV,type='l')
tInRV <- ts(InRV)
mod.fit <- arfima::arfima(z=tInRV,order = c(0,0,0))
summary(mod.fit)
acf(x= residuals(mod.fit)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit)[[1]],lag.max = 100)

mod.fit_1 <- arfima::arfima(z=tInRV,order = c(1,0,0))
summary(mod.fit_1)
acf(x= residuals(mod.fit_1)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_1)[[1]],lag.max = 100)


mod.fit_2 <- arfima::arfima(z=tInRV,order = c(0,0,1))
summary(mod.fit_2)
acf(x= residuals(mod.fit_2)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_2)[[1]],lag.max = 100)


mod.fit_3 <- arfima::arfima(z=tInRV,order = c(1,0,1))
summary(mod.fit_3)
acf(x= residuals(mod.fit_3)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_3)[[1]],lag.max = 100)


AIC_RV <- t(c(AIC(mod.fit),AIC(mod.fit_1),AIC(mod.fit_2),AIC(mod.fit_3)))
BIC_RV <- t(c(BIC(mod.fit),BIC(mod.fit_1),BIC(mod.fit_2),BIC(mod.fit_3)))
IC <- rbind(AIC_RV,BIC_RV);IC

#综合以上结果，可以看出RV-AFRIMA(1,d,1)是其中拟合的最好的模型

##对RV/SD-ARFIMA模型建立模型,根据AIC和BIC信息准则判断拟合模型精度

plot(RVSD,type='l')
tRVSD <- ts(RVSD)
mod.fit <- arfima::arfima(z=tRVSD,order = c(0,0,0),quiet = T)
summary(mod.fit)
acf(x= residuals(mod.fit)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit)[[1]],lag.max = 100)

mod.fit_1 <- arfima::arfima(z=tRVSD,order = c(1,0,0))
summary(mod.fit_1)
acf(x= residuals(mod.fit_1)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_1)[[1]],lag.max = 100)


mod.fit_2 <- arfima::arfima(z=tRVSD,order = c(0,0,1))
summary(mod.fit_2)
acf(x= residuals(mod.fit_2)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_2)[[1]],lag.max = 100)


mod.fit_3 <- arfima::arfima(z=tRVSD,order = c(1,0,1))
summary(mod.fit_3)
acf(x= residuals(mod.fit_3)[[1]],lag.max = 100)
pacf(x=residuals(mod.fit_3)[[1]],lag.max = 100)


AIC_RV <- t(c(AIC(mod.fit),AIC(mod.fit_1),AIC(mod.fit_2),AIC(mod.fit_3)))
BIC_RV <- t(c(BIC(mod.fit),BIC(mod.fit_1),BIC(mod.fit_2),BIC(mod.fit_3)))
IC <- rbind(AIC_RV,BIC_RV);IC

#综合以上结果，可以看出RV-AFRIMA(1,d,1)是其中拟合的最好的模型

##对GARCH类均值模型拟合,根据AIC和BIC判断

plot(tRt1,type='l')

PP.test(tRt1)#  P值显示日对数收益率序列可视为平稳序列
for (k in 1:4) print( Box.test(tRt1, lag = 6*k, type = "Ljung-Box")$p.value)
print(acf(tRt1,lag=24))

x.fit <- arima(tRt1,order = c(0,0,0))
summary(x.fit)
ts.diag(x.fit)
arch.test(x.fit)#通过两个检验的P值都可以看出具有异方差
spec1<-ugarchspec(
  mean.model = list(armaOrder=c(0,0),include.mean=F),
  variance.model = list(garchOrder=c(1,1),model="sGARCH"),
  distribution.model="norm")


##滚动时间窗口预测,用GARCH模型预测

x <- tRt1
spec1<-ugarchspec(
  mean.model = list(armaOrder=c(0,0),include.mean=F),
  variance.model = list(garchOrder=c(1,1),model="sGARCH"),
  distribution.model="std")
fit1 <- ugarchfit(spec1,data = tRt1,methods='ML');fit1

pred_GAR<-numeric(200)
window_size <- 600
for (i in 1:200){
  current_data<-x[i:(i+window_size-2)]
  current_fit <- ugarchfit(spec1,data = current_data,methods='ML')
  current_forecast<-ugarchforecast(current_fit,data=current_data,n.ahead=1)
  pred_GAR[i]<-(current_forecast@forecast$sigmaFor)^2
  print(i)
}
#plot(RV[601:800] ,col='blue')
# lines(pred_GAR,col='red')


##滚动时间窗口预测,用EGARCH模型预测

x <- tRt1
spec2<-ugarchspec(
  mean.model = list(armaOrder=c(0,0),include.mean=F),
  variance.model = list(garchOrder=c(1,1),model="eGARCH"),
  distribution.model="std")
fit2 <- ugarchfit(spec2,data = tRt1,methods='ML');fit2

pred_EGA<-numeric(200)
window_size <- 600
for (i in 1:200){
  current_data<-x[i:(i+window_size-2)]
  current_fit <- ugarchfit(spec2,data = current_data,methods='ML')
  current_forecast<-ugarchforecast(current_fit,data=current_data,n.ahead=1)
  pred_EGA[i]<-(current_forecast@forecast$sigmaFor)^2
  print(i)
}
# plot(RV[601:800] ,col='blue')
# lines(pred_EGA,col='red')


##滚动时间窗口预测,用GJR-GARCH模型预测

x <- tRt1
spec3<-ugarchspec(
  mean.model = list(armaOrder=c(0,0),include.mean=F),
  variance.model = list(garchOrder=c(1,1),model="gjrGARCH"),
  distribution.model="std")
fit3 <- ugarchfit(spec3,data = tRt1,methods='ML');fit3

pred_GJR<-numeric(200)
window_size <- 600
for (i in 1:200){
  current_data<-x[i:(i+window_size-2)]
  current_fit <- ugarchfit(spec3,data = current_data,methods='ML')
  current_forecast<-ugarchforecast(current_fit,data=current_data,n.ahead=1)
  pred_GJR[i]<-(current_forecast@forecast$sigmaFor)^2
  print(i)
}
# plot(RV[601:800] ,col='blue')
# lines(pred_GJR,col='red')

##通过滚动窗口的方法，对已实现波动率进行预测

##RV-ARFIMA模型滚动时间窗口的方法进行预测

x <- RV
pred_RV<-numeric(200)
window_size<-600
for (i in 1:200){
  current_data<-ts(x[i:(i+window_size-1)])
  current_fit <- arfima::arfima(z = current_data,order = c(1,0,1),quiet = TRUE)
  current_forecast<-predict(current_fit,n.ahead = 1)
  pred_RV[i]<-current_forecast[[1]]$Forecast
  print(i)
}

plot(x[601:800] ,col='blue')
lines(pred_RV,col='red')

##InRV-ARFIMA模型滚动时间窗口的方法进行预测

x <- InRV
pred<-numeric(200)
window_size<-600
for (i in 1:200){
  current_data<-ts(x[i:(i+window_size-1)])
  current_fit <- arfima::arfima(z = current_data,order = c(1,0,1),quiet = TRUE)
  #quiet参数可以让他不乱测
  current_forecast<-predict(current_fit,n.ahead = 1)
  pred[i]<-current_forecast[[1]]$Forecast
}
pred_IRV <- exp(pred)
plot(RV[601:800] ,col='blue')
lines(pred_IRV,col='red')


##RV/SD-ARFIMA模型滚动时间窗口的方法进行预测

x <- RVSD
pred_RVSD<-numeric(200)
window_size<-600
for (i in 1:200){
  current_data<-ts(x[i:(i+window_size-1)])
  current_fit <- arfima::arfima(z = current_data,order = c(0,0,0),quiet = TRUE)
  current_forecast<-predict(current_fit,n.ahead = 1)
  pred_RVSD[i]<-current_forecast[[1]]$Forecast
}

plot(RV[601:800] ,col='blue')
lines(pred_RVSD,col='red')


##波动率预测画图 ###图形参数准备

pre_df <- tibble(
  index = 1:length(pred_IRV), 
  pred_IRV = pred_IRV, 
  pred_RV = pred_RV, 
  pred_RVSD = pred_RVSD, 
)
pre_df1 <- tibble(
  index = 1:length(pred_GAR),
  pred_GAR=pred_GAR,
  pred_EGA=pred_EGA,
  pred_GJR=pred_GJR
)

pre_df <-  reshape2::melt(pre_df, id = 'index')
pre_df1 <- reshape2::melt(pre_df1,id='index')

names(pre_df) <-  c('date', 'group', 'y')
names(pre_df1) <- c('date', 'group', 'y')
pre_poi <- tibble(date = 1:length(RV[601:800]), RV = RV[601:800], group = "RV")

###ARFIMA类模型预测图

ggplot(pre_df, aes(x = date, y = y)) + 
  geom_line(aes(group = group, linetype = group, colour = group), 
            size = .8) + 
  geom_point(data = pre_poi, aes(x = date, y = RV, fill = group), 
             size = 1.8, shape = 23, colour = wes_palette('Darjeeling1')[5]) + 
  scale_colour_manual(values = wes_palette('Darjeeling1')[1:3]) + 
  scale_fill_manual(values = wes_palette('Darjeeling1')[5]) + 
  scale_x_continuous(breaks = seq(0, 200, 10)) + 
  scale_y_continuous(breaks = seq(0, 14, 1)) + 
  labs(title = "ARFIMA类模型波动率预测图", subtitle = '2019/07/12-2023/10/28', x = 'time', y = 'volatility') + 
  theme(
    legend.position = c(0.1, 0.8), 
    legend.title = element_blank(),
    plot.subtitle = element_text(hjust = 1), 
    plot.title = element_text(hjust = 0.5), 
  )

###GARCH类模型预测图

ggplot(pre_df1, aes(x = date, y = y)) + 
  geom_line(aes(group = group, linetype = group, colour = group), 
            size = .8) + 
  geom_point(data = pre_poi, aes(x = date, y = RV, fill = group), 
             size = 1.8, shape = 23, colour = wes_palette('Darjeeling1')[5]) + 
  scale_colour_manual(values = wes_palette('Darjeeling1')[1:3]) + 
  scale_fill_manual(values = wes_palette('Darjeeling1')[5]) + 
  scale_x_continuous(breaks = seq(0, 200, 10)) + 
  scale_y_continuous(breaks = seq(0, 14, 1)) + 
  labs(title = "GARCH类模型波动率预测图", subtitle = '2019/07/12-2023/10/28', x = 'time', y = 'volatility') + 
  theme(
    legend.position = c(0.1, 0.8), 
    legend.title = element_blank(),
    plot.subtitle = element_text(hjust = 1), 
    plot.title = element_text(hjust = 0.5), 
  )


##MCS模型预测效果

RV_pred <- RV[601:800]
loc <- c('SE1','SE2','AE1',"AE2","QLIKE","R2LOG")
sta <- c("Tmax","TR")
for (i in loc){
  for (j in sta){
    SE11 <- LossVol(RV_pred,pred_GAR, which = i)
    SE12 <- LossVol(RV_pred,pred_EGA, which = i)
    SE13 <- LossVol(RV_pred,pred_GJR, which = i)
    SE14 <- LossVol(RV_pred,pred_RV,  which = i)
    SE15 <- LossVol(RV_pred,pred_IRV, which = i)
    SE16 <- LossVol(RV_pred,pred_RVSD,which = i)
    df <- data.frame(SE11,SE12,SE13,SE14,SE15,SE16)
    lo <- as.matrix(df)
    mcs_result <- MCSprocedure(Loss=lo,alpha = 0.15,B=8000,cl=NULL,
                               ram.allocation = TRUE,
                               statistic = j,k=NULL,min.k = 5,
                               verbose = F)
    if (j =='Tmax'){
      tmax <- t(t(mcs_result@show[,3]))
      print(sprintf('在 %s 这种方法和在 %s 这种损失函数下的模型评估',j,i))
      print(tmax)
    }
    else {
      tr <-   t(t(mcs_result@show[,6]))
      print(sprintf('在 %s 这种方法和在 %s 这种损失函数下的模型评估',j,i))
      print(tr)
    }
  }
}


#程序每一次跑出来结果都不同，不过是大致准确的


