library(patchwork)
library(MTS)
library(fMultivar)
library(fda)
library(rgl)
library(foreach)
#计算q0
q0f<-function(x){
  log(x)*exp(-x^2/2)
}
q0<-(4/sqrt(2*pi))*integrate(q0f,0,Inf)$value

#Principle components scores
normal_fpc_score<-function(data=normal_Y,window_width){
  #determine the window
  normal_Y<-data
  
  normal_Y_mean<-apply(normal_Y,1,mean)
  
  G_data<-matrix(data=NA,nrow = 
                   nrow(normal_Y)^2,ncol=ncol(normal_Y))
  
  for (date in 1:ncol(normal_Y)){
    Gdate<-(normal_Y[,date]-normal_Y_mean)%*%
      (t(normal_Y[,date]-normal_Y_mean))
    G_data[,date]<-cbind(Gdate)
  }
  
  #covariance matrix
  G_mean<-apply(G_data,1,mean)
  G_mean<-matrix(G_mean,nrow=nrow(normal_Y))
  G_mean[lower.tri(G_mean)] <- t(G_mean)[lower.tri(G_mean)]
  diag(G_mean)<-NA
  
  G_mean_estimate<-G_mean
  for (i in 1:nrow(G_mean_estimate)){
    for (j in 1:ncol(G_mean_estimate)) {
      if (i == j) {
        # 对角线元素替换为附近的值
        if (i > 1) {
          G_mean_estimate[i, j] <- G_mean_estimate[i - 1, j]
        }
        else {
          G_mean_estimate[i, j] <- G_mean_estimate[i, j+1]
        }
      }
    }
  }
  
  G_mean_estimate[lower.tri(G_mean_estimate)] <- t(G_mean_estimate)[lower.tri(G_mean_estimate)]
  
  #eigen functions
  fpc<-eigen(G_mean_estimate)
  positive_index<-which(fpc$values>=0)
  fpc_positive<-list(
    "fpc_value"=fpc$values[positive_index],
    "fpc_fun"=fpc$vectors[,positive_index]*(1/sqrt(delta))
  )
  variance_explained<-cumsum(fpc_positive$fpc_value)/sum(fpc_positive$fpc_value)
  
  fpc_score<-function(i,k){
    delta*sum(((normal_Y[,i]-normal_Y_mean)*fpc_positive$fpc_fun[,k]))
  }
  
  fpc_score<-Vectorize(fpc_score)
  day<-c(1:ncol(normal_Y));fpc_fun_number<-c(1:length(positive_index));
  normal_fpc_score<-outer(day,fpc_fun_number,"fpc_score")
  result<-list("1d_mean"=normal_Y_mean,"2d_mean"=G_mean_estimate,
               "positive_elements"=fpc_positive,"variance_explained"=variance_explained,
               "fpc_scores"=normal_fpc_score)
  return(result)
}

#fpc scores forecasting
VAR_FPC_SCORES<-function(x,Knumber){#x是经过处理的主成分数据,数据类型为matrix，每一行代表该天的主成分得分
  if(Knumber==1){
    # 使用arima()函数拟合ARIMA模型
    arima_model <- auto.arima(x)
    prediction<- forecast(arima_model, h = 1)
    prediction<-prediction$mean[[1]]
  }else{
    VARorder<-VARorder(x,7);#定阶，最高为7
    VARorder<-VARorder$aicor;
    model_fpc_scores<-VAR(x,VARorder);#估计模型
    prediction<-VARpred(model_fpc_scores,1);
    prediction<-prediction$pred
  }
  #利用估计的模型进行一步预测
  return(prediction)
}
