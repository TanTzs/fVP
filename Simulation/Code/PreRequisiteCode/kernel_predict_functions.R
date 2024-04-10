#Kernel predictor

#高精度计算dnorm
library(Rmpfr)
dnorm_mpfr<-function(x){
  x<-mpfr(x,50)
  dnorm<-sqrt(2*pi)^(-1)*exp(-x^2/2)
}
dnorm(c(1,2,4,6))
Vectorize(dnorm_mpfr)

#做一步向前预测
kernel_estimator<-function(data,hv){
  data<-data_learning
  distence_vector<-apply(matrix((data[,-ncol(data)]-data[,ncol(data)])^2,nrow=nrow(data)),2,sum)#代表每一列的距离
  #通过交叉验证选择hv
  weight<-dnorm_mpfr(distence_vector/hv)#得到对数权重
  predict<-rep(0,nrow(data))
  for(i in 1:nrow(data)){
    predict[i]<-as.numeric(sum(data[i,][-1]*weight)/sum(weight))
  }
  prediction<-as.numeric(predict)
  return(prediction)
}

# 做交叉验证
kernel_predictor_CV<-function(data,r,initial_hv){#交叉验证
  #容差设置
  if(initial_hv==0){
    distence_vector<-apply(matrix((data[,-ncol(data)]-data[,ncol(data)])^2,nrow=nrow(data)),2,sum)#代表每一列的距离
    q1 <- quantile(distence_vector, 0.25) # 下四分位数
    q3 <- quantile(distence_vector, 0.75) # 上四分位数
    q<-q3-q1
    initial_hv<-0.9*min(sd(distence_vector),q/1.34)/length(distence_vector)^0.2
  }else{
    initial_hv=initial_hv
  }
  hv_vector<-sort(distence_vector)[1:10]
  error_result<-numeric(length(hv_vector))
  error<-numeric(r)#用于检验的观测曲线个数
  obj<-function(hv){
    for(i in (ncol(data)-r+1):ncol(data)){#从第N-r+1列一直检验到N列
      data_CV<-matrix(data[,1:(i-1)],nrow=nrow(data))#用于预测的历史数据
      error[i-ncol(data)+r]<-sum((data_CV[,ncol(data_CV)]-kernel_estimator(as.matrix(data_CV[,1:ncol(data_CV)],nrow=nrow(data_CV)),hv))^2)#第i个元素代表第i列的误差，0除外
    }
    return(sum(error))
  }
  for(i in 1:length(hv_vector)){
    hv<-hv_vector[i]
    error_result[i]<-obj(hv)
  }
  result<-list(
    'par'=hv_vector[which.min(error_result)],
    'error_result'=error_result
  )
  return(result)
  # return(nlminb(initial_hv,obj,control = list(iter.max = 100,rel.tol=1e-15,abs.tol=1e-15),lower = 0))
}



