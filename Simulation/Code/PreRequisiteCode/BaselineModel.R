# BaseCIRnojump
setwd(Pre_code_path)
source('Likelihood_function.R',encoding = 'UTF-8')

## Determine the stock code and index code for forecasting.
Stock_vector<-list.files(Data_path)
sorted_indices <- order(sapply(Stock_vector,function(x) as.numeric(substr(x, 25, nchar(x)-4))))
Stock_vector<-Stock_vector[sorted_indices]
IV_vector<-list.files(IV_path)
sorted_indices <- order(sapply(IV_vector,function(x) as.numeric(substr(x, 12, nchar(x)-4))))
IV_vector<-IV_vector[sorted_indices]

## Predicting based on the baseline model
for(index in 1:length(Stock_vector)){
  stock_code<-Stock_vector[index]
  IV_code<-IV_vector[index]
  
  # Importing data and pre-computations
  setwd(Data_path)
  opening_price <- read.table(stock_code, header = F)
  opening_price <- matrix(opening_price[seq(1,nrow(opening_price),by=48),1],nrow=1)[1,1:1500]
  normal_data <- read.table(stock_code, header = T)|>
    as.matrix()|>
    matrix(nrow = 48, byrow = F)
  setwd(IV_path)
  IV<-read.table(IV_code,header=F)[,1]
  normal_data <- rbind(opening_price,normal_data)
  normal_data_range<-apply(normal_data,2,max)-apply(normal_data,2,min)
  normal_r2<-apply(normal_data,2,diff)^2
  RV<-(apply(normal_r2,2,sum))

  # Parameter Setting for Predicting
  step_length <- 250 
  window_width <- 1250 
  start=(ncol(normal_data)-step_length-window_width+1)
  delta<-1/nrow(normal_data)

  fore_num<-step_length
  Lw=window_width #用于估计模型的时间长度
  days_num=Lw+fore_num #用于估计的数据和用于对比的数据总共的窗宽
  
  logreturn_5min<-apply(normal_data,2,diff)
  logreturn_daily<-unlist(as.vector((normal_data[49,])-normal_data[1,]))
  RV_daily<-tail(RV,days_num)
  logreturn_daily_square<-logreturn_daily^2
  hlt_square_ori<-tail(normal_data_range^2,days_num)
  data_daily<-data.frame(
    "daily_Return"=logreturn_daily,
    "daily_RV"=RV_daily
  )
  RV_test<-tail(RV,fore_num)
  IV_test<-tail(IV,fore_num)
  
  # Predicting
  
  ## GARCH Model
  forecast_GARCH=numeric(fore_num)
  for (day_id in 1:fore_num){
    data_est=as.matrix(data_daily[day_id:(day_id+Lw-1),])
    fit_RV=nlminb(inig, llfg, data=data_est, lower=lowerg, upper=upperg)
    theta=fit_RV[["par"]]
    forecast_GARCH[day_id]=iterg(theta,data_est)
    print(c(index,'GARCH',day_id))
  }
  
  ## RGARCH Model
  forecast_RGARCH=numeric(fore_num)
  for (day_id in 1:fore_num){
    data_est=as.matrix(data_daily[day_id:(day_id+Lw-1),]) 
    fit_RV=nlminb(ini0, llf0, data=data_est, lower=lower0, upper=upper0)
    theta=fit_RV[["par"]]
    forecast_RGARCH[day_id]=iter0(theta,data_est)
    print(c(index,'RGARCH',day_id))
  }
  
  ## MEM Model
  forecast_MEM=numeric(fore_num)
  for (day_id in 1:fore_num){
    data=logreturn_daily[(day_id):(day_id-1+Lw)]
    data_square=logreturn_daily_square[day_id:(day_id-1+Lw)]
    data_RV=RV_daily[day_id:(day_id+Lw-1)]
    hlt_square=hlt_square_ori[day_id:(day_id-1+Lw)]
    
    fit_RV=nlminb(ini_MEM,llf_MEM,data=data,data_square=data_square,
                  lower=lower_MEM, upper=upper_MEM,
                  data_RV=data_RV,hlt_square=hlt_square)
    theta_result=fit_RV[["par"]]
    
    forecast_MEM[day_id]=iter_MEM(theta_result,data_square=data_square,
                                  data=data,hlt_square = hlt_square,
                                  data_RV=data_RV)#
    print(c(index,'MEM',day_id))
  }
  
  ## HAR Model
  forecast_HAR=numeric(fore_num)
  for (day_id in 1:fore_num){
    data<-data_daily
    data_HAR=matrix(0,Lw,4)
    for (k in 22:Lw){
      data_HAR[k,1]=data[(k+day_id-1),2]
      data_HAR[k,2]=data[(k+day_id-2),2]
      data_HAR[k,3]=mean(data[(k+day_id-6):(k+day_id-3),2])
      data_HAR[k,4]=mean(data[(k+day_id-23):(k+day_id-7),2])
    }
    coef=lm(data_HAR[,1]~data_HAR[,2:4])[["coefficients"]]
    forecast_HAR[day_id]=coef[1]+coef[2]*data_HAR[Lw,2]+coef[3]*data_HAR[Lw,3]+coef[4]*data_HAR[Lw,4]
    print(c(index,'HAR',day_id))
  }

  ## HEAVY Model
  forecast_heavy=numeric(fore_num)
  for (day_id in 1:fore_num){
    data_est=logreturn_daily[day_id:(day_id-1+Lw)]
    data_RV=RV_daily[day_id:(day_id+Lw-1)]
    fit_RV=nlminb(ini_heavy, llf_heavy, data=data_est,
                  lower=lower_heavy, upper=upper_heavy,data_RV=data_RV)
    theta_result=fit_RV[["par"]]
    forecast_heavy[day_id]=iter_heavy(theta_result,data_est,data_RV)#
    print(c(index,'HEAVY',day_id))
  }
  
  # Result Save
  funf_result<-list()
  funf_result[['RV_test']]<-IV_test
  funf_result[["different methods"]][["GARCH"]]<-forecast_GARCH
  funf_result[["different methods"]][["RGARCH"]]<-forecast_RGARCH
  funf_result[["different methods"]][["HAR_RV"]]<-forecast_HAR
  funf_result[["different methods"]][["HEAVY"]]<-forecast_heavy
  funf_result[["different methods"]][["MEM"]]<-forecast_MEM

  save(funf_result,file =
         paste0(Result_path,'\\',
                substr(Stock_vector[index],1,nchar(Stock_vector[index])-4),
                "_funf_result"))
}

