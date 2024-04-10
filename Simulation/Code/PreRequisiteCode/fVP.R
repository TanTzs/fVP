# fVP
setwd(Pre_code_path)
source('kernel_predict_functions.R',encoding = "UTF-8")
source('sample_functions.R',encoding = "UTF-8")

## fVP
stock_vector<-list.files(paste(Data_path))
sorted_indices <- order(sapply(stock_vector,function(x) as.numeric(substr(x, 25, nchar(x)-4))))
stock_vector<-stock_vector[sorted_indices]
fpc_number_vector<-matrix(NA,ncol=3,nrow=length(stock_vector))
fpc_number_vector[,1]<-1:length(stock_vector)
colnames(fpc_number_vector)<-c('index','fpc_number','variance_explained')

for(index in 1:length(stock_vector)){
  stock_code<-stock_vector[index]
  setwd(Data_path)
  # Importing data and pre-computations
  opening_price <- read.table(stock_code, header = F)
  opening_price <- matrix(opening_price[seq(1,nrow(opening_price),by=48),1],nrow=1)[1,1:1500]
  normal_data <- read.table(stock_code, header = T)|>
    as.matrix()|>
    matrix(nrow = 48, byrow = F)
  normal_data <- rbind(opening_price,normal_data)
  Z<-(normal_data)|>apply(2,diff)
  delta<-1/nrow(Z)
  Z<-1/sqrt(delta)*Z
  # Parameter Setting for Predicting
  step_length <- 250 
  window_width <- 1250 
  start=(ncol(normal_data)-step_length-window_width+1)
  
  # Cross-validate Using the First Window
  ob_time<-seq(delta:1,by=delta)
  Z1<-Z[,start:(start+window_width-1)];
  Z1[0>Z1&Z1>=-5*10^(-15)]<-max(Z1[Z1<=-5*10^(-15)])
  Z1[0<=Z1&Z1<=5*10^(-15)]<-min(Z1[Z1>=5*10^(-15)])
  data_cvtest<-log((Z1^2))-q0
  data_cvtest_mean<-apply(data_cvtest,1,mean)
  cv_list<-list()
  for (date in 1:ncol(data_cvtest)){
    Gdate_cvtest<-(data_cvtest[,date]-data_cvtest_mean)%*%
      (t(data_cvtest[,date]-data_cvtest_mean))
    cv_list[[date]]<-Gdate_cvtest
  }
  
  # Determine the Number of Principal Components
  cumulative_variance_explained<-normal_fpc_score(data_cvtest,
                                                  window_width)$variance_explained
  fpc_number<-max(which(cumulative_variance_explained>=
                          0.95)[1])
  fpc_number_vector[,2][index]=fpc_number
  fpc_number_vector[,3][index]<-cumulative_variance_explained[fpc_number]
  
  # Rolling forecast
  V_result_fpckr<-numeric(step_length)
  V_result_fpcVAR<-numeric(step_length)
  
  for(i in 1:step_length){
    print(paste('fVP',stock_code,i))
    #functional principle components scores
    DATA<-Z[,i:(i+window_width-1)];
    DATA[0>DATA&DATA>=-5*10^(-15)]<-max(DATA[DATA<=-5*10^(-15)])
    DATA[0<=DATA&DATA<=5*10^(-15)]<-min(DATA[DATA>=5*10^(-15)])
    
    if(type=='CIRjump'){
      # determine the threshold
      BV<-(pi/2)*apply(abs(DATA[-nrow(DATA),]*DATA[-1,]),2,sum)
      RV<-apply(DATA,2,sum)
      threshold<-3*sqrt(pmin(BV,RV))*delta^(3/8)

      # truncation
      for(col in 1:ncol(Z1)){
        DATA[,col][abs(DATA[,col])>threshold[col]]<-threshold[col]
      }
    }else{
      NULL
    }
    
    DATA<-log(DATA^2)-q0
    
    fpc<-normal_fpc_score(DATA,window_width)
    
    fpc_score<-function(i,k){
      delta*sum(((DATA[,i]-fpc[["1d_mean"]])*fpc[["positive_elements"]]$fpc_fun[,k]))#计算得到FPC分数
    }
    fpc_score<-Vectorize(fpc_score)
    day<-c(1:ncol(DATA));fpc_fun_number<-c(1:fpc_number);
    predictor_fpc_scores<-outer(day,fpc_fun_number,"fpc_score")
    data_learning<-t(predictor_fpc_scores)
    
    #fpcVAR approach
    predict_fpcVAR_scores<-data.frame(t(data_learning))%>%
      VAR_FPC_SCORES(Knumber=fpc_number)

    V_result_fpcVAR[i]<-
      sum(delta*exp(fpc$'1d_mean'+
                      predict_fpcVAR_scores%*%
                      t(as.matrix(fpc$
                                    positive_elements$
                                    fpc_fun[,1:fpc_number],
                                  nrow=fpc_number))))
    
  # # fpckr approach
  #   if(i%%90==1){
  #     hv<-kernel_predictor_CV(data_learning,10,initial_hv=0)$par
  #   }else{
  #     hv=hv
  #   }
  #   predict_fpckr_scores<-kernel_estimator(data_learning,hv)[1:fpc_number]
  # 
  #   V_result_fpckr[i]<-
  #     sum(delta*exp(fpc$'1d_mean'+
  #                     predict_fpckr_scores%*%
  #                     t(as.matrix(fpc$
  #                                   positive_elements$
  #                                   fpc_fun[,1:fpc_number],
  #                                 nrow=fpc_number))))
  
  }
  result_name<-substr(stock_vector[index],
                      1, nchar(stock_vector[index]) - 4)
  
  name0<-paste0(result_name,"_funf_result")
  load(paste(Result_path,name0,sep="\\"))
  # funf_result[["different methods"]][["fpckr"]]<-V_result_fpckr
  funf_result[["different methods"]][["fpcVAR"]]<-V_result_fpcVAR
  save(funf_result,file=paste0(Result_path,'\\',
                               name0))
}


