setwd(Result_path)
Result_vector<-list.files()
load(Result_vector[1])
method_name<-names(funf_result$`different methods`)
loss_fun<-c('MAE','RMSE','HMAE','HRMSE','AMAPE')
lossfun_Result<-array(NA,dim=c(length(loss_fun),length(method_name),length(Result_vector)))# lossfun,method, trajectories
dimnames(lossfun_Result)<-list(loss_fun,method_name,Result_vector)

for(index in 1:length(Result_vector)){
  name<-Result_vector[index]
  load(paste(Result_path,name,sep='\\')) 
  funf_result_new<-funf_result
  
  # MAE
  for(method_index in 1:length(method_name)){
    method=method_name[method_index]
    lossfun_Result[1,method_index,index]<-10000*mean(abs(funf_result$RV_test-funf_result$`different methods`[[method]]))
  }
  
  # RMSE
  for(method_index in 1:length(method_name)){
    method=method_name[method_index]
    lossfun_Result[2,method_index,index]<-10000*sqrt(mean((funf_result$RV_test-funf_result$'different methods'[[method]])^2))
  }
  
  # HMAE
  for(method_index in 1:length(method_name)){
    method=method_name[method_index]
    lossfun_Result[3,method_index,index]<-mean(abs((funf_result$'different methods'[[method]])/(funf_result$RV_test)-1))
  }
  
  # HRMSE
  for(method_index in 1:length(method_name)){
    method=method_name[method_index]
    lossfun_Result[4,method_index,index]<-sqrt(mean(((funf_result$'different methods'[[method]])/(funf_result$RV_test)-1)^2))
  }
  
  # AMAPE
  for(method_index in 1:length(method_name)){
    method=method_name[method_index]
    lossfun_Result[5,method_index,index]<-mean(abs((funf_result$RV_test-funf_result$`different methods`[[method]])/
                                                     (funf_result$RV_test+funf_result$`different methods`[[method]])))
  }
}