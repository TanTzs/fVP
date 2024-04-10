# Working Directory
Base_wd<-'D:\\Academic\\Projects\\'  # Can be configured according to the download path
type<-'CIRnojump'  # can choose 'Smooth', 'CIRnojump', or 'CIRjump'

# Automatic settings
Pre_code_path<-paste0(Base_wd,'fVP\\Simulation\\Code\\PreRequisiteCode')
setwd(Pre_code_path)
source('AutoSetting.R',encoding = "UTF-8")

# Forecasting

## Baseline Model
source('BaselineModel.R',encoding='UTF-8')
## fVP
setwd(Pre_code_path)
source('fVP.R',encoding='UTF-8')


# Result Analysis
setwd(Pre_code_path)
source('ResultAnalysis.R',encoding='UTF-8') # Get Array 'lossfun_Result'
##  Array 'lossfun_Result':
## Each dimension represents the loss function, prediction method, and sample trajectory respectively;
MAEmat<-lossfun_Result[1,,] ##  Any dimension can be chosen for observation, for example, observing the MAE loss function.

## Plot
density_est<-apply(MAEmat,1,density)#need sufficient amount of sample trajectory data, please generate more simulated data trajectories
#parameter for plot
library(RColorBrewer)
max_x<-c();max_y<-c();min_x<-c()
for(i in 1:length(density_est)){
  density_test<-density_est[[i]]
  max_x <- c(max_x,max(density_test$x));min_x<-c(min_x,min(density_test$x))
  max_y <- c(max_y,max(density_test$y))
}
max_x<-max(max_x);max_y<-max(max_y);min_x<-min(min_x)
color_density<-brewer.pal(round(length(density_est)),'Set1')
plot(density_est[[1]],main='',
     col=color_density[1],
     xlim=c(min_x, max_x), ylim=c(0, max_y), xlab='MAE', ylab="Density")
lwd<-c()
for(i in 2:length(density_est)){
  if(names(density_est)[i]=='fpcVAR'){
    lines(density_est[[i]],col=color_density[i],lwd=3)
    lwd<-c(lwd,3)
  }else{
    lines(density_est[[i]],col=color_density[i]) 
    lwd<-c(lwd,1)
  }
}
legend("topright", legend=names(density_est), col=color_density, lwd = c(1,lwd),lty=1)

