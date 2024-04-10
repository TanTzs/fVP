# Auto-Setting
if(type=='Smooth'){
  IV_path<-paste0(Base_wd,'fVP\\Simulation\\Data\\Smooth\\IV')
  Data_path<-paste0(Base_wd,'fVP\\Simulation\\Data\\Smooth\\Price')
  Pre_code_path<-paste0(Base_wd,'fVP\\Simulation\\Code\\PreRequisiteCode')
  Result_path<-paste0(Base_wd,'fVP\\Simulation\\','Result\\Smooth')
}else if(type=='CIRnojump'){
  IV_path<-paste0(Base_wd,'fVP\\Simulation\\Data\\CIR\\IV')
  Data_path<-paste0(Base_wd,'fVP\\Simulation\\Data\\CIR\\Price')
  Pre_code_path<-paste0(Base_wd,'fVP\\Simulation\\Code\\PreRequisiteCode')
  Result_path<-paste0(Base_wd,'fVP\\Simulation\\','Result\\CIRnojump')
}else if(type=='CIRjump'){
  IV_path<-paste0(Base_wd,'fVP\\Simulation\\Data\\CIR\\IV')
  Data_path<-paste0(Base_wd,'fVP\\Simulation\\Data\\CIR\\Jump_price')
  Pre_code_path<-paste0(Base_wd,'fVP\\Simulation\\Code\\PreRequisiteCode')
  Result_path<-paste0(Base_wd,'fVP\\Simulation\\','Result\\CIRjump')
}

# Packages Setting
library(gmm)
library(dplyr)
library(KernSmooth)
library(forecast)