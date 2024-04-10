#本文件储存了需要使用的似然函数与上下界
###估计初始值及上下界###
#GARCH
inig=c(0,0.1,0.1,0.2)
lowerg=c(-Inf,0,0,0)
upperg=c(Inf,Inf,Inf,Inf)

#RGARCH
# ini0=c(0,0.1,0.1,0.2,0.000,0.000,0.00,0.1,1)
#RGARCH
ini0=c(0.0004,0.04,0.001,0.2,0.05,0.000500,0.0004,0.0001,0.01)
lower0=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0)
upper0=c(Inf, Inf ,Inf ,Inf ,Inf ,Inf ,Inf ,Inf,Inf)

#GARCH-MIDAS
K=20
inigm<- c(0,0.1,0.8,0,0.1,0,1.5)#把mbar初值改成了0
lowergm<-c(0,0,0,0,0,0,0)
uppergm<-c(Inf,1,1,Inf,Inf,300,300)

#Heavy
ini_heavy<- c(1,1,0.5)#omega,alpha,beta
lower_heavy<-c(0,0,0)
upper_heavy<-c(Inf,Inf,1)

#MEM
# ini_MEM<-c(0.2,0.3,0.2,0.007,0.1,0.2)
ini_MEM<-c(2,0.04,0.02,0.004,0.1,0.2)
lower_MEM<-c(10^(-7),10^(-7),10^(-7),-Inf,10^(-7),10^(-7))
upper_MEM<-c(Inf,1,1,Inf,Inf,Inf)

################################
#####GARCH(1,1)的Likelihood#####
################################
llfg=function(param,data){#未知波动率
  
  #获取参数
  mu=param[1]
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  t=nrow(data)
  h=numeric(t)
  z=numeric(t)
  
  #计算初值
  # h[1]=1
  # h[1]=omega/(1-beta-gamma)
  h[1]=t^(-1/2)*sum(data[,2][1:floor(t^(1/2))]^2) #omega/(1-beta-gamma)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  
  #计算似然函数值
  QL=0
  llvalue1=0
  if(is.na(beta)&is.na(gamma)){
    QL=-1e10
  }else{
    for (i in 2:t){
      
      #迭代计算
      h[i]=omega+beta*h[i-1]+gamma*(z[i-1]*sqrt(h[i-1]))^2
      z[i]=(data[i,1]-mu)/sqrt(h[i])
      
      #计算似然
      llvalue1=-1/2*log(h[i])-1/2*z[i]^2
      QL=QL+llvalue1
    }
  }
  #返回结果
  return(-QL)
}

geth_llfg=function(param,data){
  #根据估计的参数值和可观测值返回h全时期的值
  
  #获取参数
  mu=param[1]
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  t=nrow(data)
  h=numeric(t)
  z=numeric(t)
  
  #计算h
  # h[1]=1
  h[1]=omega/(1-beta-gamma)
  # h[1]=t^(-1/2)*sum(data[,2][1:floor(t^(1/2))]^2) #omega/(1-beta-gamma)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  for (i in 2:t){
    h[i]=omega+beta*h[i-1]+gamma*(z[i-1]*sqrt(h[i-1]))^2
    z[i]=(data[i,1]-mu)/sqrt(h[i])
  }
  
  #返回结果
  return(h)
}

pllfg=function(mu,data,ht){
  #计算GARCH的似然---需知晓波动率
  
  t=nrow(data)
  h=ht
  z=numeric(t)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  QL=0
  llvalue1=0
  for (i in 2:t){
    z[i]=(data[i,1]-mu)/sqrt(h[i])
    llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(2*pi)
    QL=QL+llvalue1
  }
  return(QL)
}

iterg=function(param,data){
  #根据参数和数据计算下期一期的波动率
  
  n=nrow(data)
  h=numeric(n+1)
  h[1]=1
  # h[1]=param[2]/(1-param[3]-param[4])
  # h[1]=n^(-1/2)*sum(data[,2][1:floor(n^(1/2))]^2) #param[2]/(1-param[3]-param[4]) #取1是不是太武断了？？？
  for (i in 1:n){
    h[i+1]=param[2]+param[3]*h[i]+param[4]*(data[i,1]-param[1])^2
  }
  x_pred = h[n+1]
  return(x_pred)
}



###################################
####Realized GARCH的Likelihood#####
###################################
llf0<-function(param,data){#未知波动率
  
  #获取参数
  mu=param[1]
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  xi=param[5]
  tau1=param[6]
  tau2=param[7]
  sigmau=param[8]
  phi=param[9]
  
  t=nrow(data)
  h=numeric(t)
  z=numeric(t)
  tau=numeric(t)
  g=numeric(t)
  
  #计算初值
  # h[1]=1
  # h[1]=exp((omega+gamma*phi)/(1-beta-gamma*phi))
  h[1]=t^(-1/2)*sum(data[,2][1:floor(t^(1/2))]^2) #exp((omega+gamma*phi)/(1-beta-gamma*phi))
  
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  tau[1]=tau1*z[1]+tau2*(z[1]^2-1)
  g[1]=log(data[1,2])-xi-tau[1]-phi*log(h[1])
  
  #计算似然函数值
  QL=0
  llvalue1=0
  llvalue2=0
  if(is.na(beta)|is.na(gamma)|is.na(phi)){
    QL<-(-1e10)
  }else{
    if(abs(beta+phi*gamma)>=1){
      QL<-(-1e10)
    }else{
      for (i in 2:t){
        #迭代计算
        h[i]=exp(omega+beta*log(h[i-1])+gamma*log(data[i-1,2]))
        z[i]=(data[i,1]-mu)/sqrt(h[i])
        tau[i]=tau1*z[i]+tau2*(z[i]^2-1)
        g[i]=log(data[i,2])-xi-tau[i]-phi*log(h[i])
        
        #计算似然
        llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(sigmau^2)
        llvalue2=-(g[i])^2/(2*sigmau^2)
        QL=QL+llvalue1+llvalue2
      }
    }
  }
  #返回结果
  return(-QL)
}

geth_llf0=function(param,data){
  #根据估计的参数值和可观测值返回h全时期的值
  
  #获取参数
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  phi=param[9]
  t=nrow(data)
  h=numeric(t)
  
  #计算h
  # h[1]=1
  # h[1]=exp((omega+gamma*phi)/(1-beta-gamma*phi))
  h[1]=t^(-1/2)*sum(data[,2][1:floor(t^(1/2))]^2) #exp((omega+gamma*phi)/(1-beta-gamma*phi))
  for (i in 2:t){
    h[i]=exp(omega+beta*log(h[i-1])+gamma*log(data[i-1,2]))
  }
  
  #返回结果
  return(h)
}

pllf0=function(param,data,ht){
  #计算RGARCH的偏似然---需知晓波动率
  
  mu=param[1]
  xi=param[2]
  phi=param[3]
  tau1=param[4]
  tau2=param[5]
  sigmau=param[6]
  t=nrow(data)
  h=ht
  z=numeric(t)
  tau=numeric(t)
  g=numeric(t)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  tau[1]=tau1*z[1]+tau2*(z[1]^2-1)
  g[1]=log(data[1,2])-xi-tau[1]-phi*log(h[1])
  QL=0
  llvalue1=0
  llvalue2=0
  for (i in 2:t){
    z[i]=(data[i,1]-mu)/sqrt(h[i])
    tau[i]=tau1*z[i]+tau2*(z[i]^2-1)
    g[i]=log(data[i,2])-xi-tau[i]-phi*log(h[i])
    llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(sigmau^2)
    llvalue2=-(g[i])^2/(2*sigmau^2)-log(2*pi)
    QL=QL+llvalue1+llvalue2
  }
  return(QL)
}

iter0=function(param,data){
  #根据参数和数据计算下期一期的波动率
  
  n=nrow(data)
  h=numeric(n+1)
  # h[1]=1
  # h[1]=exp((param[2]+param[4]*param[9])/(1-param[3]-param[4]*param[9]))
  h[1]=n^(-1/2)*sum(data[,2][1:floor(n^(1/2))]^2) #exp((param[2]+param[4]*param[9])/(1-param[3]-param[4]*param[9]))
  for (i in 1:n){
    h[i+1]=exp(param[2]+param[3]*log(h[i])+param[4]*log(data[i,2]))
  }
  x_pred = h[n+1]
  return(x_pred)
}

# ####################################
# #####Realized GARCH的Likelihood#####
# ####################################
# llf0=function(param,data){#未知波动率
#   
#   #获取参数
#   mu=param[1]
#   omega=param[2]
#   beta=param[3]
#   gamma=param[4]
#   xi=param[5]
#   tau1=param[6]
#   tau2=param[7]
#   sigmau=param[8]
#   phi=param[9]
#   
#   t=nrow(data)
#   h=numeric(t)
#   z=numeric(t)
#   tau=numeric(t)
#   g=numeric(t)
#   
#   #计算初值
#   # h[1]=1
#   h[1]=exp((omega+gamma*phi)/(1-beta-gamma*phi))
#   z[1]=(data[1,1]-mu)/sqrt(h[1])
#   tau[1]=tau1*z[1]+tau2*(z[1]^2-1)
#   g[1]=log(data[1,2])-xi-tau[1]-phi*log(h[1])
#   
#   #计算似然函数值
#   QL=0
#   llvalue1=0
#   llvalue2=0
#   for (i in 2:t){
#     
#     #迭代计算
#     h[i]=exp(omega+beta*log(h[i-1])+gamma*log(data[i-1,2]))
#     z[i]=(data[i,1]-mu)/sqrt(h[i])
#     tau[i]=tau1*z[i]+tau2*(z[i]^2-1)
#     g[i]=log(data[i,2])-xi-tau[i]-phi*log(h[i])
#     
#     #计算似然
#     llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(sigmau^2)
#     llvalue2=-(g[i])^2/(2*sigmau^2)
#     QL=QL+llvalue1+llvalue2
#   }
#   
#   #返回结果
#   return(-QL)
# }
# 
# geth_llf0=function(param,data){
#   #根据估计的参数值和可观测值返回h全时期的值
#   
#   #获取参数
#   omega=param[2]
#   beta=param[3]
#   gamma=param[4]
#   phi=param[9]
#   t=nrow(data)
#   h=numeric(t)
#   
#   #计算h
#   # h[1]=1
#   h[1]=exp((omega+gamma*phi)/(1-beta-gamma*phi))
#   for (i in 2:t){
#     h[i]=exp(omega+beta*log(h[i-1])+gamma*log(data[i-1,2]))
#   }
#   
#   #返回结果
#   return(h)
# }
# 
# pllf0=function(param,data,ht){
#   #计算RGARCH的偏似然---需知晓波动率
#   
#   mu=param[1]
#   xi=param[2]
#   phi=param[3]
#   tau1=param[4]
#   tau2=param[5]
#   sigmau=param[6]
#   t=nrow(data)
#   h=ht
#   z=numeric(t)
#   tau=numeric(t)
#   g=numeric(t)
#   z[1]=(data[1,1]-mu)/sqrt(h[1])
#   tau[1]=tau1*z[1]+tau2*(z[1]^2-1)
#   g[1]=log(data[1,2])-xi-tau[1]-phi*log(h[1])
#   QL=0
#   llvalue1=0
#   llvalue2=0
#   for (i in 2:t){
#     z[i]=(data[i,1]-mu)/sqrt(h[i])
#     tau[i]=tau1*z[i]+tau2*(z[i]^2-1)
#     g[i]=log(data[i,2])-xi-tau[i]-phi*log(h[i])
#     llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(sigmau^2)
#     llvalue2=-(g[i])^2/(2*sigmau^2)-log(2*pi)
#     QL=QL+llvalue1+llvalue2
#   }
#   return(QL)
# }
# 
# iter0=function(param,data){
#   #根据参数和数据计算下期一期的波动率
#   #获取参数
#   
#   n=nrow(data)
#   h=numeric(n+1)
#   # h[1]=1
#   h[1]=exp((param[2]+param[4]*param[9])/(1-param[3]-param[4]*param[9]))
#   for (i in 1:n){
#     h[i+1]=exp(param[2]+param[3]*log(h[i])+param[4]*log(data[i,2]))
#   }
#   x_pred = h[n+1]
#   return(x_pred)
# }



###################################
#####MDH_RealGARCH的Likelihood#####
###################################

#################################
#####GARCH-MIDAS的likelihood#####
#################################

kernel_phi_vec<-function(K,omega1,omega2){
  
  #需要的核函数
  vec = numeric(K)
  sum = 0
  for (i in 1:K){
    sum = sum + ((i/K)**(omega1-1))*((1-i/K)**(omega2-1))
  }
  for (i in 1:K){
    tmp1 = ((i/K)**(omega1-1))*((1-i/K)**(omega2-1))
    tmp2 = tmp1/sum
    vec[i] = tmp2
  }
  return(rev(vec))
}


llfgm<-function(param,data_5min,data_RV,index0,K){
  #GARCH-MIDAS
  #似然函数本尊
  mu<-param[1]
  alpha = param[2]
  beta = param[3]
  mbar = param[4]
  theta = param[5]
  omega1 = param[6]
  omega2 = param[7]
  vec_tmp = kernel_phi_vec(K,param[6],param[7])
  t_day = length(data_RV)
  g = matrix(0,ncol=t_day,nrow = nrow(data_5min))
  g[1,K+1] = 1#初始化这个即可
  m = numeric(t_day)
  m[K+1] = mbar + theta * vec_tmp%*%data_RV[(index0):(index0-1+K)]#从K+1开始算，
  #最后得到(250-K)个数据
  h = matrix(0,ncol=t_day,nrow = nrow(data_5min))
  h[1,K+1] = g[1,K+1]*m[1]
  QL<-0
  llvalue1<-0
  for(t in (K+1):t_day){#这样给了多少天，会有前K天的数据只用来计算别的
    # t=K+1
    data=data_5min[,t]#取出第t天的数据
    for (i in 1:length(data)){#这个都是固定了天数t得到的结果，对观测点i
      # i=1
      m[t]<- (mbar + theta * vec_tmp%*%data_RV[(index0+t-K-1):(index0-2+t)])
      if(i==1){
        if(t==K+1){
          g[i,t]=g[1,K+1]
        }else{
          g[i,t]<-(1-alpha-beta) + alpha * ((data_5min[length(data),(t-1)]-mu)**2)/m[t] + 
            beta * g[length(data),(t-1)]
        }
      }else{
        g[i,t]<- (1-alpha-beta) + alpha * ((data[i-1]-mu)**2)/m[t] + beta * g[i-1,t]
      }
      h[i,t]<- m[t]*g[i,t]#
      z = (data[i]-mu)/sqrt(h[i,t])
      llvalue1<-llvalue1+(-1/2*log(h[i,t])-1/2*z^2)
    }
    QL<-QL+llvalue1
  }
  return(-QL)
}

penalty_function <- function(param,data_5min,data_RV,index0,K) {
  # 如果 x 不满足约束条件，返回一个大的惩罚值
  if(is.na(param[2])|is.na(param[3])){
    return(1e6)
  }else{
    if ((param[2]+param[3]) >= 1) {
      return(1e6)  # 可以根据具体问题设置适当的惩罚值
    } else {
      return(0)
    }
  }
}

llfgm_with_penalty<-function(param,data_5min,data_RV,index0,K){
  llfgm(param,data_5min,data_RV,index0,K)+ penalty_function(param,data_5min,data_RV,index0,K)
}

itergm<-function(param,data_5min,data_RV,index0,K){
  param<-theta_result
  data_5min<-data_est
  mu<-param[1]
  alpha = param[2]
  beta = param[3]
  mbar = param[4]
  theta = param[5]
  omega1 = param[6]
  omega2 = param[7]
  vec_tmp = kernel_phi_vec(K,param[6],param[7])
  t_day = length(data_RV)
  g = matrix(0,ncol=t_day,nrow = nrow(data_5min))
  g[1,K+1] = 1#初始化这个即可
  m = numeric(t_day+1)
  m[K+1] = mbar + theta * vec_tmp%*%data_RV[(index0):(index0-1+K)]#从K+1开始算，
  #最后得到(250-K)个数据
  h = matrix(0,ncol=t_day,nrow = nrow(data_5min))
  h[1,K+1] = g[1,K+1]*m[1]
  for(t in (K+1):t_day){#这样给了多少天，会有前K天的数据只用来计算别的
    # t=K+1
    data=data_5min[,t]#取出第t天的数据
    for (i in 1:length(data)){#这个都是固定了天数t得到的结果，对观测点i
      # i=2
      m[t]<- (mbar + theta * vec_tmp%*%data_RV[(index0+t-K-1):(index0-2+t)])
      if(i==1){
        if(t==K+1){
          g[i,t]=g[1,K+1]
        }else{
          g[i,t]<-(1-alpha-beta) + alpha * ((data_5min[length(data),(t-1)]-mu)**2)/m[t] + 
            beta * g[length(data),(t-1)]
        }
      }else{
        g[i,t]<- (1-alpha-beta) + alpha * ((data[i-1]-mu)**2)/m[t] + beta * g[i-1,t]
      }
      h[i,t]<- m[t]*g[i,t]#
    }
  }
  #跑完循环后，g矩阵已经全部得到了
  m[t_day+1]<-mbar + theta * vec_tmp%*%data_RV[(index0+t_day-K):(index0-1+t_day)]
  result_day_vol<-m[t_day+1]*(48+(1-((alpha+beta)^48)/(1-(alpha+beta)))*g[48,t_day])
  return(result_day_vol)
}

#################################
#####HEAVY model的likelihood#####
#################################
llf_heavy<-function(param,data,data_RV){#data是日度数据
  T<-length(data_RV);
  omega<-param[1];alpha<-param[2];beta<-param[3];
  rt<-data
  ht<-numeric(T)
  ht[1]<-T^(-1/2)*sum(data[1:floor(T^(1/2))]^2)
  lt<-numeric(T)#用来保存似然值的
  for(t in 2:T){#对天数循环
    ht[t]<-omega+alpha*data_RV[t-1]+beta*ht[t-1]
    lt[t]<-(-1/2*(log(ht[t])+rt[t]^2/ht[t]))
  }
  return(-sum(lt))
}

iter_heavy<-function(param,data,data_RV){
  T<-length(data_RV)
  omega<-param[1];alpha<-param[2];beta<-param[3];
  ht<-numeric(T+1)
  ht[1]<-T^(-1/2)*sum(data[1:floor(T^(1/2))]^2)
  for(t in 1:T){
    ht[t+1]<-omega+alpha*data_RV[t]+beta*ht[t]
  }
  return(ht[T+1])
}

#################################
#####MEM model的likelihood#####
#################################
llf_MEM<-function(param,data,data_square,
                  data_RV,hlt_square){#data是日对数收益率,data_square是日对数收益率平方
  omega<-param[1]
  alpha<-param[2]
  beta<-param[3]
  delta<-param[4]
  phi<-param[5]
  psi<-param[6]
  QL<-0#似然值
  day_num<-length(data_RV)
  ht<-numeric(days_num)
  # ht[1]<-1
  # h[1]=omega/(1-beta-gamma)
  ht[1]<-day_num^(-1/2)*sum(data[1:floor(day_num^(1/2))]^2)
  lt<-numeric(days_num)#用来保存似然值的
  if(is.na(alpha)|is.na(beta)){
    QL<-(-1e10)
  }else{
    if(alpha+beta>=1|(omega-delta^2/(4*alpha))<=0){
      QL<-(-1e10)
    }else{
      for(t in 1:day_num){
        if(t==1){
          ht[t]<-ht[1]
        }else{
          ht[t]<-(omega+alpha*data_square[t-1]+beta*ht[t-1]+delta*data[t-1]+
                    phi*hlt_square[t-1]+psi*data_RV[t-1]) 
        }
        lt[t]<--(log(ht[t])+data_square[t]/ht[t])
      }
      QL<-sum(lt)
    }
  }
  return(-QL)
}

iter_MEM<-function(param,data_square,data,hlt_square,data_RV){
  omega<-param[1]
  alpha<-param[2]
  beta<-param[3]
  delta<-param[4]
  phi<-param[5]
  psi<-param[6]
  day_num<-length(data_RV)
  ht<-numeric(length(data_RV+1))
  # ht[1]<-1
  # h[1]=omega/(1-beta-gamma)
  ht[1]<-day_num^(-1/2)*sum(data[1:floor(day_num^(1/2))]^2)
  for(t in 1:day_num){
    ht[t+1]<-omega+alpha*data_square[t]+beta*ht[t]+
      delta*data[t]+phi*hlt_square[t]+psi*data_RV[t]
  }
  return(ht[day_num+1])
}
