library(openxlsx)
library(actuar)
library(psych)
library(MASS)
library(mnormt)
library(mFilter)
library(tseries)
library(urca)
library(forecast)
library(tmvtnorm)
library(zoo)           
library(xts)            
library(plyr)
library(SDMTools) 

library(signal)

library(WaveletComp)


# GDP_deOLS=priors$GDP_deOLS
# 
# 
# for(i in 1:priors$N)
# {
#   a=adf.test(GDP_deOLS[,i],k = 1)$p.value
#   b=pp.test(GDP_deOLS[,i])$p.value
#   c=kpss.test(GDP_deOLS[,i])$p.value
#   
#   if(a>0.1|b>0.1|c<0.0105)
#   {
#     print(i)
#   }
# }
# #6,9不平稳
# matplot(priors$GDP_deOLS[,c(6,9)],type='l')

#需要考虑的问题：第一，主瓣与旁瓣的比值；
#第二，哪些频率段对边际频率段有较大的影响；
#总归一句话，旁瓣衰退速率

#模拟数据的三种滤波方法比较
#压缩效应与相位移动效应

#生成数据
Time=256

simulate_signal=matrix(0,nrow=Time,ncol=1)

for(j in 1:Time)
{
  simulate_signal[j,]=sin(2*pi*j/32)-0.5*sin(2*pi*j/8)
}


#+0.5*sin(2*pi*j/16)
#基于傅立叶变换和Hamming窗函数的滤波方法(IACOBUCCI,NOULLEZ,2005)

Fourier_seq=matrix(0,nrow=floor(Time/2)+1,ncol=1)

Hamming_seq=matrix(0,nrow=floor(Time/2)+1,ncol=1)

Fourier_reconstruct_seq=matrix(0,nrow=Time,ncol=1)

#Fourier变换

  for(k in 0:floor(Time/2))
  {
    for(j in 0:(Time-1))
    {
      Fourier_seq[k+1,]=Fourier_seq[k+1,]+simulate_signal[j+1,]*exp(-1i*2*pi*j*k/Time)
    }
  }
  
  v_l=1/32
  v_h=1/8
  
  #Hamming窗
  for(k in 0:floor(Time/2))
  {
    if(k>floor(Time*v_l)&&k<ceiling(Time*v_h))
    {
      Hamming_seq[k+1,]=Fourier_seq[k+1,]
      
    }else if(k==floor(Time*v_l)||k==ceiling(Time*v_h))
    {
      Hamming_seq[k+1,]=(25/46+(1-25/46)/2)*Fourier_seq[k+1,]
      
    }else if(k==floor(Time*v_l)-1||k==ceiling(Time*v_h)+1)
    {
      Hamming_seq[k+1,]=((1-25/46)/2)*Fourier_seq[k+1,]
    }
  }
  
  
  #重构
  for(j in 0:(Time-1))
  {
    for(k in 1:floor(Time/2))
    {
      Fourier_reconstruct_seq[j+1,]=Fourier_reconstruct_seq[j+1,]+
        (Hamming_seq[k+1,]*exp(1i*2*pi*j*k/Time)+
           Conj(Hamming_seq[k+1,])*exp(-1i*2*pi*j*k/Time))
    }
    Fourier_reconstruct_seq[j+1,]=(Hamming_seq[1,]+Fourier_reconstruct_seq[j+1,])/Time
  }

#小波变换


#1/16划分下的尺度（周期）的小波变换
simulate_signal = data.frame(simulate_signal = simulate_signal)

Wavelet_model = analyze.wavelet(my.data = simulate_signal, my.series = "simulate_signal",
                        loess.span = 0,
                        dt = 1, dj = 1/16,
                        lowerPeriod = 1,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)

simulate_signal=as.matrix(simulate_signal)
#时频图
wt.image(Wavelet_model, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 5))


#重构
reconstruct_period=c(8:32)

Wavelet_reconstruct_seq=reconstruct(Wavelet_model,
                                    plot.waves = F,
                                    col=c("black","red"),
                                    show.legend = F,
                                    legend.coords = "",
                                    sel.lower =reconstruct_period[1],
                                    sel.upper = reconstruct_period[length(reconstruct_period)],
                                    only.sig =F,only.coi = F)



#EWA小波，消除边界效应

#母小波
mother_morlet=function(t)
{
  return(pi^(-1/4)*exp(6i*t)*exp(-t^2/2))
}

#时间尺度（与时间点一致）
Scale=matrix(0,nrow=log2(Time)*16+1,ncol=1)

for(i in 0:(log2(Time)*16+1))
{
  Scale[i+1,]=6/(2*pi)*2^(i/16)

}

#女儿小波
Wavelet_tao_scale=matrix(0,nrow=Time,ncol=log2(Time)*16+1)


for(i in 1:Time)
{
  for(j in 1:(log2(Time)*16+1))
  {
    for(t in 1:Time)
    {
      Wavelet_tao_scale[i,j]=Wavelet_tao_scale[i,j]+simulate_signal[t,]*Conj(mother_morlet((t-i)/Scale[j]))/(Scale[j]^(1/2))
    }
  }
}

#erf函数
erf=function(x)
{
  return(2*pnorm(sqrt(2)*x) - 1)
}


#EWA修正
EWA_wavelet_seq=matrix(0,nrow=length((log2(8)*16+1):(log2(32)*16+1)),ncol=Time)
EWA_reconstruct_wavelet=matrix(0,nrow=Time,ncol=length((log2(8)*16+1):(log2(32)*16+1)))


for(t in 1:Time)
{
  for(s in (log2(8)*16+1):(log2(32)*16+1))
  {     
    
    EWA_wavelet_seq[s-log2(8)*16,t]=Wavelet_tao_scale[t,s]*
                                    (1+erf(t/(sqrt(2)*(Scale[s])))/2)^(-1)*
                                    (1+erf((Time-t)/(sqrt(2)*(Scale[s])))/2)^(-1)
    
    EWA_reconstruct_wavelet[t,s-log2(8)*16] = 1/16/(pi^(-1/4))* 0.776*
                                          (Re(EWA_wavelet_seq[s-log2(8)*16,t])/sqrt(Scale[s]))
                                             

  }     
}







#matplot((t(Mod(Wavelet_tao_scale)^2)/Scale[1:129])[11:129,],type='l')





EWA_reconstruct_seq=rowSums(EWA_reconstruct_wavelet)
EWA_reconstruct_seq_mean= mean(EWA_reconstruct_seq)
EWA_reconstruct_seq_sd= sd(EWA_reconstruct_seq)

EWA_reconstruct_seq= (EWA_reconstruct_seq- EWA_reconstruct_seq_mean) *
                     sd(simulate_signal)/EWA_reconstruct_seq_sd




#自适应EWA小波

#通过边缘周期段的单独重构猜测该周期段的能量
edge_Wave




myWave_ad=matrix(0,nrow=length(c(49:81)),ncol=256)
rec.waves_ad=matrix(0,nrow=256,ncol=length(c(49:81)))



for(t in 1:256)
{
  for(s in 49:81)
  {     

    if(s %in% c(50:57))
    {
      myWave_ad[s-48,t]=Wavelet_model$Wave[s,t]*(1+erf(t/(sqrt(2)*Scale[s]))/2)^(-1)*
        (1+erf((256-t)/(sqrt(2)*Scale[s]))/2)^(-1)
      
    }
    
    else if(s %in% c(56:77))
    {
      myWave_ad[s-48,t]=Wavelet_tao_scale[t,s]*(1+erf(t/(sqrt(2)*Scale[s]))/2)^(-1)*
        (1+erf((256-t)/(sqrt(2)*Scale[s]))/2)^(-1)*0.44/0.65*0.776
      
    }
    else if(s %in% c(78:81))
    {
      myWave_ad[s-48,t]=Wavelet_tao_scale[t,s]*(1+erf(t/(sqrt(2)*Scale[s]))/2)^(-1)*
        (1+erf((256-t)/(sqrt(2)*Scale[s]))/2)^(-1)*0.44/0.65
    }
    
    
    rec.waves_ad[t,s-48] = (Re(myWave[s-48,t])/sqrt(Scale[s]))*
      1/16/(pi^(-1/4))* 0.776
    
  }
}



rec.x_ad=rowSums(rec.waves_ad)
rec.x_ad_mean= mean(rec.x_ad)
rec.x_ad_sd= sd(rec.x_ad)

rec.x_ad= (rec.x_ad- rec.x_ad_mean) * sd(simulate_signal)/rec.x_ad_sd

# par(mfrow=c(2,2)) 
plot(simulate_signal,type='l',main="自适应EWA小波变换",xlab='',ylab='')
lines(rec.x_ad,col='red',lty=2)






plot(Wavelet_model$Power.avg[41:89],type='l',xaxt='n')
plot(Wavelet_model$Power.avg[41:49],type='l',xaxt='n')

plot(Wavelet_model$Power.avg[49:55],type='l',xaxt='n')

plot(Wavelet_model$Power.avg[78:84],type='l',xaxt='n')
plot(Wavelet_model$Power.avg[78:81],type='l',xaxt='n')

(Wavelet_model$Power[42:55,-c(1:7,256:(256-7))])
Wavelet_model$Power[42:55,13]

Wavelet_model$Period[1:log2(Time)*16+1]
c(49,65,81)


mean((simulate_signal-Re(Fourier_reconstruct_seq))^2)


mean((simulate_signal-Wavelet_reconstruct_seq$series$simulate_signal.r)^2)

mean((simulate_signal-EWA_reconstruct_seq)^2)



plot(simulate_signal,type='l',main="窗函数滤波",xlab='',ylab='')
lines(Re(Fourier_reconstruct_seq),lty=2,col='red')


plot(simulate_signal,type='l',main="小波变换",xlab='',ylab='')
lines(Wavelet_reconstruct_seq$series$simulate_signal.r,lty=2,col='red')


plot(simulate_signal,type='l',main="EWA小波变换",xlab='',ylab='')
lines(EWA_reconstruct_seq,col='red',lty=2)


