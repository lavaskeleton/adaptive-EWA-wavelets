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
library(WaveletComp)
library(signal)




GDP_deOLS=priors$GDP_deOLS


for(i in 1:priors$N)
{
  a=adf.test(GDP_deOLS[,i],k = 1)$p.value
  b=pp.test(GDP_deOLS[,i])$p.value
  c=kpss.test(GDP_deOLS[,i])$p.value
  
  if(a>0.1|b>0.1|c<0.0105)
  {
    print(i)
  }
}
#6,9不平稳
matplot(priors$GDP_deOLS[,c(6,9)],type='l')


#模拟数据的三种滤波方法比较
#压缩效应与相位移动效应

#生成数据
province_number=1
Time=256
u_j=matrix(0,nrow=Time,ncol=1)

for(j in 1:Time)
{
  u_j[j,]=sin(2*pi*j/32)-0.15*sin(2*pi*j/8)
}
#+0.5*sin(2*pi*j/16)
#基于傅立叶变换和Hamming窗函数的滤波方法(IACOBUCCI,NOULLEZ,2005)

U_k=matrix(0,nrow=floor(Time/2)+1,ncol=province_number)
V_k=matrix(0,nrow=floor(Time/2)+1,ncol=province_number)
v_j=matrix(0,nrow=Time,ncol=province_number)
#Fourier变换
for(i in 1:province_number)
{
  for(k in 0:floor(Time/2))
  {
    for(j in 0:(Time-1))
    {
      U_k[k+1,i]=U_k[k+1,i]+u_j[j+1,i]*exp(-1i*2*pi*j*k/Time)
    }
  }
  
  v_l=1/32
  v_h=1/8
  #Hamming窗
  for(k in 0:floor(Time/2))
  {
    if(k>floor(Time*v_l)&&k<ceiling(Time*v_h))
    {
      V_k[k+1,i]=U_k[k+1,i]
    }else if(k==floor(Time*v_l)||k==ceiling(Time*v_h))
    {
      V_k[k+1,i]=(25/46+(1-25/46)/2)*U_k[k+1,i]
    }else if(k==floor(Time*v_l)-1||k==ceiling(Time*v_h)+1)
    {
      V_k[k+1,i]=((1-25/46)/2)*U_k[k+1,i]
    }
  }
  
  
  #重构
  for(j in 0:(Time-1))
  {
    for(k in 1:floor(Time/2))
    {
      v_j[j+1,i]=v_j[j+1,i]+
        (V_k[k+1,i]*exp(1i*2*pi*j*k/Time)+
           Conj(V_k[k+1,i])*exp(-1i*2*pi*j*k/Time))
    }
    v_j[j+1,i]=(V_k[1,i]+v_j[j+1,i])/Time
  }
}






#小波变换


#1/16划分下的尺度（周期）的小波变换
my.data <- data.frame(x = u_j)
my.w <- analyze.wavelet(my.data = my.data, my.series = "x",
                        loess.span = 0,
                        dt = 1, dj = 1/16,
                        lowerPeriod = 1,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)


#时频图
wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 5))


#重构
reconstruct_period=c(8:32)

x=reconstruct(my.w, plot.waves = F,col=c("black","red"),show.legend = F,
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
Scale=matrix(0,nrow=129,ncol=1)

for(i in 0:129)
{
  Scale[i+1,]=6/(2*pi)*2^(i/16)

}

#女儿小波
Wavelet_tao_scale=matrix(0,nrow=256,ncol=129)


for(i in 1:256)
{
  for(j in 1:129)
  {
    for(t in 1:256)
    {
      Wavelet_tao_scale[i,j]=Wavelet_tao_scale[i,j]+u_j[t,]*Conj(mother_morlet((t-i)/Scale[j]))/(Scale[j]^(1/2))
    }
  }
}

#erf函数
erf=function(x)
{
  return(2*pnorm(sqrt(2)*x) - 1)
}


#EWA修正
myWave=matrix(0,nrow=length(c(49:81)),ncol=256)
rec.waves=matrix(0,nrow=256,ncol=length(c(49:81)))
my.w$Period

for(t in 1:256)
{
  for(s in 49:81)
  {     
    
    myWave[s-48,t]=Wavelet_tao_scale[t,s]*(1+erf(t/(sqrt(2)*(Scale[s])))/2)^(-1)*
                   (1+erf((256-t)/(sqrt(2)*(Scale[s])))/2)^(-1)
    
    rec.waves[t,s-48] = (Re(myWave[s-48,t])/sqrt(Scale[s]))*
                        1/16/(pi^(-1/4))* 0.776   

  }     
}







#matplot((t(Mod(Wavelet_tao_scale)^2)/Scale[1:129])[11:129,],type='l')





rec.x=rowSums(rec.waves)
rec.x_mean= mean(rec.x)
rec.x_sd= sd(rec.x)

rec.x= (rec.x- rec.x_mean) * sd(u_j)/rec.x_sd




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
      myWave_ad[s-48,t]=my.w$Wave[s,t]*(1+erf(t/(sqrt(2)*Scale[s]))/2)^(-1)*
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

rec.x_ad= (rec.x_ad- rec.x_ad_mean) * sd(u_j)/rec.x_ad_sd

# par(mfrow=c(2,2)) 
plot(u_j,type='l',main="自适应EWA小波变换",xlab='',ylab='')
lines(rec.x_ad,col='red',lty=2)






plot(my.w$Power.avg[41:89],type='l',xaxt='n')
plot(my.w$Power.avg[41:49],type='l',xaxt='n')

plot(my.w$Power.avg[49:55],type='l',xaxt='n')

plot(my.w$Power.avg[78:84],type='l',xaxt='n')
plot(my.w$Power.avg[78:81],type='l',xaxt='n')

(my.w$Power[42:55,-c(1:7,256:(256-7))])
my.w$Power[42:55,13]

my.w$Period[1:129]
c(49,65,81)


mean((u_j-Re(v_j))^2)

mean((x$series$x-x$series$x.r)^2)

mean((u_j-rec.x)^2)



plot(u_j,type='l',main="窗函数滤波",xlab='',ylab='')
lines(Re(v_j),lty=2,col='red')


plot(x$series$x,type='l',main="小波变换",xlab='',ylab='')
lines(x$series$x.r,lty=2,col='red')


plot(u_j,type='l',main="EWA小波变换",xlab='',ylab='')
lines(rec.x,col='red',lty=2)


