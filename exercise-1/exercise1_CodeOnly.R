source('../functions/sim_functions.R')
knitr::opts_chunk$set(fig.align = 'center')

#Ricker parameters
log.a=1.5 #note - this parameter tends to vary from... ~0.2 to ~ 3 in real populations (use exp(x) to translate this to raw max. recruits/spawner); note at log(0) = 1, this implies essentially an extinction vortex (as all values of S_t > 0 will lead to negative R/S)
smax=5000 #this parameter can vary massively, from hundreds to millions
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50 #time series length

#For this workshop, we've made a function to simulate spawner-recruit dynamics, if you are very keen you can dig into it in the exercise-1/sim_functions.R file.
df.st=salmon_sim(log.a=log.a,smax=smax,sigma=sigma,N=N,form='static')

#visualize your simulated time-series
par(mfrow=c(2,1))
plot(df.st$R~df.st$S,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5),xlim=c(0,max(df.st$S)),ylim=c(0,max(df.st$R)),xlab='spawners',ylab='recruits') #plot data
abline(c(0,1),lty=5) #1:1 line to indicate where recruits = spawners
#expectation based on ricker parameters:
S_p=seq(0,max(df.st$S))
pred=exp(log.a-S_p/smax)*S_p
lines(pred~S_p,lwd=2) #Spawner Recruit curve

#linearized
plot(df.st$logRS~df.st$S,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5),xlim=c(0,max(df.st$S)),xlab='spawners',ylab='log(R/S)') #plot data
pred2=log.a-S_p/smax
lines(pred2~S_p,lwd=2) #Spawner Recruit curve

library(samEst)
#ricker_TMB requires data-frame or list with S = a vector of spawning abundances, and logRS = a vector of log(recruits/spawner)

st.fit1=ricker_TMB(data=df.st,silent=T)#note silent = T just suppresses the updates on the likelihood gradient, the default includes priors - you can set priors_flag = 0 will turn off priors on Ricker parameters

#st.fit1=ricker_stan(data=df.st,ac=F) #alternative Stan variant

#We can view the MLE estimates of the key parameters by calling them up from the function
c(log.a, st.fit1$logalpha)

c(smax, st.fit1$Smax)

c(sigma, st.fit1$sigma)

#compare the estimated curve to the true curve:
plot(df.st$R~df.st$S,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5),xlim=c(0,max(df.st$S)),ylim=c(0,max(df.st$R)),xlab='spawners',ylab='recruits') #plot data
abline(c(0,1),lty=5) #1:1 line to indicate where recruits = spawners
#expectation based on true ricker parameters:
S_p=seq(0,max(df.st$S))
pred=exp(log.a-S_p/smax)*S_p
lines(pred~S_p,lwd=2) #spawner-recruit curve (from simulation)
pred.est=exp(st.fit1$logalpha-S_p/st.fit1$Smax)*S_p
lines(pred.est~S_p,col='darkred',lwd=2) #spawner-recruit curve (predicted)


true.smsy=smsyCalc(log.a,b=1/smax)
c(true.smsy, st.fit1$Smsy) #true and estimated smsy


true.umsy=umsyCalc(log.a)
c(true.umsy, st.fit1$Umsy)


st.fit2=ricker_TMB(data=df.st,silent=T,Smax_mean=10000,Smax_sd=500) #setting prior to be higher for Smax, with low variance

plot(df.st$R~df.st$S,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5),xlim=c(0,max(df.st$S)),ylim=c(0,max(df.st$R)),xlab='spawners',ylab='recruits') #plot data
abline(c(0,1),lty=5) #1:1 line to indicate where recruits = spawners
#expectation based on true ricker parameters:
S_p=seq(0,max(df.st$S))
pred=exp(log.a-S_p/smax)*S_p
lines(pred~S_p,lwd=2) #spawner-recruit curve (simulated relationship)
pred.est1=exp(st.fit1$logalpha-S_p/st.fit1$Smax)*S_p
lines(pred.est1~S_p,col='darkred',lwd=2) #spawner-recruit curve (predicted relationship)
pred.est2=exp(st.fit2$logalpha-S_p/st.fit2$Smax)*S_p
lines(pred.est2~S_p,col='navy',lwd=2) #spawner-recruit curve (new Smax prior)

#Ricker parameters
log.a=1.5 #note - this parameter tends to vary from ~0.2 to ~ 3 in real populations (use exp(x) to translate this to raw max. recruits/spawner); note at log(0) = 1, this implies essentially an extinction vortex (as all values of S_t > 0 will lead to negative expect R/S)
smax=5000 #this parameter can vary massively, from hundreds to millions
sigma=0.6 #most stocks range from ~0.3 to 1.5
rho = 0.9 #autocorrelation parameter - the scale of 'memory' in recruitment deviations

df.ac=salmon_sim(log.a=log.a,smax=smax,sigma=sigma,N=N,rho=rho,form='autocorr')


par(mfrow=c(2,1))
plot(df.ac$R~df.ac$S,bty='l',type='n',xlim=c(0,max(df.ac$S)),ylim=c(0,max(df.ac$R)),xlab='spawners',ylab='recruits') 
lines(df.ac$R~df.ac$S,lwd=0.5,col=adjustcolor('darkgray',alpha.f=0.5))
points(df.ac$R~df.ac$S,pch=21,bg=viridis::viridis(length(df.ac$S))) #plot data
abline(c(0,1),lty=5) #1:1 line to indicate where recruits = spawners

#expectation based on true Ricker parameters:
S_p=seq(0,max(df.ac$S))
pred=exp(log.a-S_p/smax)*S_p
lines(pred~S_p,lwd=2) #spawner-recruit curve

#residual plot
plot(df.ac$eps,bty='l',type='l',xlab='year of simulation',ylab='residual productivity')
abline(h=0,lty=5)
points(df.ac$eps,pch=21,bg=viridis::viridis(length(df.ac$S)))


par(mfrow=c(2,1))
plot(df.st$R~df.st$S,bty='l',type='n',pch=21,xlim=c(0,max(df.st$S)),ylim=c(0,max(df.st$R)),xlab='spawners',ylab='recruits') #plot data
lines(df.st$R~df.st$S,lwd=0.5,col=adjustcolor('darkgray',alpha.f=0.5))
points(df.st$R~df.st$S,pch=21,bg=viridis::viridis(length(df.st$S)))
abline(c(0,1),lty=5) #1:1 line to indicate where recruits = spawners
#expectation based on true ricker parameters:
S_p=seq(0,max(df.st$S))
pred=exp(log.a-S_p/smax)*S_p
lines(pred~S_p,lwd=2) #spawner-recruit curve

#residual plot
plot(df.st$eps,bty='l',type='l',xlab='year of simulation',ylab='residual productivity')
abline(h=0,lty=5)
points(df.st$eps,pch=21,bg=viridis::viridis(length(df.st$S)))

ac.fit.ac=ricker_TMB(data=df.ac,ac=TRUE,silent=T)

#ac.fit.ac=ricker_stan(data=df.st,ac=T) #alternative stan variant

knitr::kable(data.frame("parameter" = c("rho", "Smax", "log_a", "Smsy", "Umsy"), 
                        "true" = c(rho, smax, log.a, smsyCalc(log.a,1/smax), umsyCalc(log.a)), 
                        "estimate" = c(ac.fit.ac$rho, ac.fit.ac$Smax, ac.fit.ac$logalpha, 
                                       ac.fit.ac$Smsy, ac.fit.ac$Umsy)), 
             digits = 2)

static_sr_plot(data=df.ac,mod=ac.fit.ac,plot.params=TRUE) #note - plot.params will print out the main Ricker parameters 

ac.fit.st=ricker_TMB(data=df.ac,ac=FALSE,silent=T)

knitr::kable(data.frame("parameter" = c("rho", "Smax", "log_a", "Smsy", "Umsy"), 
                        "true" = c(rho, smax, log.a, smsyCalc(log.a,1/smax), umsyCalc(log.a)), 
                        "estimate-autcorr" = c(ac.fit.ac$rho, ac.fit.ac$Smax, ac.fit.ac$logalpha, 
                                       ac.fit.ac$Smsy, ac.fit.ac$Umsy),
             "estimate-static" = c(ac.fit.st$rho, ac.fit.st$Smax, ac.fit.st$logalpha, 
                            ac.fit.st$Smsy, ac.fit.st$Umsy)), 
           digits = 2)

st.fit.ac=ricker_TMB(data=df.st,ac=TRUE,silent=T)

st.fit.ac$rho


knitr::kable(data.frame("model" = c("autocorr-sim/autocorrelated-est", "autocorr-sim/stable-est"), 
                        "AICc" = c(ac.fit.ac$AICc, ac.fit.st$AICc), 
                        "BIC" = c(ac.fit.ac$BIC,ac.fit.st$BIC), 
                        "AICc_Weight" = model_weights(c(ac.fit.ac$AICc,ac.fit.st$AICc),form='TMB'), 
                        "BIC_Weight" = model_weights(c(ac.fit.ac$BIC,ac.fit.st$BIC),form='TMB')), 
           digits = 2)

#stable fit with autocorrelation vs. stable fit with static model

knitr::kable(data.frame("model" = c("stable-sim/autocorrelated-est", "stable-sim/stable-est"), 
                        "AICc" = c(st.fit.ac$AICc, st.fit1$AICc), 
                        "BIC" = c(st.fit.ac$BIC,st.fit1$BIC), 
                        "AICc_Weight" = model_weights(c(st.fit.ac$AICc,st.fit1$AICc),form='TMB'), 
                        "BIC_Weight" = model_weights(c(st.fit.ac$BIC,st.fit1$BIC),form='TMB')), 
           digits = 2)


#Ricker parameters
log.a0=1.5 #initial productivity
p.change=-0.75 #proprtional change in time-varying parameter, -0.5 = -50% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50

df.lin.prod=salmon_sim.tv(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par="a",tv.form="linear")

ns.lin.prod.fit=ricker_rw_TMB(data=df.lin.prod,tv.par='a',silent=T)
#ns.lin.prod.fit=ricker_rw_stan(data=df.lin.prod,tv.par='a')

#estimate in red
plot(df.lin.prod$loga.t,type='l',lwd=2,ylab='log(alpha) parameter',xlab='simulation year',ylim=c(min(c(ns.lin.prod.fit$logalpha,df.lin.prod$loga.t)),max(c(ns.lin.prod.fit$logalpha,df.lin.prod$loga.t))))
lines(ns.lin.prod.fit$logalpha,lwd=2,col='darkred') #predicted


rw_sr_plot(data=df.lin.prod,mod=ns.lin.prod.fit)

#Ricker parameters
log.a0=1.5 #static productivity
p.change=-0.75 #proprtional change in time-varying parameter, -0.75 = -75% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50

df.lin.smax=salmon_sim.tv(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='b',tv.form='linear')

ns.lin.smax.fit=ricker_rw_TMB(data=df.lin.smax,tv.par='b',silent=T)
#ns.lin.smax.fit=ricker_rw_stan(data=df.lin.smax,tv.par='b')

plot(df.lin.smax$smax.t,type='l',lwd=2,ylab='Smax parameter',xlab='simulation year',ylim=c(min(c(ns.lin.smax.fit$Smax,df.lin.smax$smax.t)),max(c(ns.lin.smax.fit$Smax,df.lin.smax$smax.t))))
lines(ns.lin.smax.fit$Smax,lwd=2,col='darkred') #predicted

rw_sr_plot(data=df.lin.smax,mod=ns.lin.smax.fit)

#Ricker parameters
log.a0=1.5 #static productivity
p.change=-0.75 #proprtional change in time-varying parameter, -0.75 = -75% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50

iter=50

par(mfrow=c(2,1))
#productivity change
plot(df.lin.prod$loga.t,type='l',lwd=2,ylab='log(alpha) parameter',xlab='simulation year',ylim=c(min(c(ns.lin.prod.fit$logalpha,df.lin.prod$loga.t)),max(c(ns.lin.prod.fit$logalpha,df.lin.prod$loga.t))))
for(i in 1:iter){
df.p=salmon_sim.tv(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par="a",tv.form="linear")
est.p=ricker_rw_TMB(data=df.p,tv.par='a',silent=T)
lines(est.p$logalpha,lwd=2,col=adjustcolor('darkred',alpha.f=0.2)) #predicted
}


#capacity change
plot(df.lin.smax$smax.t,type='l',lwd=2,ylab='Smax parameter',xlab='simulation year',ylim=c(min(c(ns.lin.smax.fit$Smax,df.lin.smax$smax.t)),max(c(ns.lin.smax.fit$Smax,df.lin.smax$smax.t))))
for(i in 1:iter){
df.s=salmon_sim.tv(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par="b",tv.form="linear")
est.s=ricker_rw_TMB(data=df.s,tv.par='b',silent=T)
lines(est.s$Smax,lwd=2,col=adjustcolor('darkred',alpha.f=0.2)) #predicted
}



#Ricker parameters
log.a0=1.5 #initial productivity
p.change=-0.75 #proportional change in time-varying parameter, -0.75 = -75% (on the log-scale)
smax0=5000 #initial Smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50

df.rw.prod=salmon_sim.tv(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='rw')

ns.prod.fit3=ricker_rw_TMB(data=df.rw.prod,tv.par='a',silent=T)
#ns.prod.fit3=ricker_rw_stan(data=df.rw.prod,tv.par='a')

plot(df.rw.prod$loga.t,type='l',lwd=2,ylab='productivity parameter',xlab='simulation year',ylim=c(min(c(df.rw.prod$loga.t,ns.prod.fit3$logalpha)),max(c(df.rw.prod$loga.t,ns.prod.fit3$logalpha))))
lines(ns.prod.fit3$logalpha,lwd=2,col='darkred') #predicted change

#Ricker parameters
log.a0=1.5 #initial productivity
p.change=-0.5#proprtional change in time-varying parameter, -0.5 = -50% (on the log-scale)
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50

df.rw.smax=salmon_sim.tv(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='b',tv.form='rw')

ns.smax.fit3=ricker_rw_TMB(data=df.rw.smax,tv.par='b',silent=T)
#ns.smax.fit3=ricker_rw_stan(data=df.rw.smax,tv.par='b')

plot(df.rw.smax$smax.t,type='l',lwd=2,ylab='Smax parameter',xlab='simulation year',ylim=c(min(c(df.rw.smax$smax.t,ns.smax.fit3$Smax)),max(c(df.rw.smax$smax.t,ns.smax.fit3$Smax))))
lines(ns.smax.fit3$Smax,lwd=2,col='darkred')


ns.lin.smax.prod.fit=ricker_rw_TMB(data=df.lin.smax,tv.par='a',silent=T)

ns.lin.prod.smax.fit=ricker_rw_TMB(data=df.lin.prod,tv.par='b',silent=T)

par(mfrow=c(2,1))
plot(ns.lin.smax.prod.fit$logalpha,type='l',lwd=2,ylab='productivity parameter',xlab='simulation year',col='darkred')
lines(df.lin.smax$loga.t,lwd=2) #true prod change
plot(ns.lin.prod.smax.fit$Smax,type='l',lwd=2,ylab='Smax parameter',xlab='simulation year',col='darkred')
lines(df.lin.prod$smax.t,lwd=2) #true prod change


#tv-productivity generating data - tv productivity and tv smax fit
stable.lin.prod.fit=ricker_TMB(data=df.lin.prod,ac=T,silent=T)

knitr::kable(data.frame("model" = c("sim-prod/est-stable","sim-prod/est-prod", "sim-prod/est-smax"), 
                        "AICc" = c(stable.lin.prod.fit$AICc,ns.lin.prod.fit$AICc, ns.lin.prod.smax.fit$AICc), 
                        "AICc_weight" = model_weights(c(stable.lin.prod.fit$AICc,ns.lin.prod.fit$AICc,ns.lin.prod.smax.fit$AICc),
                                                      form='TMB')), 
           digits = 2)

#tv-smax generating data - tv productivity and tv smax fit
stable.lin.smax.fit=ricker_TMB(data=df.lin.smax,ac=T,silent=T)

knitr::kable(data.frame("model" = c("sim-smax/est-stable","sim-smax/est-prod", "sim-smax/est-smax"), 
                        "AICc" = c(stable.lin.smax.fit$AICc,ns.lin.smax.prod.fit$AICc, ns.lin.smax.fit$AICc), 
                        "AICc_weight" = model_weights(c(stable.lin.smax.fit$AICc,ns.lin.smax.prod.fit$AICc,ns.lin.smax.fit$AICc),form='TMB')), 
           digits = 2)

#Ricker parameters
log.a0=1.5 #initial productivity
p.change=-0.5#proprtional change in time-varying parameter, -0.75 = -75% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50

df.reg.prod=salmon_sim.tv(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='regime',reg.length=8)

hmm.fit.prod=ricker_hmm_TMB(data=df.reg.prod,tv.par='a',k_regime=2,silent=T)


hmm.fit.prod$logalpha

hmm.fit.prod$qij

round(hmm.fit.prod$probregime,2)

hmm.fit.prod$regime

plot(df.reg.prod$loga.t,type='l',lwd=2,ylim=c(min(c(df.reg.prod$loga.t,hmm.fit.prod$logalpha)),max(c(df.reg.prod$loga.t,hmm.fit.prod$logalpha))))
lines(hmm.fit.prod$logalpha[hmm.fit.prod$regime],lwd=2,col='darkred') #predicted state


hmm.fit.prod$logalpha.t=hmm.fit.prod$logalpha%*%hmm.fit.prod$probregime

rw.fit.reg=ricker_rw_TMB(data=df.reg.prod,tv.par='a',silent=T)

plot(df.reg.prod$loga.t,type='l',lwd=2,ylim=c(min(c(df.reg.prod$loga.t,hmm.fit.prod$logalpha.t,rw.fit.reg$logalpha)),max(c(df.reg.prod$loga.t,hmm.fit.prod$logalpha.t,rw.fit.reg$logalpha))),ylab='productivity through time')
lines(as.numeric(hmm.fit.prod$logalpha.t),lwd=2,col='darkred') #predicted state
lines(rw.fit.reg$logalpha,lwd=2,col='navy') #predicted state * probibility 

#Ricker parameters
log.a0=1.5 #initial productivity
p.change=-0.75#proprtional change in productivity, -0.75 = -75% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
p.change2=-0.5#proprtional change in capacity, -0.75 = -75% (on the log-scale), 0.5 = +50%

sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50

df.reg.b=salmon_sim.tv(log.a0=log.a0,p.change=p.change,p.change2=p.change2,smax0=smax0,sigma=sigma,N=N,tv.par='both',tv.form='regime',reg.length=8)

hmm.fit.b=ricker_hmm_TMB(data=df.reg.b,tv.par='both',k_regime=2,silent=T)

hmm_sr_plot(data=df.reg.b,mod=hmm.fit.b)

par(mfrow=c(2,1))
plot(df.reg.b$loga.t,type='l',lwd=2,ylim=c(min(c(df.reg.b$loga.t,hmm.fit.b$logalpha)),max(c(df.reg.b$loga.t,hmm.fit.b$logalpha))),ylab='productivity',xlab='year')
lines(hmm.fit.b$logalpha[hmm.fit.b$regime],lwd=2,col='darkred')
plot(df.reg.b$smax.t,type='l',lwd=2,ylim=c(min(c(df.reg.b$smax.t,hmm.fit.b$Smax)),max(c(df.reg.b$smax.t,hmm.fit.b$Smax))),ylab='Smax',xlab='year')
lines(hmm.fit.b$Smax[hmm.fit.b$regime],lwd=2,col='darkred')

df=read.csv('../datasets/harrison_chinook.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/nicola_chinook.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/kitsumkalum_chinook.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/upper_skeena_chinook.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/middle_yukon_chinook.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/northern_yukon_chinook.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/nordenskiold_chinook.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/chilko_sockeye.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/cultus_sockeye.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/quesnel_sockeye.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/fraser_pinks.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/fraser_canyon_coho.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/mid_fraser_coho.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/south_thompson_coho.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)

df=read.csv('../datasets/north_thompson_coho.csv')
df$logRS=log(df$R/df$S)

f1=ricker_TMB(data=df,ac=T,silent=T)
static_sr_plot(data=df,mod=f1,plot.params=TRUE)

#ricker_rw_TMB(data=df,tv.par='a',silent=T)
#ricker_hmm_TMB(data=df,tv.par='a',silent=T)
