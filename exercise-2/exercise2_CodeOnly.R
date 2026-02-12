source('../functions/sim_functions.R')
source('../functions/plot_functions.R')
knitr::opts_chunk$set(fig.align = 'center')

library(samEst)
log.a=seq(0.1,2,by=0.1)
smax0=5000
log.a0=1
smax=seq(500,500000)
smsy.a=smsyCalc(log.a,1/smax0)
smsy.b=smsyCalc(log.a0,1/smax)
umsy=umsyCalc(log.a)

par(mfrow=c(2,2))
plot(smsy.a~log.a,type='l',lwd=2,bty='l',ylab ='Smsy')
plot(smsy.b~smax,type='l',lwd=2,bty='l',ylab ='Smsy', xlab ='Smax')
plot(umsy~log.a,type='l',lwd=2,bty='l',ylim=c(0,1),ylab='Umsy')
plot(rep(umsyCalc(log.a0),length(smax))~smax,type='l',lwd=2,bty='l',ylab='Umsy', xlab='Smax',ylim=c(0,1))

#Ricker parameters
log.a0=1.2 #initial productivity
p.change=-0.85 #proprtional change in time-varying parameter, -0.5 = -50% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50

df.lin.prod=salmon_sim.tv(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='linear')

ns.prod.fit1=ricker_rw_TMB(data=df.lin.prod,tv.par='a',silent=T)
#ns.prod.fit1=ricker_rw_stan(data=df.lin.prod,tv.par='a')
true.umsy=umsyCalc(df.lin.prod$loga.t)
true.smsy=smsyCalc(df.lin.prod$loga.t,1/df.lin.prod$smax.t)

#estimate in red
par(mfrow=c(2,1))
plot(true.umsy,type='l',lwd=2,ylab='Umsy',xlab='simulation year',ylim=c(min(c(true.umsy,ns.prod.fit1$Umsy)),max(c(true.umsy,ns.prod.fit1$Umsy))))
lines(ns.prod.fit1$Umsy,lwd=2,col='darkred')
plot(true.smsy,type='l',lwd=2,ylab='Smsy',xlab='simulation year',ylim=c(min(c(true.smsy,ns.prod.fit1$Smsy)),max(c(true.smsy,ns.prod.fit1$Smsy))))
lines(ns.prod.fit1$Smsy,lwd=2,col='darkred')

#Ricker parameters
log.a=c(0.5,1,1.5)
smax0=c(5000)
smsy=smsyCalc(log.a,1/smax0)
umsy=umsyCalc(log.a)

#Harvest control rule parameters
U.min=0.025
eg.scalar=1
upper.tar.scalar=2
U.scalar=1

hcr_plot_example(U.min=U.min,eg.limit=eg.scalar*smsy,upper.tar=upper.tar.scalar*smsy,target.ER=U.scalar*umsy,smsy=smsy,smax0=smax0)


log.a0=1.5 #initial productivity
p.change=-0.75 #proprtional change in time-varying parameter, -0.75 = -75% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50 #length of simulation

#common parameters between different assessment model forms to keep comparability:

eg.scalar=1 #escapement goal as a fraction of Smsy - ie. 1 = 1*Smsy
assess.freq=10 #years between reference point re-assessments
upper.tar.scalar=2 #run size, as a fraction of Smsy, to maximize harvest rate - ie. 2 = 2*Smsy
U.scalar=1 #max. harvest target, as a fraction of Umsy
U.min=0.025 #minimum harvest, applied below escapement goal


mse.prod.st=salmon_sim.tv_hcr(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='linear',hcr.form='stable',hcr.par=NULL,assess.freq=assess.freq,eg.scalar=eg.scalar,upper.tar.scalar=upper.tar.scalar,U.scalar=U.scalar,U.min=U.min)

mse.prod.ns=salmon_sim.tv_hcr(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='linear',hcr.form='rw',hcr.par='a',assess.freq=assess.freq,eg.scalar=eg.scalar,upper.tar.scalar=upper.tar.scalar,U.scalar=U.scalar,U.min=U.min)


par(mfrow=c(2,1))
plot(mse.prod.st$umsy.true,ylim=c(0,1),type='l',lwd=2,xlab='year of simulation',ylab='Umsy')
lines(mse.prod.st$umsy.est,lwd=2,col='navy')
lines(mse.prod.ns$umsy.est,lwd=2,col='darkred')
text(x=par('usr')[2]*0.1,y=par('usr')[4]*0.15,'Stable',col='navy')
text(x=par('usr')[2]*0.1,y=par('usr')[4]*0.1,'Non-stationary',col='darkred')

plot(mse.prod.st$smsy.true,ylim=c(min(c(mse.prod.st$smsy.true,mse.prod.st$smsy.est,mse.prod.ns$smsy.est)),max(c(mse.prod.st$smsy.true,mse.prod.st$smsy.est,mse.prod.ns$smsy.est))),type='l',lwd=2,xlab='year of simulation',ylab='Smsy')
lines(mse.prod.st$smsy.est,lwd=2,col='navy')
lines(mse.prod.ns$smsy.est,lwd=2,col='darkred')
text(x=par('usr')[2]*0.1,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.8,'Stable',col='navy')
text(x=par('usr')[2]*0.1,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.85,'Non-stationary',col='darkred')

tp=c(1,round(nrow(mse.prod.st)/2),nrow(mse.prod.st)) #time slice at t=1, half way, final

eg.limit.st=eg.scalar*mse.prod.st$smsy.est #escapement goals through time for stable simulation
eg.limit.ns=eg.scalar*mse.prod.ns$smsy.est #escapement goals through time for non-stationary simulation
upper.tar.st=upper.tar.scalar*mse.prod.st$smsy.est #upper limit beyond which harvest is static at Umsy
upper.tar.ns=upper.tar.scalar*mse.prod.ns$smsy.est #upper limit beyond which harvest is static at Umsy

target.ER.st=U.scalar*mse.prod.st$umsy.est #max. harvest rate
target.ER.ns=U.scalar*mse.prod.ns$umsy.est #max. harvest rate

par(mfrow=c(2,1))
#stable Ricker model plot
hcr_plot_est(sim.df=mse.prod.st,tp=tp,eg.limit=eg.limit.st,upper.tar=upper.tar.st,target.ER=target.ER.st,U.min=U.min,title='Stable Ricker model')
#nonstationary Ricker model plot
hcr_plot_est(sim.df=mse.prod.ns,tp=tp,eg.limit=eg.limit.ns,upper.tar=upper.tar.ns,target.ER=target.ER.ns,U.min=U.min,title='Non-stationary Ricker model')

#scaled -catch 
scaled.catch.st=mean(mse.prod.st$catch/max(c(mse.prod.st$catch,mse.prod.ns$catch)))
scaled.catch.ns=mean(mse.prod.ns$catch/max(c(mse.prod.st$catch,mse.prod.ns$catch)))

#proprtion of time spawners exceed escapement goal
prop.spn.st=sum(ifelse(mse.prod.st$S>mse.prod.st$smsy.true*eg.scalar,1,0))/nrow(mse.prod.st)
prop.spn.ns=sum(ifelse(mse.prod.ns$S>mse.prod.ns$smsy.true*eg.scalar,1,0))/nrow(mse.prod.ns)


plot(c(0,1)~c(0.6,2.4),type='n',ylim=c(0,1),xaxt='n',ylab='Scaled to max. observed (C) or Proportion of years (S)',xlab='')
lines(c(scaled.catch.st,prop.spn.st)~c(0.95,1.95),lwd=2,col='navy')
points(c(scaled.catch.st,prop.spn.st)~c(0.95,1.95),pch=21,bg='navy',cex=2)
lines(c(scaled.catch.ns,prop.spn.ns)~c(1.05,2.05),lwd=2,col='darkred')
points(c(scaled.catch.ns,prop.spn.ns)~c(1.05,2.05),pch=21,bg='darkred',cex=2)
mtext(side=1,'Scaled Catch',at=1,line=1)
mtext(side=1,'Prop. years above EG',at=2,line=1)
text(x=par('usr')[1]*1.5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,'Stable',col='navy')
text(x=par('usr')[1]*1.5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.2,'Non-stationary',col='darkred')

log.a0=1.5 #initial productivity
p.change=-0.75 #proprtional change in time-varying parameter, -0.75 = -75% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50 #length of simulation

#common parameters between different assessment model forms to keep comparability:

eg.scalar=1 #escapement goal as a fraction of Smsy - ie. 1 = 1*Smsy
assess.freq=10 #years between reference point re-assessments
upper.tar.scalar=2 #run size, as a fraction of Smsy, to maximize harvest rate - ie. 2 = 2*Smsy
U.scalar=1 #max. harvest target, as a fraction of Umsy
U.min=0.025 #minimum harvest, applied below escapement goal

sims.st=list() #list to hold outputs for simulations with stable assessment models
sims.ns=list() #list to hold outputs for simulations with stable assessment models

iter=50 #nynber of simulations to run
m.catch.st=numeric(iter);m.catch.ns=numeric(iter)
sd.catch=numeric(iter);sd.catch.ns=numeric(iter)
m.spn.st=numeric(iter);m.spn.ns=numeric(iter)
sd.spn.st=numeric(iter);sd.spn.ns=numeric(iter)

for(i in 1:iter){
  sims.st[[i]]=salmon_sim.tv_hcr(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='linear',hcr.form='stable',hcr.par=NULL,assess.freq=assess.freq,eg.scalar=eg.scalar,upper.tar.scalar=upper.tar.scalar,U.scalar=U.scalar,U.min=U.min)
  
  sims.ns[[i]]=salmon_sim.tv_hcr(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='linear',hcr.form='rw',hcr.par='a',assess.freq=assess.freq,eg.scalar=eg.scalar,upper.tar.scalar=upper.tar.scalar,U.scalar=U.scalar,U.min=U.min)


  #calculate sim specific mean & sd of catch and spawner abundance
   m.spn.st[i]= sum(ifelse(sims.st[[i]]$S>sims.st[[i]]$smsy.true*eg.scalar,1,0))/nrow(sims.st[[i]])
  m.spn.ns[i]= sum(ifelse(sims.ns[[i]]$S>sims.ns[[i]]$smsy.true*eg.scalar,1,0))/nrow(sims.ns[[i]])
 
}
m.catch.st=as.numeric(lapply(sims.st, function(x) exp(mean(log(x$catch))))) #geometric mean catch per simulation
m.catch.ns=as.numeric(lapply(sims.ns, function(x) exp(mean(log(x$catch)))))
sc.catch.st=m.catch.st/max(c(m.catch.st,m.catch.ns))
sc.catch.ns=m.catch.ns/max(c(m.catch.st,m.catch.ns))


plot(c(0,1)~c(0.6,2.4),type='n',ylim=c(0,1),xaxt='n',ylab='Scaled to max. observed (C) or Proportion of years (S)',xlab='')
lines(c(mean(sc.catch.st)-sd(sc.catch.st),mean(sc.catch.st)+sd(sc.catch.st))~rep(0.95,2),col='navy')
lines(c(mean(m.spn.st)-sd(m.spn.st),mean(m.spn.st)+sd(m.spn.st))~rep(1.95,2),col='navy')
lines(c(mean(sc.catch.ns)-sd(sc.catch.ns),mean(sc.catch.ns)+sd(sc.catch.ns))~rep(1.05,2),col='darkred')
lines(c(mean(m.spn.ns)-sd(m.spn.ns),mean(m.spn.ns)+sd(m.spn.ns))~rep(2.05,2),col='darkred')
lines(c(mean(sc.catch.st),mean(m.spn.st))~c(0.95,1.95),lwd=2,col='navy')
points(c(mean(sc.catch.st),mean(m.spn.st))~c(0.95,1.95),pch=21,bg='navy',cex=2)
lines(c(mean(sc.catch.ns),mean(m.spn.ns))~c(1.05,2.05),lwd=2,col='darkred')
points(c(mean(sc.catch.ns),mean(m.spn.ns))~c(1.05,2.05),pch=21,bg='darkred',cex=2)
mtext(side=1,'Scaled Catch',at=1,line=1)
mtext(side=1,'Prop. years above EG',at=2,line=1)
text(x=par('usr')[1]*1.5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,'Stable',col='navy')
text(x=par('usr')[1]*1.5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.2,'Non-stationary',col='darkred')

log.a0=1.5 #initial productivity
p.change=-0.75 #proprtional change in time-varying parameter, -0.75 = -75% (on the log-scale), 0.5 = +50%
smax0=5000 #initial smax
sigma=0.6 #most stocks range from ~0.3 to 1.5
N=50 #length of simulation

#common parameters between different assessment model forms to keep comparability:

eg.scalar=1 #escapement goal as a fraction of Smsy - ie. 1 = 1*Smsy
assess.freq=10 #years between reference point re-assessments
upper.tar.scalar=2 #run size, as a fraction of Smsy, to maximize harvest rate - ie. 2 = 2*Smsy
U.scalar=1 #max. harvest target, as a fraction of Umsy
U.min=0.025 #minimum harvest, applied below escapement goal

sims.st=list() #list to hold outputs for simulations with stable assessment models
sims.ns=list() #list to hold outputs for simulations with stable assessment models
sims.mix=list() #list to hold outputs for simulations with stable assessment models

iter=50 #nynber of simulations to run
m.catch.st=numeric(iter);m.catch.ns=numeric(iter);m.catch.mix=numeric(iter)
m.spn.st=numeric(iter);m.spn.ns=numeric(iter);m.spn.mix=numeric(iter)

for(i in 1:iter){
  sims.st[[i]]=salmon_sim.tv_hcr(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='linear',hcr.form='stable',hcr.par=NULL,assess.freq=assess.freq,eg.scalar=eg.scalar,upper.tar.scalar=upper.tar.scalar,U.scalar=U.scalar,U.min=U.min)
  
  sims.ns[[i]]=salmon_sim.tv_hcr(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='linear',hcr.form='rw',hcr.par='a',assess.freq=assess.freq,eg.scalar=eg.scalar,upper.tar.scalar=upper.tar.scalar,U.scalar=U.scalar,U.min=U.min)
  
  sims.mix[[i]]=salmon_sim.tv_hcr(log.a0=log.a0,p.change=p.change,smax0=smax0,sigma=sigma,N=N,tv.par='a',tv.form='linear',hcr.form='mixed-rw',hcr.par='a',assess.freq=assess.freq,eg.scalar=eg.scalar,upper.tar.scalar=upper.tar.scalar,U.scalar=U.scalar,U.min=U.min)


  #calculate sim specific mean & sd of catch and spawner abundance
  m.spn.st[i]= sum(ifelse(sims.st[[i]]$S>sims.st[[i]]$smsy.true*eg.scalar,1,0))/nrow(sims.st[[i]])
  m.spn.ns[i]= sum(ifelse(sims.ns[[i]]$S>sims.ns[[i]]$smsy.true*eg.scalar,1,0))/nrow(sims.ns[[i]])
  m.spn.mix[i]= sum(ifelse(sims.mix[[i]]$S>sims.mix[[i]]$smsy.true*eg.scalar,1,0))/nrow(sims.mix[[i]])
 
}
m.catch.st=as.numeric(lapply(sims.st, function(x) exp(mean(log(x$catch))))) #geometric mean catch per simulation
m.catch.ns=as.numeric(lapply(sims.ns, function(x) exp(mean(log(x$catch)))))
m.catch.mix=as.numeric(lapply(sims.mix, function(x) exp(mean(log(x$catch)))))
sc.catch.st=m.catch.st/max(c(m.catch.st,m.catch.ns,m.catch.mix))
sc.catch.ns=m.catch.ns/max(c(m.catch.st,m.catch.ns,m.catch.mix))
sc.catch.mix=m.catch.mix/max(c(m.catch.st,m.catch.ns,m.catch.mix))

plot(c(0,1)~c(0.6,2.4),type='n',ylim=c(0,1),xaxt='n',ylab='Scaled to max. observed (C) or Proportion of years (S)',xlab='')
lines(c(mean(sc.catch.st)-sd(sc.catch.st),mean(sc.catch.st)+sd(sc.catch.st))~rep(0.95,2),col='navy')
lines(c(mean(m.spn.st)-sd(m.spn.st),mean(m.spn.st)+sd(m.spn.st))~rep(1.95,2),col='navy')
lines(c(mean(sc.catch.ns)-sd(sc.catch.ns),mean(sc.catch.ns)+sd(sc.catch.ns))~rep(1.05,2),col='darkred')
lines(c(mean(m.spn.ns)-sd(m.spn.ns),mean(m.spn.ns)+sd(m.spn.ns))~rep(2.05,2),col='darkred')
lines(c(mean(sc.catch.mix)-sd(sc.catch.mix),mean(sc.catch.mix)+sd(sc.catch.mix))~rep(1,2),col='goldenrod')
lines(c(mean(m.spn.mix)-sd(m.spn.mix),mean(m.spn.mix)+sd(m.spn.mix))~rep(2,2),col='goldenrod')
lines(c(mean(sc.catch.st),mean(m.spn.st))~c(0.95,1.95),lwd=2,col='navy')
points(c(mean(sc.catch.st),mean(m.spn.st))~c(0.95,1.95),pch=21,bg='navy',cex=2)
lines(c(mean(sc.catch.ns),mean(m.spn.ns))~c(1.05,2.05),lwd=2,col='darkred')
points(c(mean(sc.catch.ns),mean(m.spn.ns))~c(1.05,2.05),pch=21,bg='darkred',cex=2)
lines(c(median(sc.catch.mix),median(m.spn.mix))~c(1,2),lwd=2,col='goldenrod')
points(c(median(sc.catch.mix),median(m.spn.mix))~c(1,2),pch=21,bg='goldenrod',cex=2)
mtext(side=1,'Scaled Catch',at=1,line=1)
mtext(side=1,'Prop. years above EG',at=2,line=1)
text(x=par('usr')[1]*1.5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,'Stable',col='navy')
text(x=par('usr')[1]*1.5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.2,'Non-stationary',col='darkred')
text(x=par('usr')[1]*1.5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.25,'Mixed',col='goldenrod')

tp=c(1,round(nrow(sims.st[[1]])/2),nrow(sims.st[[1]])) #time slice at t=1, half way, final

eg.limit.st=eg.scalar*sims.st[[1]]$smsy.est #escapement goals through time for stable simulation
eg.limit.ns=eg.scalar*sims.ns[[1]]$smsy.est #escapement goals through time for non-stationary simulation
eg.limit.mix=eg.scalar*sims.mix[[1]]$smsy.est #escapement goals through time for non-stationary simulation

upper.tar.st=upper.tar.scalar*sims.st[[1]]$smsy.est #upper limit beyond which harvest is static at Umsy
upper.tar.ns=upper.tar.scalar*sims.ns[[1]]$smsy.est #upper limit beyond which harvest is static at Umsy
upper.tar.mix=upper.tar.scalar*sims.mix[[1]]$smsy.est #upper limit beyond which harvest is static at Umsy

target.ER.st=U.scalar*sims.st[[1]]$umsy.est #max. harvest rate
target.ER.ns=U.scalar*sims.ns[[1]]$umsy.est #max. harvest rate
target.ER.mix=U.scalar*sims.mix[[1]]$umsy.est #max. harvest rate

par(mfrow=c(2,2))
hcr_plot_est(sim.df=sims.st[[1]],tp=tp,eg.limit=eg.limit.st,upper.tar=upper.tar.st,target.ER=target.ER.st,U.min=U.min,title='Stable Ricker model')
hcr_plot_est(sim.df=sims.mix[[1]],tp=tp,eg.limit=eg.limit.mix,upper.tar=upper.tar.mix,target.ER=target.ER.mix,U.min=U.min,title='Mixed models')
hcr_plot_est(sim.df=sims.ns[[1]],tp=tp,eg.limit=eg.limit.ns,upper.tar=upper.tar.ns,target.ER=target.ER.ns,U.min=U.min,title='Non-stationary models')
