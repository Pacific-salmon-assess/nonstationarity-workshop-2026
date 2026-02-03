#functions
salmon_sim=function(log.a,smax,sigma,N,form=c('static','autocorr'),rho=NA){
  if(form=='static'){
  #age structure
  #this accounts for different maturations as is typical for most salmon species
  A=5 #max age class
  age=c(3,4,5) #ages at return - set to 3,4,5
  age.p=c(0.25,0.5,0.25) # proportions returning at each age class = note if modifying this, it must sum to 1. This is fixed in this example to keep things simple but one could vary with some noise from year-to-year using e.g. dirichlet distribution
  
  umsy=samEst::umsyCalc(log.a)
  K = log.a*smax #carrying capacity
  #simulation parameters
  N=50 #number of simulated years for our spawner-recruit curve
  L=N+A*2 #total simulation length, add 2x max age for starting cohorts to seed the simulation and for final incomplete brood years - these will be dropped later
 
  U=runif(L,umsy*0.25,umsy*1.25) #harvest rate - uniformly drawn with constraints around 25% and 1.25x Umsy
  
  #simulation vectors
  logRS=numeric(L) #productivity, log(R/S), for each brood cohort
  S=numeric(L);S[1:c(A)]=runif(A,K*0.1,K*1.25) #initial spawners, randomized between 0.5-1.5*K 
  R=numeric(L) #recruits in each brood cohort year
  Rs=numeric(L) #run size, ie. the returning abundance in each year (recruits x age structure)
  
  #first recruitments and run sizes from randomized spawners to initialize simulation
  logRS[1:c(A)]=rnorm(A, mean=(log.a - S[1:c(A)]/smax),sd=sigma)
  R[1:c(A)]=exp(logRS[1:c(A)])*S[1:c(A)] #transform to recruits by converting log(R/S) to R/S times spawners
  for(t in 1:c(A)){ #initialize first 2 generations of run sizes
    for(a in 1:length(age)){
      Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a] #add proportion of recruits by age class
    }
  }
  #this has setup the the initial run sizes for the simulation
  
  for(t in c(A+1):c(L-A)){ #from years max age + 1
    #for each year, we will harvest the run, get the remaining spawning abundance
    S[t]=Rs[t]*(1-U[t]) # spawning abundance (escapement) is run size minus harvest rate
    
    #recruitment for this brood cohort (year t) - mean is equal to the expectation from the ricker model based on input parameters
    logRS[t]=rnorm(1, mean=(log.a - S[t]/smax),sd=sigma)
    
    #convert to number of recruits for brood cohort
    R[t]=exp(logRS[t])*S[t]
    
    #stagger recruits into future run sizes based on age maturation schedule:
    for(a in 1:length(age)){
      Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a]
    }
  }
  #We'll trim out the initialization period (first generation) and last generation (incomplete cohorts)
  S=S[c(A+1):c(L-A)]
  R<-R[c(A+1):c(L-A)]
  logRS<-logRS[c(A+1):c(L-A)]
  
  df=data.frame(S=S,R=R,logRS=logRS,eps=logRS-(log.a-S/smax),by=seq(1:length(S)))
  }
  
  if(form=='autocorr'){
    if(is.na(rho)==T){print(' you must specify rho (-1 to 1 range)')}
    if(rho< -1|rho >1){print('invalid rho parameter, must be -1 to 1')}
    #age structure
    #this accounts for different maturations as is typical for most salmon species
    A=5 #max age class
    age=c(3,4,5) #ages at return - set to 3,4,5
    age.p=c(0.25,0.5,0.25) # proportions returning at each age class = note if modifying this, it must sum to 1. This is fixed in this example to keep things simple but one could vary with some noise from year-to-year using e.g. dirichlet distribution
    
    umsy=samEst::umsyCalc(log.a)
    K = log.a*smax #carrying capacity
    
    #simulation parameters
    N=50 #number of simulated years for our spawner-recruit curve
    L=N+A*2 #total simulation length, add 2x max age for starting cohorts to seed the simulation and for final incomplete brood years - these will be dropped later
   
    U=runif(L,umsy*0.25,umsy*1.25) #harvest rate - uniformly drawn (you can change this but be careful of exceeding, 0-1 boundaries)
    
    #U=plogis(normal(L,0,0.5)) - alternative harvest rate protocol using a logit-transformed normal draw, less variable than above, can comment the above line with # and this off to switch 
    
    #simulation vectors
    logRS=numeric(L) #productivity, log(R/S), for each brood cohort
    S=numeric(L);S[1:c(A)]=runif(A,K*0.25,K*1.25) #initial spawners, randomized between 0.5-1.5*K 
    R=numeric(L) #recruits in each brood cohort year
    Rs=numeric(L) #run size, ie. the returning abundance in each year (recruits x age structure)
    eps=numeric(L) #residual productivity - epsilon
    
    #first recruitments and run sizes from randomized spawners to initialize simulation
    #first draw has no autocorrelation, all subsequent draws will include the autocorrelation component
    logRS[1]=rnorm(1, mean=(log.a - S[1]/smax),sd=sigma)
    eps[1]=logRS[1]-(log.a - S[1]/smax) #first residual, realized log(R/S) minus expected (mean) log(R/S)
    for(t in 2:A){ #all subsequent have expectation shifted by rho*eps[t-1]
      logRS[t]=rnorm(1, mean=(log.a - S[t]/smax + rho*eps[t-1]),sd=sigma)
      eps[t]=logRS[t]-(log.a - S[t]/smax)
    }
    R[1:c(A)]=exp(logRS[1:c(A)])*S[1:c(A)] #transform to recruits by converting log(R/S) to R/S times spawners
    for(t in 1:c(A)){ #initialize first 2 generations of run sizes
      for(a in 1:length(age)){
        Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a] #add proportion of recruits by age class
      }
    }
    #this has setup the the initial run sizes for the simulation
    
    for(t in c(A+1):c(L-A)){ #from years max age + 1
      #for each year, we will harvest the run, get the remaining spawning abundance
      S[t]=Rs[t]*(1-U[t]) # spawning abundance (escapement) is run size minus harvest rate
      
      #recruitment for this brood cohort (year t) - mean is equal to the expectation from the ricker model based on input parameters
      logRS[t]=rnorm(1, mean=(log.a - S[t]/smax+ rho*eps[t-1]),sd=sigma)
      eps[t]=logRS[t]-(log.a - S[t]/smax)
      #convert to number of recruits for brood cohort
      R[t]=exp(logRS[t])*S[t]
      
      #stagger recruits into future run sizes based on age maturation schedule:
      for(a in 1:length(age)){
        Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a]
      }
    }
    #We'll trim out the initialization period (first generation) and last generation (incomplete cohorts)
    S=S[c(A+1):c(L-A)]
    R<-R[c(A+1):c(L-A)]
    logRS<-logRS[c(A+1):c(L-A)]
    eps<- eps[c(A+1):c(L-A)]
    
    df=data.frame(S=S,R=R,logRS=logRS,eps=eps,by=seq(1:length(S)))
  }
  
  return(df)
}

salmon_sim.tv=function(log.a0,smax0,sigma,N,tv.par=c('a','b','both'),p.change,p.change2,tv.form=c('linear','rw','regime'),reg.length=NULL){
  
    #age structure
    #this accounts for different maturations as is typical for most salmon species
    A=5 #max age class
    age=c(3,4,5) #ages at return - set to 3,4,5
    age.p=c(0.25,0.5,0.25) # proportions returning at each age class = note if modifying this, it must sum to 1. This is fixed in this example to keep things simple but one could vary with some noise from year-to-year using e.g. dirichlet distribution
    
    umsy=samEst::umsyCalc(log.a0)
    K=log.a0*smax0
    #simulation parameters
    N=50 #number of simulated years for our spawner-recruit curve
    L=N+A*2 #total simulation length, add 2x max age for starting cohorts to seed the simulation and for final incomplete brood years - these will be dropped later
    
    U=runif(L,umsy*0.25,umsy*1.25) #harvest rate - uniformly drawn with constraints around 25% and 1.25x Umsy
    
    if(tv.par=='a'){
      smax=rep(smax0,L)
      if(tv.form=='linear'){
        log.a=numeric(L)
        log.a[c(1:A)]=log.a0
        log.a[c(A+1):c(L-A)]=seq(log.a0,log.a0*(1+p.change),length.out=N)
        log.a[c(L-A+1):L]=log.a0*(1+p.change)
      }
      if(tv.form=='rw'){
        log.a=numeric(L)
        log.a[c(1:A)]=log.a0
        p.trend=p.change/N
        for(t in 1:c(N+A)){
          log.a[A+t]=log.a[A+t-1]*(1+rnorm(1,p.trend,0.05+abs(p.trend)*2))
        }
      }
      if(tv.form=='regime'){
        if(is.null(reg.length)==T){print('must specify regime length')}
        log.a=numeric(L)
        log.a=rep(rep(c(log.a0,log.a0*(1+p.change)),each=reg.length),length.out=L)
      }
    } 
    if(tv.par=='b'){
      log.a=rep(log.a0,L)
      if(tv.form=='linear'){
        smax=numeric(L)
        smax[c(1:A)]=smax0
        smax[c(A+1):c(L-A)]=seq(smax0,smax0*(1+p.change),length.out=N)
        smax[c(L-A+1):L]=smax0*(1+p.change)
      }
      if(tv.form=='rw'){
        smax=numeric(L)
        smax[c(1:A)]=smax0
        p.trend=p.change/N
        for(t in 1:c(N+A)){
          smax[A+t]=smax[A+t-1]*(1+rnorm(1,p.trend,0.05+abs(p.trend)*2))
        }
      }
      if(tv.form=='regime'){
        if(is.null(reg.length)==T){print('must specify regime length')}
        smax=numeric(L)
        smax=rep(rep(c(smax0,smax0*(1+p.change)),each=reg.length),length.out=L)
      }
    } 
    if(tv.par=='both'){
      if(is.null(reg.length)==T){print('must specify regime length')}
      smax=numeric(L)
      log.a=numeric(L)
      log.a=rep(rep(c(log.a0,log.a0*(1+p.change)),each=reg.length),length.out=L)
      smax=rep(rep(c(smax0,smax0*(1+p.change2)),each=reg.length),length.out=L)
      }
    #simulation vectors
    logRS=numeric(L) #productivity, log(R/S), for each brood cohort
    S=numeric(L);S[1:c(A)]=runif(A,K*0.1,K*1.25) #initial spawners, randomized between 0.5-1.5*K 
    R=numeric(L) #recruits in each brood cohort year
    Rs=numeric(L) #run size, ie. the returning abundance in each year (recruits x age structure)
    
    #first recruitments and run sizes from randomized spawners to initialize simulation
    logRS[1:c(A)]=rnorm(A, mean=(log.a[1:c(A)] - S[1:c(A)]/smax[1:c(A)]),sd=sigma)
    R[1:c(A)]=exp(logRS[1:c(A)])*S[1:c(A)] #transform to recruits by converting log(R/S) to R/S times spawners
    for(t in 1:c(A)){ #initialize first 2 generations of run sizes
      for(a in 1:length(age)){
        Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a] #add proportion of recruits by age class
      }
    }
    #this has setup the the initial run sizes for the simulation
    
    for(t in c(A+1):c(L-A)){ #from years max age + 1
      #for each year, we will harvest the run, get the remaining spawning abundance
      S[t]=Rs[t]*(1-U[t]) # spawning abundance (escapement) is run size minus harvest rate
      
      #recruitment for this brood cohort (year t) - mean is equal to the expectation from the ricker model based on input parameters
      logRS[t]=rnorm(1, mean=(log.a[t] - S[t]/smax[t]),sd=sigma)
      
      #convert to number of recruits for brood cohort
      R[t]=exp(logRS[t])*S[t]
      
      #stagger recruits into future run sizes based on age maturation schedule:
      for(a in 1:length(age)){
        Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a]
      }
    }
    #We'll trim out the initialization period (first generation) and last generation (incomplete cohorts)
    S=S[c(A+1):c(L-A)]
    R<-R[c(A+1):c(L-A)]
    logRS<-logRS[c(A+1):c(L-A)]
    smax<- smax[c(A+1):c(L-A)]
    log.a<- log.a[c(A+1):c(L-A)]
    
    df=data.frame(S=S,R=R,logRS=logRS,eps=logRS-(log.a-S/smax),by=seq(1:length(S)),loga.t=log.a,smax.t=smax)
 
  return(df)
}

salmon_sim.tv_hcr=function(log.a0,smax0,sigma,N,tv.par=c('a','b','both'),p.change,p.change2,tv.form=c('linear','rw','regime'),reg.length=NULL,hcr.form=c('stable','rw','hmm'),hcr.par=c('a','b','both'),assess.freq=10,eg.scalar=1,upper.tar.scalar=1.5,U.scalar=1,U.min=0.025){
  
  #age structure
  #this accounts for different maturations as is typical for most salmon species
  A=5 #max age class
  age=c(3,4,5) #ages at return - set to 3,4,5
  age.p=c(0.25,0.5,0.25) # proportions returning at each age class = note if modifying this, it must sum to 1. This is fixed in this example to keep things simple but one could vary with some noise from year-to-year using e.g. dirichlet distribution
  
  K=log.a0*smax0
  #simulation parameters
  N=50 #number of simulated years for our spawner-recruit curve
  L=N+A*5 #total simulation length, add 4x max age for starting cohorts to seed the simulation and for final incomplete brood years - these will be dropped later
  
  
  if(tv.par=='a'){
    smax=rep(smax0,L)
    if(tv.form=='linear'){
      log.a=numeric(L)
      log.a[1:c(A*4)]=log.a0
      log.a[c(A*4+1):c(L-A)]=seq(log.a0,log.a0*(1+p.change),length.out=N)
      log.a[c(L-A+1):L]=log.a0*(1+p.change)
    }
    if(tv.form=='rw'){
      log.a=numeric(L)
      log.a[1:c(A*4)]=log.a0
      p.trend=p.change/N
      for(t in 1:c(N)){
        log.a[A*4+t]=log.a[A*4+t-1]*(1+rnorm(1,p.trend,0.05+abs(p.trend)*2))
      }
    }
    if(tv.form=='regime'){
      if(is.null(reg.length)==T){print('must specify regime length')}
      log.a=numeric(L)
      log.a[1:c(A*4)]=log.a0
      log.a[c(A*4+1):L]=rep(rep(c(log.a0,log.a0*(1+p.change)),each=reg.length),length.out=c(L-A*4))
    }
  } 
  if(tv.par=='b'){
    log.a=rep(log.a0,L)
    if(tv.form=='linear'){
      smax=numeric(L)
      smax[1:c(A*4)]=smax0
      smax[c(A*4+1):c(L-A)]=seq(smax0,smax0*(1+p.change),length.out=N-A)
      smax[c(L-A+1):L]=smax0*(1+p.change)
    }
    if(tv.form=='rw'){
      smax=numeric(L)
      smax[1:c(A*4)]=smax0
      p.trend=p.change/N
      for(t in 1:c(N)){
        smax[A*4+t]=smax[A*4+t-1]*(1+rnorm(1,p.trend,0.05+abs(p.trend)*2))
      }
    }
    if(tv.form=='regime'){
      if(is.null(reg.length)==T){print('must specify regime length')}
      smax=numeric(L)
      smax[1:c(A*4)]=smax0
      smax=rep(rep(c(smax0,smax0*(1+p.change)),each=reg.length),length.out=c(L-A*4))
    }
  } 
  if(tv.par=='both'){
    if(is.null(reg.length)==T){print('must specify regime length')}
    smax=numeric(L)
    log.a=numeric(L)
    log.a[1:c(A*4)]=log.a0
    smax[1:c(A*4)]=smax0
    log.a=rep(rep(c(log.a0,log.a0*(1+p.change)),each=reg.length),length.out=L)
    smax=rep(rep(c(smax0,smax0*(1+p.change2)),each=reg.length),length.out=L)
  }
  #simulation vectors
  logRS=numeric(L) #productivity, log(R/S), for each brood cohort
  S=numeric(L);S[1:c(A)]=runif(A,K*0.1,K*1.25) #initial spawners, randomized between 0.5-1.5*K 
  R=numeric(L) #recruits in each brood cohort year
  Rs=numeric(L) #run size, ie. the returning abundance in each year (recruits x age structure)
  U=numeric(L) #exploitation rate
  C=numeric(L) #number of fish harvested through time
  Smsy.true=samEst::smsyCalc(log.a,1/smax) #true spawners at MSY
  Umsy.true=samEst::umsyCalc(log.a) #true U at MSY
  Smsy.est=numeric(L) #estimated Smsy through time
  Umsy.est=numeric(L) #estimated Umsy through time
  
  #first recruitments and run sizes from randomized spawners to initialize simulation
  logRS[1:c(A)]=rnorm(A, mean=(log.a[1:c(A)] - S[1:c(A)]/smax[1:c(A)]),sd=sigma)
  R[1:c(A)]=exp(logRS[1:c(A)])*S[1:c(A)] #transform to recruits by converting log(R/S) to R/S times spawners
  for(t in 1:c(A)){ #initialize first 2 generations of run sizes
    for(a in 1:length(age)){
      Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a] #add proportion of recruits by age class
    }
  }
  #this has setup the the initial run sizes for the simulation
  
  #run it for the first 3 generations to generate initial S-R data for assessments - random exploitation based on initial Umsy
  for(t in c(A+1):c(A*4)){ #from years max age + 1
    
    U[t]=plogis(rnorm(1,qlogis(Umsy.true[t]),0.6))
      
    #for each year, we will harvest the run, get the remaining spawning abundance
    S[t]=Rs[t]*(1-U[t]) # spawning abundance (escapement) is run size minus harvest rate
    C[t]=Rs[t]*U[t] #catch in each year
    
    #recruitment for this brood cohort (year t) - mean is equal to the expectation from the ricker model based on input parameters
    logRS[t]=rnorm(1, mean=(log.a[t] - S[t]/smax[t]),sd=sigma)
    
    #convert to number of recruits for brood cohort
    R[t]=exp(logRS[t])*S[t]
    
    #stagger recruits into future run sizes based on age maturation schedule:
    for(a in 1:length(age)){
      Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a]
    }
  }
  
  if(hcr.form=='stable'){
    df=data.frame(S=S[1:c(A*4)],by=seq(1,c(A*4)),logRS=logRS[1:c(A*4)])
    est=samEst::ricker_TMB(data=df,silent=T,ac=T)
    Smsy.est[1:c(A*4)]=est$Smsy
    Umsy.est[1:c(A*4)]=est$Umsy 
  }
  if(hcr.form=='rw'){
    df=data.frame(S=S[1:c(A*4)],by=seq(1,c(A*4)),logRS=logRS[1:c(A*4)])
    est=samEst::ricker_rw_TMB(data=df,tv.par=hcr.par,silent=T)
    Smsy.est[1:c(A*4)]=mean(est$Smsy[c(A*4-assess.freq):c(A*4)])
    Umsy.est[1:c(A*4)]=mean(est$Umsy[c(A*4-assess.freq):c(A*4)]) 
  }
  if(hcr.form=='hmm'){
    df=data.frame(S=S[1:c(A*4)],by=seq(1,c(A*4)),logRS=logRS[1:c(A*4)])
    est=samEst::ricker_hmm_TMB(data=df,tv.par=hcr.par,silent=T)
    Smsy.est[1:c(A*4)]=est$Smsy[est$regime[c(A*4)]]
    Umsy.est[1:c(A*4)]=est$Umsy[est$regime[c(A*4)]]
  }
  
  for(t in c(A*4+1):L){ #from years max age + 1
    if(t %% assess.freq == 0){
      if(hcr.form=='stable'){
        df=data.frame(S=S[1:t-1],by=seq(1,t-1),logRS=logRS[1:t-1])
        est=samEst::ricker_TMB(data=df,silent=T,ac=T)
        Smsy.est[t]=est$Smsy
        Umsy.est[t]=est$Umsy 
      }else if(hcr.form=='rw'){
        df=data.frame(S=S[1:t-1],by=seq(1,t-1),logRS=logRS[1:t-1])
        est=samEst::ricker_rw_TMB(data=df,tv.par=hcr.par,silent=T)
        Smsy.est[t]=mean(est$Smsy[c(t-assess.freq):c(t-1)])
        Umsy.est[t]=mean(est$Umsy[c(t-assess.freq):c(t-1)]) 
      }else   if(hcr.form=='hmm'){
        df=data.frame(S=S[1:t-1],by=seq(1,t-1),logRS=logRS[1:t-1])
        est=samEst::ricker_hmm_TMB(data=df,tv.par=hcr.par,silent=T)
        Smsy.est[t]=mean(est$Smsy[est$regime[c(t-assess.freq):c(t-1)]])
        Umsy.est[t]=mean(est$Umsy[est$regime[c(t-assess.freq):c(t-1)]])
      }
    }else{
      Umsy.est[t]=Umsy.est[t-1]
      Smsy.est[t]=Smsy.est[t-1]
    }
    
    if(Rs[t]<Smsy.est[t]*eg.scalar){
      U[t]=plogis(rnorm(1,qlogis(U.min),0.4))
    }else if(Rs[t]>Smsy.est[t]*eg.scalar&Rs[t]<Smsy.est[t]*upper.tar.scalar){
      U.it=U.min+Umsy.est[t]*U.scalar*((Rs[t]-Smsy.est[t]*eg.scalar)/(Smsy.est[t]*upper.tar.scalar-Smsy.est[t]*eg.scalar))
      U.it=min(U.it,Umsy.est[t]*U.scalar)
      U[t]=plogis(rnorm(1,qlogis(U.it),0.4))
    }else if(Rs[t]>Smsy.est[t]*upper.tar.scalar){
      U[t]=plogis(rnorm(1,qlogis(Umsy.est[t]*U.scalar),0.4))
    }
    
    #for each year, we will harvest the run, get the remaining spawning abundance
    S[t]=Rs[t]*(1-U[t]) # spawning abundance (escapement) is run size minus harvest rate
    C[t]=Rs[t]*U[t] #catch in each year
    
    #recruitment for this brood cohort (year t) - mean is equal to the expectation from the ricker model based on input parameters
    logRS[t]=rnorm(1, mean=(log.a[t] - S[t]/smax[t]),sd=sigma)
    
    #convert to number of recruits for brood cohort
    R[t]=exp(logRS[t])*S[t]
    
    #stagger recruits into future run sizes based on age maturation schedule:
    for(a in 1:length(age)){
      Rs[t+age[a]]=Rs[t+age[a]]+R[t]*age.p[a]
    }
  }
  #We'll trim out the initialization period (first generations and initial fishery dynamics without HCR) and last generation (incomplete cohorts)
  S=S[c(A*4+1):c(L-A)]
  R<-R[c(A*4+1):c(L-A)]
  logRS<-logRS[c(A*4+1):c(L-A)]
  smax<- smax[c(A*4+1):c(L-A)]
  log.a<- log.a[c(A*4+1):c(L-A)]
  U<- U[c(A*4+1):c(L-A)]
  C<- C[c(A*4+1):c(L-A)]
  Rs<- Rs[c(A*4+1):c(L-A)]
  Smsy.true<- Smsy.true[c(A*4+1):c(L-A)]
  Umsy.true<- Umsy.true[c(A*4+1):c(L-A)]
  Smsy.est<- Smsy.est[c(A*4+1):c(L-A)]
  Umsy.est<- Umsy.est[c(A*4+1):c(L-A)]
  
  
  df=data.frame(S=S,R=R,logRS=logRS,eps=logRS-(log.a-S/smax),by=seq(1:length(S)),loga.t=log.a,smax.t=smax,run.size=Rs,catch=C,harvest.rate=U,smsy.true=Smsy.true,smsy.est=Smsy.est,umsy.true=Umsy.true,umsy.est=Umsy.est)
  
  return(df)
}