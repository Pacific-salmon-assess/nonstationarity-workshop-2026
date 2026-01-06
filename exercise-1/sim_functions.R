#functions
salmon_sim=function(log.a,smax,sigma,N,form,rho=NA){
  if(form=='static'){
  #age structure
  #this accounts for different maturations as is typical for most salmon species
  A=5 #max age class
  age=c(3,4,5) #ages at return - set to 3,4,5
  age.p=c(0.25,0.5,0.25) # proportions returning at each age class = note if modifying this, it must sum to 1. This is fixed in this example to keep things simple but one could vary with some noise from year-to-year using e.g. dirichlet distribution
  
  umsy=samEst::umsyCalc(log.a)
  
  #simulation parameters
  N=50 #number of simulated years for our spawner-recruit curve
  L=N+A*2 #total simulation length, add 2x max age for starting cohorts to seed the simulation and for final incomplete brood years - these will be dropped later
 
  U=runif(L,umsy*0.25,umsy*1.25) #harvest rate - uniformly drawn (you can change this but be careful of exceeding, 0-1 boundaries)
 
  #U=plogis(normal(L,0,0.5)) - alternative harvest rate protocol using a logit-transformed normal draw, less variable than above, can comment the above line with # and this off to switch 
  
  #simulation vectors
  logRS=numeric(L) #productivity, log(R/S), for each brood cohort
  S=numeric(L);S[1:c(A)]=runif(A,K*0.25,K*1.5) #initial spawners, randomized between 0.5-1.5*K 
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
  
  df=data.frame(S=S,R=R,logRS=logRS,true.resids=logRS-(log.a-S/smax))
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
    
    #simulation parameters
    N=50 #number of simulated years for our spawner-recruit curve
    L=N+A*2 #total simulation length, add 2x max age for starting cohorts to seed the simulation and for final incomplete brood years - these will be dropped later
   
    U=runif(L,umsy*0.25,umsy*1.25) #harvest rate - uniformly drawn (you can change this but be careful of exceeding, 0-1 boundaries)
    
    #U=plogis(normal(L,0,0.5)) - alternative harvest rate protocol using a logit-transformed normal draw, less variable than above, can comment the above line with # and this off to switch 
    
    #simulation vectors
    logRS=numeric(L) #productivity, log(R/S), for each brood cohort
    S=numeric(L);S[1:c(A)]=runif(A,K*0.25,K*1.5) #initial spawners, randomized between 0.5-1.5*K 
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
    
    df=data.frame(S=S,R=R,logRS=logRS,true.resids=eps)
  }
  
  return(df)
}