#' @title Function to calculate Response to selection for susceptibility and infectivity over generations  (Smaragda Tsairidou 7.9.16)  
#' @description Calculates response to selection by calculating the change in the mean and variance of susceptibility and infectivity over generations 
#' @param pop population data
#' @param nReps number of replicates
#' @param nGens number of generations
#' @param rr current replicate
#' @param kk current generation
#' @param r.g selection accuracy for susceptibility
#' @param r.f selection accuracy for infectivity
#' @param SigG.g genetic variance for susceptibility
#' @param SigG.f genetic variance for infectivity
#' @param SigE.g environmental variance for susceptibility
#' @param SigE.f environmental variance for infectivity
#' @param i.g selection intensity for susceptibility
#' @param i.f selection intensity for infectivity
#' @param km.g the Bulmer effect factor for susceptibility
#' @param km.f the Bulmer effect factor for infectivity
#' @param sires number of sires
#' @param dpsire number of dams per sire
#' @param N population size
#' @param gr.size group size
#' @param ngroups number of groups
#' @param Data_SimSel data from selection
#' @param Results_SimSel results from selection
#' @return Dataframe with new population parameters
#' @importFrom utils write.table
#' @importFrom utils read.table
SimSel<-function(pop, nReps, nGens, rr, kk, r.g, r.f, SigG.g, SigG.f, SigE.g, SigE.f, i.g, i.f, km.g, km.f, 
                 sires, dpsire, N, gr.size, ngroups, Data_SimSel, Results_SimSel)
{
  
  rhoG = rhoE = 0
  
  #---------------Generating BVs and Phenotypes for the base generation----------------# 
  
  if (kk==1) {
    
    print(paste("Generating generation ", kk, sep="")) 
    
    #parents
    covG <- rhoG*sqrt(SigG.g)*sqrt(SigG.f)                     
    sigmaG <- matrix(c(SigG.g,covG,covG,SigG.f),ncol=2)          
    auxBV.s <- mvrnorm(n = sires, mu=c(0,0), Sigma=sigmaG)
    auxBV.d <- mvrnorm(n = N, mu=c(0,0), Sigma=sigmaG)
    BV.s <- data.frame(ID=1:sires,ng.off=NA,Ag=auxBV.s[,1],Af=auxBV.s[,2])
    BV.d <- data.frame(ID=1:N,sire=sort(rep(1:sires,dpsire)),Ag=auxBV.d[,1],Af=auxBV.d[,2])
    
    #offspring                        
    for (i in 1:sires){
      auxBVmend <- mvrnorm(n = dpsire, mu=c(0,0), Sigma=0.5*sigmaG)                        
      pop[pop$sire==i,]$Ag <- 0.5*BV.s$Ag[i]+0.5*BV.d[BV.d$sire==i,]$Ag+auxBVmend[,1]     
      pop[pop$sire==i,]$Af <- 0.5*BV.s$Af[i]+0.5*BV.d[BV.d$sire==i,]$Af+auxBVmend[,2]
    } 
    SigG.g_off_new<- -999
    SigG.f_off_new<- -999
    mu.g <- 0
    mu.f <- 0 
    R.g <- 0
    R.f <- 0
    
    covE <- rhoE*sqrt(SigE.g)*sqrt(SigE.f)
    sigmaE <- matrix(c(SigE.g,covE,covE,SigE.f),ncol=2)
    Eaux <- mvrnorm(n = N, mu=c(0,0), Sigma=sigmaE)             
    pop$eg <- Eaux[,1]; pop$ef<-Eaux[,2]
    pop$g <- mu.g+pop$Ag+pop$eg                                    
    pop$f <- mu.f+pop$Af+pop$ef          
    
    susc<-mean(pop$g)
    inf<-mean(pop$f)
    
  } else if (kk==2) { 
    
    print(paste("Generating generation ", kk, sep=""))
    
    R.g = i.g * r.g * sqrt(SigG.g)             
    R.f = i.f * r.f * sqrt(SigG.f)                
    
    mu.g<-Results_SimSel$mu.g + (-R.g)
    mu.f<-Results_SimSel$mu.f + (-R.f)
    
    SigG.g_sire_new <- ((1 - ((r.g^2)*km.g)) * 0.25*SigG.g)
    SigG.f_sire_new <- ((1 - ((r.f^2)*km.f)) * 0.25*SigG.f)
    SigG.g_off_new <- (0.25*SigG.g_sire_new + 0.25*SigG.g + 0.5*SigG.g)
    SigG.f_off_new <- (0.25*SigG.f_sire_new + 0.25*SigG.f + 0.5*SigG.f)
    
    pop$Ag <- rnorm(n = N, mean=0, sd=sqrt(SigG.g_off_new))   
    pop$Af <- rnorm(n = N, mean=0, sd=sqrt(SigG.f_off_new)) 
    
    pop$eg <- rnorm(n = N, mean = 0, sd = sqrt(SigE.g))
    pop$ef <- rnorm(n = N, mean = 0, sd = sqrt(SigE.f))
    pop$g <- mu.g+pop$Ag+pop$eg                                         
    pop$f <- mu.f+pop$Af+pop$ef
    
    susc<-mean(pop$g)
    inf<-mean(pop$f)
    
  } else {
    
    
    print(paste("Generating generation ",kk, sep=""))
    
    R.g = i.g * r.g * sqrt(Results_SimSel$SigG.g_off_new)             
    R.f = i.f * r.f * sqrt(Results_SimSel$SigG.f_off_new)             
    
    mu.g<-Results_SimSel$mu.g + (-R.g)
    mu.f<-Results_SimSel$mu.f + (-R.f)   
    
    SigG.g_sire_new <- ((1 - ((r.g^2)*km.g)) * 0.25*Results_SimSel$SigG.g_off_new)
    SigG.f_sire_new <- ((1 - ((r.f^2)*km.f)) * 0.25*Results_SimSel$SigG.f_off_new)
    
    SigG.g_off_new <-(0.25*SigG.g_sire_new + 0.25*Results_SimSel$SigG.g_off_new + 0.5*SigG.g)
    SigG.f_off_new <-(0.25*SigG.f_sire_new + 0.25*Results_SimSel$SigG.f_off_new + 0.5*SigG.f) 
    
    pop$Ag <- rnorm(n = N, mean=0, sd=sqrt(SigG.g_off_new))   
    pop$Af <- rnorm(n = N, mean=0, sd=sqrt(SigG.f_off_new)) 
    
    pop$eg <- rnorm(n = N, mean = 0, sd = sqrt(SigE.g))
    pop$ef <- rnorm(n = N, mean = 0, sd = sqrt(SigE.f))
    
    pop$g <- mu.g+pop$Ag+pop$eg                                         
    pop$f <- mu.f+pop$Af+pop$ef  
    
    susc<-mean(pop$g)
    inf<-mean(pop$f)
    
  }
  
  Data_SimSel<-list(repl=rr, gen=kk, mu.g=mu.g, mu.f=mu.f, Rsusc=R.g, Rinf=R.f, mean_susc=susc, mean_inf=inf, 
                    SigG.g_off_new=SigG.g_off_new, SigG.f_off_new=SigG.f_off_new, pop=pop)
  
  return(Data_SimSel)
  
}





