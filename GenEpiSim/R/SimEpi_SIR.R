#' @title Stochastc SIR model (Osvaldo A.) 
#' @description Generating data from an stochastc SIR model
#' @param pop the population data
#' @param nReps number of replicates
#' @param nGens number of generations
#' @param repl current replicate
#' @param gen current generation
#' @param r.g selection accuracy for susceptibility
#' @param r.f selection accuracy for infectivity
#' @param TRbeta average transmission coefficient beta
#' @param RRgamma recovery rate gamma
#' @param SigG.g genetic variance for susceptibility
#' @param SigG.f genetic variance for infectivity
#' @param SigE.g environmental variance for susceptibility
#' @param SigE.f environmental variance for infectivity
#' @param sires number of sires
#' @param dpsire number of dams per sire
#' @param N population size
#' @param gr.size group size
#' @param ngroups number of groups
#' @param freqDep T/F for frequency or density dependent force of infection
#' @param reprnum T/F modify calculation of infection rate to calculate reproduction number R0
#' @param SIRdata SIR data
#' @return list with SIR data frame as time series (SIRts) and on a population level (pop)
#' @importFrom stats rnorm
#' @importFrom stats rexp
#
SIRdataGen<-function(pop, nReps, nGens, repl, gen, r.g, r.f, TRbeta, RRgamma, 
                     SigG.g, SigG.f, SigE.g, SigE.f, sires, dpsire, N, gr.size, ngroups, freqDep, reprnum, SIRdata) 
{
  
  if (freqDep) freqDepTerm <-gr.size else  freqDepTerm<-1  
  
  
  #--------Assigning index cases and groups (Random family assignment)----#
  
  #index cases
  pop[,"s"] <- sample(c(rep(0,ngroups),rep(1,N-ngroups)),replace=F)       #################### s==0 is the index case index
  
  #groups:
  pop$group <- NA
  pop[pop$s==0,]$group <- 1:ngroups
  pop[pop$s==1,]$group <- sample(rep(1:ngroups,gr.size-1),replace=F)
  pop[pop$s==0,"Tinf"] <- 0
  
  
  #--------------------Generating epidemics-------------------------#
  
  #Infected indicator
  pop$I <- 1-pop$s             #indicator if the individual was infected or not, where s the the index cases (initialisation - at the individual level)
  
  pop[pop$s==1,"Tinf"] <--1    #tau the "time-to-infection"
  pop[pop$s==0,"Tinf"] <-0     #index case: infected at the start
  pop$last.Tinf<-0
  
  #time to recovery and its indicator  
  pop$R<-0  ; pop$Trec<-(-999)
  
  
  #---------------------Transmission rate vectors--------------------#
  eRates<-list()
  for (k in 1:ngroups) eRates[[k]]<-rep(NA,nr=gr.size)   #the Rate of the next event(based on beta or gamma, depending on which compartment the individual is in)
  
  SIRdata<-list()
  SIRts<-list() 
  
  for (i in 1:ngroups)
  {
    
    ##   repeat{
    pop.gr<-pop[pop$group==i,]
    
    Group<-i
    SIRts[[i]]<-data.frame(repl, gen, Group, time=0,S=gr.size-1,I=1,R=0)
    
    nextT<-0 
    ll<-1
    # epidemic ends when there is no S to infect or no I to recover
    
    while (sum(pop.gr$I) > 0 & sum(pop.gr$I==0) > 0)              #end: infected to be recovered==0 or susceptible to be infected==0
    { # event rates
      for (j in 1:gr.size){
        #(recovered individuals cannot be selected, if recovered R==1 and then it can't be selected eRates==0)
        if (pop.gr$R[j]==1) eRates[[i]][j]=0   
        else{
          if (pop.gr$I[j]==0){             #if the individual is not infected, rate of infection for each susceptible
            if (reprnum) eRates[[i]][j]=exp(pop.gr$g[j])*(TRbeta/freqDepTerm)*exp(pop.gr[pop.gr$s==0,]$f) else eRates[[i]][j]=exp(pop.gr$g[j])*(TRbeta/freqDepTerm)*sum(exp(pop.gr[pop.gr$I==1,]$f)) 
          }
          
          else {                        # constant recovery rate of each infected (no variation)
            eRates[[i]][j]=RRgamma }}  
      }
      
      nextT<-nextT+rexp(1, rate = sum(eRates[[i]]))
      IDnextT<-sample(pop.gr$ID,size=1,prob=eRates[[i]]/sum(eRates[[i]]))
      
      
      SIRts[[i]]<-rbind(SIRts[[i]],SIRts[[i]][nrow(SIRts[[i]]),])       #data with number of SIR over time
      jj<-nrow(SIRts[[i]])
      
      #ll > 1                                         #to avoid that the first event if recovered of index case
      
      if (pop.gr[pop.gr$ID==IDnextT,]$I==1){          #recovery of an infected
        pop.gr[pop.gr$ID==IDnextT,"I"] <- 0
        pop.gr[pop.gr$ID==IDnextT,"R"] <- 1
        pop.gr[pop.gr$ID==IDnextT,"Trec"] <- nextT
        
        SIRts[[i]][jj,]$time<-nextT
        SIRts[[i]][jj,]$I<-SIRts[[i]][jj,]$I-1
        SIRts[[i]][jj,]$R<-SIRts[[i]][jj,]$R+1 
      }
      
      else{                                           #infection of a susceptible         
        pop.gr[pop.gr$ID==IDnextT,"I"] <- 1
        pop.gr[pop.gr$ID==IDnextT,"Tinf"] <- nextT 
        
        SIRts[[i]][jj,]$time<-nextT
        SIRts[[i]][jj,]$S<-SIRts[[i]][jj,]$S-1
        SIRts[[i]][jj,]$I<-SIRts[[i]][jj,]$I+1
      } 
      ll <- ll +1
    }
    
    ##        if (sum(pop.gr$R) >= 2) break
    ##            }
    
    row.names(SIRts[[i]])<-1:nrow(SIRts[[i]]) 
    
    pop.gr[pop.gr$Tinf==max(pop.gr$Tinf),"last.Tinf"]<-1
    pop.gr[pop.gr$Tinf==-1,"Tinf"]<-(-999)
    
    pop[pop$group==i,]<-pop.gr
  }
  
  
  SIRdata<-list(pop=pop, SIRts=SIRts)
  
  return(SIRdata)
}




