#' @title Main function (developped by Smaragda Tsairidou 7.9.16) 
#' @description Calculates response to selection for susceptibility and infectivity and generates new population (function "SimSel" created by S. Tsairidou), and simulates epidemics based on an SIR model (function "SIRdataGen" created by O. Anacleto)
#' @param nReps number of replicates
#' @param nGens number of generations
#' @param r.g selection accuracy for susceptibility
#' @param r.f selection accuracy for infectivity
#' @param TRbeta average transmission coefficient beta
#' @param RRgamma recovery rate gamma
#' @param SigG.g genetic variance for susceptibility
#' @param SigG.f genetic variance for infectivity
#' @param SigE.g environmental variance for susceptibility
#' @param SigE.f environmental variance for infectivity
#' @param im.g selection intensity on the males for susceptibility
#' @param xm.g truncation point for the males for susceptibility
#' @param ifem.g selection intensity on the females for susceptibility
#' @param xfem.g truncation point for the females for susceptibility
#' @param im.f selection intensity on th emales for infectivity
#' @param xm.f truncation point on the males for infectivity
#' @param ifem.f selection intensity on the females for infectivity
#' @param xfem.f truncation point for the females for infectivity
#' @param sires number of sires
#' @param dpsire number of dams per sire
#' @param N population size
#' @param gr.size group size
#' @param ngroups number of groups
#' @param freqDep T/F for frequensy or density dependent force of infection
#' @param reprnum T/F for modified calculation of infection rate to calculate reproduction number R0
#' @param save.pop T/F to save all population data for each generation and each replicate
#' @return Generates data files: "Inf_Stat_" with population and epidemic descriptive statistics and assessment measures of epidemic severity and duration
#'                               "SIRtime_" with number of S, I and Rs at time-points 
#'                               "Pop_#repl_#gen_Susc_Inf_" if save.pop=T saving all data for each population in each generation and replicate
#' @export
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @importFrom MASS mvrnorm

#library(MASS)

#source("SimRespSel_inf_susc.R")
#source("SimEpi_SIR.R")

Sel_DataGen<-function(nReps, nGens, r.g, r.f, TRbeta, RRgamma, SigG.g, SigG.f, SigE.g, SigE.f, im.g, xm.g, ifem.g, xfem.g, im.f, xm.f, ifem.f, xfem.f,                       
                      sires, dpsire, N=sires*dpsire, gr.size=N/ngroups, ngroups=N/gr.size, freqDep, reprnum, save.pop)
{
  
  seeds <- read.table("seed_data.txt", header=F) 
  seeds <- as.numeric(seeds[,1])
  
  pfx<-paste("Susc_Inf_r07", sep="")
  file.Inf_Stat=paste("Inf_Stat_", pfx, ".txt", sep="")
  file.SIRtime=paste("SIRtime_", pfx, ".txt", sep="")
  
  #Selection intensity   
  i.g <- (0.5*im.g + 0.5*ifem.g)   
  km.g <- (im.g*(im.g-xm.g))  
  
  i.f <- (0.5*im.f + 0.5*ifem.f)   
  km.f <- (im.f*(im.f-xm.f))
  
  for (rr in 1:nReps){
    set.seed(seeds[rr]) 
    
    Results_SimSel<-list(); Results_SimSel[[1]]<-NA         
    Results_SIR<-list()
    Epi_stats<-list()
    
    for (kk in 2:(nGens+1)) {                                                                           #Generation 2 is the base generation
      pop <- data.frame(rr, kk, ID=(1:N), sire=sort(rep(1:sires,dpsire)), Ag=rep(NA,N), Af=rep(NA,N))          
      
      Results_SimSel[[kk]]<-SimSel(pop, nReps, nGens, rr, kk-1, r.g, r.f, SigG.g, SigG.f, SigE.g, SigE.f, i.g, i.f, km.g, km.f,   
                                   sires, dpsire, N, gr.size, ngroups, Results_SimSel=Results_SimSel[[kk-1]])                                                                                         
      
      Results_SIR[[kk]]<-SIRdataGen(pop=Results_SimSel[[kk]]$pop, nReps, nGens, repl=rr, gen=kk-1, r.g, r.f, TRbeta, RRgamma,
                                    SigG.g, SigG.f, SigE.g, SigE.f, sires, dpsire, N, gr.size, ngroups, freqDep, reprnum, SIRdata)  
      
      
      if (save.pop) {        
        #####Saving Results_SIR_pop in data file
        file.pop=paste("Pop_", rr, "_", kk-1, "_", pfx, ".txt", sep="")
        Pop1<-do.call(cbind, Results_SIR[[kk]]$pop)
        suppressWarnings(write.table(Pop1, file=file.pop, append=F, row.names=F, col.names=F, quote=FALSE, sep = "    "))
      } 
      
      #####Saving Results_SIR_SIRts in data file
      Time1<-do.call(rbind, Results_SIR[[kk]]$SIRts)
      suppressWarnings(write.table(Time1, file=file.SIRtime, append=TRUE, col.names=!file.exists(file.SIRtime), row.names=F, quote=FALSE, sep = "    "))
      
      #####Generating data file with epidemic statistics
      Epi_stats<-list()
      
      #Total proportion of Infected individuals during the entire epidemic
      Infs <- length(which(Results_SIR[[kk]]$pop$Tinf!=-999))           
      Prop_inf <- Infs / N
      
      #Mean proportion across groups of the total # Infected within group
      Recs <- list()
      for (i in 1:ngroups) { Recs[[i]]<-(max(as.numeric(Results_SIR[[kk]]$SIRts[[i]]$R))) }
      Recs <- unlist(Recs)
      w <- which(Recs>1)
      if (length(which(w!=0))) {mean_prop_Recs <- mean(Recs[w]) / gr.size} else {mean_prop_Recs <- 1 / gr.size }
      
      #Mean proportion of groups with >1 Infected individuals, i.e. groups with epidemics
      Prop_epi <- length(w) / ngroups  
      
      #Mean proportion of groups that had epidemics classified as short, medium and long
      maxtime <- list()
      for (i in 1:ngroups) { maxtime[[i]]<-(max(as.numeric(Results_SIR[[kk]]$SIRts[[i]]$time))) }
      maxtime <- unlist(maxtime)
      Short_dur <- (length(which(maxtime<=180))) / ngroups
      Medium_dur <- (length(which(maxtime>180 & maxtime<=365))) / ngroups
      Long_dur <- (length(which(maxtime>365))) / ngroups
      
      #Generating the Inf_Stat data file
      Epi_stats<-as.data.frame(do.call(cbind, Results_SimSel[[kk]][1:8]))
      epi_descriptives<-cbind(Epi_stats, Infs, Prop_inf, mean_prop_Recs, Prop_epi, Short_dur, Medium_dur, Long_dur)        
      suppressWarnings(write.table(epi_descriptives, file=file.Inf_Stat, append=TRUE, col.names=!file.exists(file.Inf_Stat), row.names=F, quote=FALSE, sep = "    "))
      
      
    }}
}







