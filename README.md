# GenEpiSim

Package developed by Smaragda Tsairidou 9.16-9.17 (Tsairidou et al. 2017)

SIR model code developed by Osvaldo Anacleto (Anacleto et al. 2015)

The package contains 3 functions:

1) Sel_DataGen.R: The main function.  
   It internally calls the other functions and returns the results. Specifically this program calculates response to selection for susceptibility and for infectivity, generates new populations after selection, and simulates epidemics following an SIR model. 
   It generates the following data files: "Inf_Stat_" with population and epidemic descriptive statistics and assessment measures of epidemic severity and duration; "SIRtime_" with number of Susceptible, Infected and Recovered individuals at time-points; and, "Pop_#repl_#gen_Susc_Inf_" if (save.pop=T), saving all data for each population in each generation and replicate. 
   The user can modify the parameter values when calling this function (see below). 

2) RespSelCalc_inf_susc.R: Function to calculate Response to selection for susceptibility and infectivity over generations. It calculates response to selection by calculating the change in the mean and variance of susceptibility and infectivity over generations. 

3) SimEpi_SIR.R: Generates data from a stochastc SIR model and produces the list with the SIR data frame as time series (SIRts), and on population level (pop). If (reprnum=T) only the infectivity of the index cases is accounted for in the calculation of the infection rate, and the output of the analysis can be used for the calculation of the reproduction number R0.    

Note: to run the main function a file containing seeds for each replicate is required in the working directory. This file is provided by the user and named "seed_data.txt".

Parameters in Sel_DataGen.R:

nReps = number of replications

nGens = number of generations

r.g = selection accuracy for susceptibility

r.f = selection accuracy for infectivity

TRbeta = average transmission coefficient beta

RRgamma = recovery rate gamma

SigG.g = genetic variance for susceptibility

SigG.f = genetic variance for infectivity

SigE.g = environmental variance for susceptibility

SigE.f = environmental variance for infectivity

im.g = selection intensity on the males for susceptibility

xm.g = deviation of the truncation point in the males for susceptibility

ifem.g = selection intensity on the females for susceptibility

xfem.g = deviation of the truncation point in the females for susceptibility

im.f = selection intensity on the males for infectivity

xm.f = deviation of the truncation point in the males for infectivity

ifem.f = selection intensity on the females for infectivity 

xfem.f = deviation of the truncation point in the females for infectivity   

sires = number of sires

dpsire = number of dams per sire

gr.size = group size

freqDep = T/F for frequensy or density dependent force of infection

reprnum = T/F if true then only the infectivity of the index cases is accounted for in the calculation of the infection rate 

save.pop = T/F if true then saves all population data for each generation and replicate

Example:

Sel_DataGen(nReps=5, nGens=2, r.g=0.7071, r.f=0.7071, TRbeta=0.02, RRgamma=0.017, SigG.g=0.5, SigG.f=0.5, SigE.g=2, SigE.f=2, im.g=0.798, xm.g=0, ifem.g=0, xfem.g=0, im.f=0.35, xm.f=-0.842, ifem.f=0, xfem.f=0, sires=200, dpsire=50, gr.size=100, freqDep=T, reprnum=T, save.pop=F)
