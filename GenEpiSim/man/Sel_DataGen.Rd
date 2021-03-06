% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sel_DataGen.R
\name{Sel_DataGen}
\alias{Sel_DataGen}
\title{Main function (developped by Smaragda Tsairidou 7.9.16)}
\usage{
Sel_DataGen(nReps, nGens, r.g, r.f, TRbeta, RRgamma, SigG.g, SigG.f, SigE.g,
  SigE.f, im.g, xm.g, ifem.g, xfem.g, im.f, xm.f, ifem.f, xfem.f, sires, dpsire,
  N = sires * dpsire, gr.size = N/ngroups, ngroups = N/gr.size, freqDep,
  reprnum, save.pop)
}
\arguments{
\item{nReps}{number of replicates}

\item{nGens}{number of generations}

\item{r.g}{selection accuracy for susceptibility}

\item{r.f}{selection accuracy for infectivity}

\item{TRbeta}{average transmission coefficient beta}

\item{RRgamma}{recovery rate gamma}

\item{SigG.g}{genetic variance for susceptibility}

\item{SigG.f}{genetic variance for infectivity}

\item{SigE.g}{environmental variance for susceptibility}

\item{SigE.f}{environmental variance for infectivity}

\item{im.g}{selection intensity on the males for susceptibility}

\item{xm.g}{truncation point for the males for susceptibility}

\item{ifem.g}{selection intensity on the females for susceptibility}

\item{xfem.g}{truncation point for the females for susceptibility}

\item{im.f}{selection intensity on th emales for infectivity}

\item{xm.f}{truncation point on the males for infectivity}

\item{ifem.f}{selection intensity on the females for infectivity}

\item{xfem.f}{truncation point for the females for infectivity}

\item{sires}{number of sires}

\item{dpsire}{number of dams per sire}

\item{N}{population size}

\item{gr.size}{group size}

\item{ngroups}{number of groups}

\item{freqDep}{T/F for frequensy or density dependent force of infection}

\item{reprnum}{T/F for modified calculation of infection rate to calculate reproduction number R0}

\item{save.pop}{T/F to save all population data for each generation and each replicate}
}
\value{
Generates data files: "Inf_Stat_" with population and epidemic descriptive statistics and assessment measures of epidemic severity and duration
                              "SIRtime_" with number of S, I and Rs at time-points 
                              "Pop_#repl_#gen_Susc_Inf_" if save.pop=T saving all data for each population in each generation and replicate
}
\description{
Calculates response to selection for susceptibility and infectivity and generates new population (function "SimSel" created by S. Tsairidou), and simulates epidemics based on an SIR model (function "SIRdataGen" created by O. Anacleto)
}
