% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimRespCalc_inf_susc.R
\name{SimSel}
\alias{SimSel}
\title{Function to calculate Response to selection for susceptibility and infectivity over generations  (Smaragda Tsairidou 7.9.16)}
\usage{
SimSel(pop, nReps, nGens, rr, kk, r.g, r.f, SigG.g, SigG.f, SigE.g, SigE.f, i.g,
  i.f, km.g, km.f, sires, dpsire, N, gr.size, ngroups, Data_SimSel,
  Results_SimSel)
}
\arguments{
\item{pop}{population data}

\item{nReps}{number of replicates}

\item{nGens}{number of generations}

\item{rr}{current replicate}

\item{kk}{current generation}

\item{r.g}{selection accuracy for susceptibility}

\item{r.f}{selection accuracy for infectivity}

\item{SigG.g}{genetic variance for susceptibility}

\item{SigG.f}{genetic variance for infectivity}

\item{SigE.g}{environmental variance for susceptibility}

\item{SigE.f}{environmental variance for infectivity}

\item{i.g}{selection intensity for susceptibility}

\item{i.f}{selection intensity for infectivity}

\item{km.g}{the Bulmer effect factor for susceptibility}

\item{km.f}{the Bulmer effect factor for infectivity}

\item{sires}{number of sires}

\item{dpsire}{number of dams per sire}

\item{N}{population size}

\item{gr.size}{group size}

\item{ngroups}{number of groups}

\item{Data_SimSel}{data from selection}

\item{Results_SimSel}{results from selection}
}
\value{
Dataframe with new population parameters
}
\description{
Calculates response to selection by calculating the change in the mean and variance of susceptibility and infectivity over generations
}
