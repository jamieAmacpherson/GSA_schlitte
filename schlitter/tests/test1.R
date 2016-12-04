#! /usr/bin/R
#===============================================================================
## compute the Schlitter entropy for a given trajectory
## (C) Jens Kleinjung 2016
#===============================================================================

#_______________________________________________________________________________
## fast covariance function
cov2 = function(x)
{
	# 'scale' is used here to centre columnes (subtract column mean)
	1/(NROW(x)-1) * crossprod(scale(x, TRUE, FALSE))
}

#_______________________________________________________________________________
## read trajectory
dat = read.table("debug_traj_C.dat");
## scale trajectory
dat.s = scale(dat, TRUE, FALSE);
write.table(dat.s, file = "debug_scale_R.dat", quote = FALSE, col.names = FALSE, row.names = FALSE);
## covariance matrix
dat.cov = cov(dat);
write.table(dat.cov, file = "covar_R.dat", quote = FALSE, col.names = FALSE);

## eigen system
dat.eigen =  eigen(dat.cov);
write.table(dat.eigen$value, file = "S_sch_eigenval_R.dat", quote = FALSE, col.names = FALSE);
write.table(dat.eigen$vector, file = "S_sch_eigenvec_R.dat", quote = FALSE, col.names = FALSE);

#_______________________________________________________________________________
## Schlitter entropy
## read physical constants
physcon = read.table("physcon.dat", header = TRUE);

## assign constants
k_B = physcon$"value"[physcon$"name" %in% "GSL_CONST_MKSA_BOLTZMANN"]; # J K^{-1}
T = 300; # K
e_sq = exp(1) * exp(1);
h_bar_sq = (physcon$"value"[physcon$"name" %in% "GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR"])**2; # J s
prefr = k_B * T * e_sq / h_bar_sq; # kg^{-1} m^{-2}
m_CH = 13.01864; # mass of unified C$^{\alpha}$H in au units
## unified atomic mass: conversion factor from au to kg
cf_au_kg = physcon$"value"[physcon$"name" %in% "GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS"];
cf_nmsq_msq = 1e-18; # conversion factor from nm^2 to m^2

## if trajectory too short, eigen values will be under-determined
eval_lim = min(length(dat.eigen$values), dim(dat)[1]);
S_sch = vector("numeric", length = eval_lim);
S_sch = sapply(dat.eigen$values[1:eval_lim], function(x) {
		return(log(1 + (prefr *  cf_au_kg * m_CH * x * cf_nmsq_msq)));
	}
)
write.table(S_sch, file = "S_sch_R.dat", quote = FALSE, col.names = FALSE);

S_sch.cumsum = cumsum(S_sch);
plot(1:eval_lim, S_sch.cumsum);
write.table(S_sch.cumsum, file = "S_sch_cumsum_R.dat", quote = FALSE, col.names = FALSE);

