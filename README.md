# SIR-patch

Code for solving master equation for the SIR model in an arbitrary number of patches/subpopulations.
Adapted from [Nguyen & Carlson (2016)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0152950), which solves the master equation for 2 interacting populations, and which is in turn adapted from [Jenkinson & Goutsias (2012)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0036160).


File 'SIR_patch.m' solves the master equation given user-specified values for:
Npatch: number of patches

s0: initial numbers of susceptible individuals (Npatch x 1 vector)

i0: initial numbers of infected individuals (Npatch x 1 vector)

rep: reproductive number of each patch (Npatch x 1 vector)

f: symmetric, doubly-stochastic matrix specifying fractions of home (diagonal terms) and foreign (off-diagonal) contacts (Npatch x Npatch matrix)


SIR_patch.m creates three auxiliary functions if they do not already exist:

genMatrix_$NPATCH: generates transition matrix

propensity_$NPATCH: calculates propensities for each reaction

mapProbs_$NPATCH: maps probability vector into probability distribution of # susceptible, # infected for each patch