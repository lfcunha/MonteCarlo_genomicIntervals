# MonteCarlo_genomicIntervals
Monte Carlo simulation to calculate the probability of overlap in a random distribution of genomic intermals. 
It can parallelize the shuffling phase - select number pf threads to use. Should be equal or less than the number of individuals.

Usage:
./ip_gene.py -np <number of threads> -nRuns <number of steps> -o <output file>
