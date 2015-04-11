# MonteCarlo_genomicIntervals
<p>Monte Carlo simulation to calculate the probability of overlap in a random distribution of genomic intervals. </p>
It can parallelize the shuffling phase - select number of threads to use. Should be equal or less than the number of individuals.

Usage:
<code>./ip_gene.py -np <number of threads> -nRuns <number of steps> -o <output file> </code>
