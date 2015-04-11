# MonteCarlo_genomicIntervals
<p>Monte Carlo simulation to calculate the probability of segment overlap in a random distribution of genomic intervals. </p>

The shuffling phase can be parallelized - select number of threads to equal or less than the number of individuals.

Dependencies:
  - Individual's bedfile containing the genomic intervals
  - refSeq gene locations
  - refseq Chromosome positions
<p>These files must be present in a MonteCarlo_genomicIntervals/data directory</p>


Usage:
<p><code>./ip_gene.py -np number_of_threads -nRuns number_of_steps -o output_file </code></p>


