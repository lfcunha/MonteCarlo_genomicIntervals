#!/usr/bin/env python
 
__author__ = 'Luis Cunha'
 
 
import sys
import os
sys.path.append(os.path.abspath("bin"))
from bx.intervals.intersection import IntervalTree
import time
import cPickle as pickle
import random
import collections
from multiprocessing import Pool
import readIntervals as ri
import json
import argparse
 
 
"""
Perform Monte Carlo simulation to permutate a list of genomic segments and calculate the proportion of intervals
overlapping the same gene.
 
Input:
   number of steps;
   segments file in bed format must be present in the data folder
 
Output:
   file with the result of each simulation per line
   result as json describing occurrence or not (1 / 0) of each number of overlaps
   e.g: observed 5, 6, 7, 9, 11 overlaps  {0:0, 1:0, 2:0, 3:0, 4:0, 5:1, 6:1, 7:1, 8:0, 9:1, 10:0, 11:1, 12:0}
 
Usage:
   intervalpermutation -n <number of steps> -o <output file>
 
"""
 
 
 
#Globals
 
EXOME_ARRAY_LENGTH = 74494878 - 1  # length of the gene position array. Value used to generate random number in this range
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REFSEQ_POSITIONS_ARRAY = BASE_DIR + "/data/positionArrayRefSeq_chrsizeX.bed"
REFSEQ_GENES_BEDFILE = BASE_DIR + "/data/refseq_genes.bed"
REFSEQ_POSITION_ARRAY = "positionArrayRefSeqX.p"
CHRM_SIZE_REFSEQ_POSITION_ARRAY = "positionArrayRefSeq_chrsizeX.p"
COUNT_SHUFFLE = {}                                #keeps track of the distribution of random permutations per chromosome
COUNT = {}
 
 
 
class IntervalTreeDict(object):
    """Build an Interval Tree data structure of gene position ranges
   """
 
    def __init__(self, bed_file):
        """
       :param bed_file:
       :return interval tree data structure of gene ranges:
       """
        self.interval_tree_dict = dict()
        error_message = "Skipping line {0} too short with only {1} column(s).Check that it conforms to bed format:\n{2}"
 
        with open(bed_file, 'r') as bed_file_Handler:
            for count, line in enumerate(bed_file_Handler):
                segmentProperties = line.split("\t")
                numberOfColumns = len(segmentProperties)
                try:
                    chromosome, segment_start, segment_end, name = segmentProperties[:4]
                except ValueError:
                    print error_message.format(count + 1, numberOfColumns, line.strip())
                    continue
 
                segment_start, segment_end = int(segment_start), int(segment_end)
 
                if chromosome in self.interval_tree_dict:
                    tree = self.interval_tree_dict[chromosome]
                else:
                    tree = IntervalTree()
                    self.interval_tree_dict[chromosome] = tree
 
                tree.add(segment_start, segment_end, tuple(segmentProperties[:4]))
 
    @property
    def interval_tree(self, chromosome, start, end):
        return self.interval_tree_dict[chromosome].find(start, end)
 
 
class InitializeIntervals:
    """Read:
       - the RefSeq gene position and gene size arrays from disk to be used on each simulation step
       - the dictionary of intervals of all individuals, from the pickle file on disk. If the file does not exits,
       - all intervals file. Pickle dump it to disk to be used in subsequent runs. If file not exist (first time usage)
       call the readIntervals module (in the bin subdirectory) and process individual interval files into a single file.
   """
 
    def __init__(self):
        self.positionArray = pickle.load(open(BASE_DIR + '/data/'+ REFSEQ_POSITION_ARRAY, 'rb'))
        self.chrArray = pickle.load(open(BASE_DIR + '/data/' + CHRM_SIZE_REFSEQ_POSITION_ARRAY, 'rb'))
        self.individuals = self.load_intervals_all_individuals()
 
    def load_intervals_all_individuals(self):
        """
       :return all interval dictionary:
       """
        if not os.path.isfile(BASE_DIR + "/data/allintervals.p"):
            if not os.path.exists("data"):
                os.makedirs("data")
            readIntervals = ri.readIntervals()
            individual = readIntervals.readIntervals()
            intervals = {}
            for id in individual:
                intervals[id] = []
                with open(BASE_DIR + "/data/" + id + ".id.bed", "r") as bef_file_handler:
                    lines = bef_file_handler.readlines()
                    for line in lines:
                        l = line.split()
                        intervals[id].append({l[0]: int(l[2]) - int(l[1])})
            pickle.dump(intervals, open(BASE_DIR + "/data/allintervals.p", "wb"))
            return intervals
        else:
            return pickle.load(open("data/allintervals.p", "rb"))
 
    @property
    def position_array(self):
        return self.positionArray
 
    @property
    def chr_array(self):
        return self.chrArray
 
    @property
    def _individuals(self):
        return self.individuals
 
 
def find_chr_from_array_position(position):
    """
   :param position:
   :return chromosome:
   """
    result = refSeq_positions_dictionary.interval_tree("1", position, position)
    try:
        chromosome = int(result[0][3].strip())
        return chromosome
    except KeyError:
        print "Position out of range: ", position
        return 0
 
 
def test_shuffled_interval_overlap(intervals):
    """Test that the shuffled intervals do not overlap each other
   If there is a single overlap, discard this while shuffle step and redo
   (discarding only this interval would introduce a bias in the probability of the position and it would not be a
   purely random shuffle)
   :param intervals:
   :return: True if no overlap, False if overlap
   """
    result = {}
    for interval in intervals.values()[0]:
        try:
            chromosome = interval[0]
            if chromosome not in result:
                result[chromosome] = {}
            result[chromosome][interval[1]] = interval[2]
        except:
            pass                                   #Do not interrupt due to any exception. Continue to the next interval
 
    for start, end in result:
        orderedItems = collections.OrderedDict(sorted(end.items()))
        start_list = []
        end_list = []
        for interval, end_ in orderedItems.items:
            start_list.append(interval)
            end_list.append(end_)
        for interval in range(0, len(a) - 1):
            if start_list[interval + 1] < end_list[interval]:
                print "Overlapped Intervals: reject and repeat shuffle\n"
                return False
    return True
 
 
def shuffle(individual_intervals):
    """Shuffle the position of every interval:
      - picks a random position along the length of the Positions Array
      - gets the exome coordinate from this position. Note that intervals are only allowed to be shuffled into coding
      regions of the exome. Specifically, those regions covered by the capture probes
      - intervals of each individual are shuffled together, with each individual shuffled in parallel in a separate
       thread
 
   :param individual_intervals:
   :return shuffled intervals:
   """
    print "running parallel process"
    random.seed()                         #set seed on each thread to avoid same random number generation on the threads
    time.sleep(random.randint(0, 9))
    intervals = {individual_intervals[0]: []}
 
    for interval in individual_intervals[1]:
        randomPosition = random.randint(0, EXOME_ARRAY_LENGTH)
        start = positionArray[randomPosition]
        end = start + int(interval.values()[0])
        chromosome = find_chr_from_array_position(randomPosition)
        intervals[individual_intervals[0]].append((chromosome, start, end))
        if chromosome not in COUNT_SHUFFLE: COUNT_SHUFFLE[chromosome] = 0
        COUNT_SHUFFLE[chromosome] += 1
 
    if test_shuffled_interval_overlap(intervals):                                            #if segments do not overlap
        """log the distribution of random intervals"""
        with open(BASE_DIR + "tmp/shuffled_intervals.json", "a") as outfile:
            json.dump(intervals, outfile)
            outfile.write("\n")
        return intervals
    else:
        shuffle(individual_intervals)
 
def inc(i):
    """Generator to increment id used in start_shuffle
   :param i:
   :return:
   """
    yield i + 1
 
 
def start_shuffle(output_file):
    """Initiate parallel threads for interval shuffling phase
   :param output_file:
   :return 0:
   """
    index=0
    starttime = time.time()
    individualIntervals = allIndividuals.items()
    try:
        print "starting parallel shuffle"
        pool = Pool(8)
        results = pool.map(shuffle, individualIntervals)
        pool.close()
        pool.join()
    except:
        os.nice(100)
        pass
    else:
        print "finished shuffling phase. Starting overlap analysis"
        elapsedtime = time.time() - starttime
        reads = {}
        persons_reads = {}
        for result in results:
            for y in result.values()[0]:
                id = str(next(inc(index)))
                reads[id] = [str(y[0]), str(y[1]), str(y[2]), str(result.keys()[0])]
                if str(result.keys()[0]) not in persons_reads: persons_reads[str(result.keys()[0])] = []
                persons_reads[str(result.keys()[0])].append(id)
                index += 1
 
        """Dictionary to keep track of occurrence of each number of overlaps: 0/1 (no/yes)"""
        local_overall_overlaps = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0,
                                  13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0}
 
        """for each interval, of each individual, get the genes in this region from the tree representation of refseq
       genes(refseq_gene_tree) and build a dictionary of gene:[list of intervals covering this gene]
       finally COUNT the number of intervals covering the gene. This is in number of intervals overlaping a gene
       """
        genes = {}
        for read in reads:
            l = reads[read]
            a = refseq_gene_tree.interval_tree(l[0], int(l[1]), int(l[2]))
            for result in a:
                b = result[3][:-1]
                if b not in genes:
                    genes[b] = []
                genes[b].append(l[3])
 
        for result in genes:
            if len(genes[result]) > 1:
                if (len(genes[result])) not in local_overall_overlaps:
                    local_overall_overlaps[len(genes[result])] = 0
                if local_overall_overlaps[len(genes[result])] == 0:
                    local_overall_overlaps[len(genes[result])] = 1
 
        with open(output_file, 'a') as outfile:
            json.dump(local_overall_overlaps, outfile)
            outfile.write("\n")
 
        print "Finished in {0:.1f}".format(elapsedtime) + " s"
 
    return 0
 
 
 
"""load global arrays and intervals"""
 
print "Loading refseq trees"
refSeq_positions_dictionary = IntervalTreeDict(REFSEQ_POSITIONS_ARRAY)
refseq_gene_tree = IntervalTreeDict(REFSEQ_GENES_BEDFILE)
print "Loading Individuals, positionsArray, and chromosomeArray"
initialize = InitializeIntervals()
allIndividuals = initialize._individuals
positionArray = initialize.positionArray
chrArray = initialize.chrArray
 
 
def usage():
    print "shuffle python shuffle5.py -n <number of steps> -o <output file>"
    sys.exit(0)
 
 
def main():
    parser = argparse.ArgumentParser(
        description='Monte Carlo simulation to permutate genomic intervals and calculate the proportion of segments overlapping the same gene')
    parser.add_argument('-o', '--output', help='output file')
    parser.add_argument('-v', dest='verbose', action='store_true')
    parser.add_argument('-nRuns', '--numberRuns', help='number of simulation steps')
    args = parser.parse_args()
    nRuns = int(args.numberRuns)
    output = BASE_DIR + "/result/" + args.output
 
    print "Starting...\n"
 
    while nRuns > 0:
        try:
            print "Iteration " + str(nRuns)
            start_shuffle(output)
            print "finished iteration " + str(nRuns) + "\n"
            nRuns -= 1
        except:                                          #If any exception occurs, continue with next step of simulation
            time.sleep(1)
            print "exception raised"
            start_shuffle(output)
 
 
if __name__ == "__main__":
    main()
