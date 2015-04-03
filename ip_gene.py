import sys,os
sys.path.append(os.path.abspath("bin"))
from bx.intervals.intersection import IntervalTree
import time
import cPickle as pickle
import argparse
import random
import collections
from multiprocessing import Pool
import scipy
#import chr_position_tree as tree
import instersection_bp as Ibp
import readIntervals as ri
import json
from max_overlap import max_overlap
import argparse
import getopt


arrayLength=74494878-1 #length of the gene position array. Value used to generate random number in this range
baseDir=os.path.dirname(os.path.abspath(__file__))
refSeq_positions_Array=baseDir + "/data/positionArrayRefSeq_chrsizeX.bed"
refSeq_genes_bedFile = baseDir + "/data/refseq_genes.bed"
count_shuffle={}
count={}


class IntervalTreeDict(object):
    def __init__(self, bed_file):
        self.interval_tree_dict = dict()
        error_message = "Skipping line {0} - too short, only {1} column(s):\n{2}"

        with open(bed_file, 'r') as bed_file_Handler:
            for count, line in enumerate(bed_file_Handler):
                split_line = line.split("\t")
                number_of_columns = len(split_line)
                try:
                    chromosome, start, end, name = split_line[:4]
                except ValueError:
                    print error_message.format(count+1, number_of_columns, line.strip())
                    continue

                start, end = int(start), int(end)
                #tree = None

                if chromosome in self.interval_tree_dict:
                    tree = self.interval_tree_dict[chromosome]
                else:
                    tree = IntervalTree()
                    self.interval_tree_dict[chromosome] = tree

                tree.add(start, end, tuple(split_line[:4]))
    @property
    def interval_tree(self, chromosome, start, end):
        return self.interval_tree_dict[chromosome].interval_tree(start, end)


class Initialize(object):
	def __init__(self):
		self.positionArray = pickle.load(open(baseDir+'/data/positionArrayRefSeqX.p', 'rb'))
        	self.chrArray = pickle.load(open(baseDir+'/data/positionArrayRefSeq_chrsizeX.p', 'rb'))
        	self.individuals = self.load_intervals_all_individuals()	
	def load_intervals_all_individuals(self):
        	if not os.path.isfile("data/allintervals.p"):
			if not os.path.exists("data"):
    				os.makedirs("data")
                	readIntervals=ri.readIntervals()
			ind=readIntervals.readIntervals()
                	intervals = {}
                	for id in ind:
                        	intervals[id] = []
                        	f = open("data/" + id + ".id.bed", "r")
                        	lines = f.readlines()
                        	f.close()
                        	for line in lines:
                                	l = line.split()
                                	intervals[id].append({l[0]:int(l[2])-int(l[1])})
                	pickle.dump(intervals, open('data/allintervals.p', 'wb'))
                	return intervals
        	else:
                	return pickle.load(open('data/allintervals.p', 'rb'))

    	@property
    	def position_array(self):
        	return self.positionArray

    	@property
    	def chr_array(self):
        	return self.chrArray

    	@property
    	def _individuals(self):
        	return self.individuals


def find_chr_from_array_position(pos):
	result= refSeq_positions_dictionary.interval_tree("1", pos, pos)
	try:
		return int(result[0][3].strip())
	except:
		print "error in find_chr_from_array_position function; pos: ", pos
		return 0

def test_non_overlap(intervals):
    print "test called"
    res={}
    for x in intervals.values()[0]:
	try:
            chr=x[0]
            if chr not in res:
                res[chr]={}
            res[chr][x[1]]=x[2]
        except:
            pass

    for f,y in res:
	print f, y
    	od = collections.OrderedDict(sorted(y.items()))
    	print od
    	a=[]
    	b=[]
    	for x,y in od.items:
		a.append(x)
		b.append(y)
    	for x in range(0, len(a)-1):
		if a[x+1]<b[x]:
			print "reject"
			return False		
    return True


def inc(i):
	yield i+1


def shuffle(individual_intervals):
	print "running parallel process " # str(inc(i))
	random.seed()
	time.sleep(random.randint(0,9))
	#intervals ={}
	#intervals[individual_intervals[0]]=[]
	intervals = {individual_intervals[0]: []}

	for x in individual_intervals[1]:
		randomPosition=random.randint(0,arrayLength)
		start = positionArray[randomPosition]
		end = start+int(x.values()[0])
		chr = find_chr_from_array_position(randomPosition)
		intervals[individual_intervals[0]].append((chr, start, end))
		if chr not in count_shuffle: count_shuffle[chr]=0
		count_shuffle[chr]+=1
		#this file can be use to plot the distribution of random positions
                #with open('tmp/random_positions.txt', "a") as outfile:
                #        outfile.write(str(randomPosition)+"\t"+str(chr)+"\t"+str(start)+"\t"+str(end)+"\n")
	#test that no interval overlaps another. If so, toss the whole permutation and start over
    	print "calling test\n"
        if test_non_overlap(intervals):
		print "done with test\n"
                #write intervals to file to later validate result
                #with open('tmp/shuffled_intervals.json', "a") as outfile:
                #       json.dump(intervals, outfile)
		#	outfile.write("\n")
                return intervals
        else:
                shuffle(individual_intervals)


def start_shuffle(output):
	i=0
	startTime  =  time.time()
	individualIntervals = allIndividuals.items()
	try:
		print "starting parallel shuffle"
		pool = Pool(8)  
		result = pool.map(shuffle, individualIntervals) 
		pool.close()
		#pool.terminate()
		pool.join()	
	except:
            #if str(ex) == "[Errno 35] Resource temporarily unavailable":
            #time.sleep(0)
	    os.nice(100)    
            pass
	else:	
		print "finished shuffling phase. Starting overlap analysis"	
		elapsedTime  =  time.time() - startTime
		#f=open("tmp/shuffle.bed", "w")
		reads={}
		persons_reads={}
		for x in result:
			for y in x.values()[0]:
				s=""
				#id=str(next(inc(i)))
				id=str(i)
				print "id", str(id)
				#s+="chr"+str(y[0])+"\t"+str(y[1])+"\t"+str(y[2])+"\t"+id+"\t"+str(x.keys()[0])
				#f.write(s)
				#f.write("\n")
				reads[id]=[str(y[0]),str(y[1]),str(y[2]),str(x.keys()[0])]
				if str(x.keys()[0]) not in persons_reads: persons_reads[str(x.keys()[0])]=[]
				persons_reads[str(x.keys()[0])].append(id)
				i+=1
		#f.close()		
		#dictionary to hold the count of of occurrence of each number of overlaps	
		local_overall_overlaps={0:0, 1:0, 2:0, 3:0, 4:0,5:0,6:0,7:0, 8:0,9:0,10:0,11:0,12:0,13:0,14:0,15:0,16:0,17:0,18:0,19:0,20:0}
		f=open("tmp/shuffle.bed", "r")
		lines=f.readlines()
		f.close()
		genes={}
		#for each interval, of each individual, get the genes in this region from the tree representation of refseq genes(refseq_gene_tree)
		#build a dictionary of gene:[list of intervals covering this gene]
		#finally count the number of intervals covering the gene. This is in number of intervals overlaping a gene
		

#		for line in lines:
#			l=line.split()
			#a=refseq_gene_tree.find(l[0][3:], int(l[1]), int(l[2]))
		for read in reads:
			l=reads[read]
			a=refseq_gene_tree.find(l[0], int(l[1]), int(l[2]))
			for x in a:
				b= x[3][:-1]
				if b not in genes:
					genes[b]=[]
				genes[b].append(l[3])
				#genes[b].append(l[4])				
		for x in genes:
        		if len(genes[x])>1:
				if (len(genes[x])) not in local_overall_overlaps:
					local_overall_overlaps[len(genes[x])]=0
				if local_overall_overlaps[len(genes[x])]==0:local_overall_overlaps[len(genes[x])]=1

		with open(output, 'a') as outfile:
			json.dump(local_overall_overlaps, outfile)
			outfile.write("\n")

		print "Finished in {0:.1f}".format(elapsedTime) + " s" 
		
#global variables holding the tree data structure for refseq positions and refseq genes

print "Loading refseq trees"
refSeq_positions_dictionary=IntervalTreeDict(refSeq_positions_Array)
refseq_gene_tree=IntervalTreeDict(refSeq_genes_bedFile)
#load individual's intervals, and set up gene and chromosome arrays
print "Loading Individuals, positionsArray, and chromosomeArray"
initialize = Initialize()
allIndividuals = initialize._individuals
positionArray = initialize.position_array
chrArray = initialize.chr_array

def usage():
        print "shuffle python shuffle5.py -n 110 -o result_genes.txt [--ovlpSize overlap_size]"
        sys.exit(0)

def main():
        parser = argparse.ArgumentParser(description='Permutate genomic segments and calculate # overlaps with minimal size provided as argument')
	parser.add_argument('--ovlpSize',  help='minimal length of overlap')
        parser.add_argument('-o', '--output')
        parser.add_argument('-v', dest='verbose', action='store_true')
        parser.add_argument('-n', '--numberRuns')
	args = parser.parse_args()
	n = int(args.numberRuns)
	output = "result/"+args.output
        basedir=os.path.dirname(os.path.abspath(__file__))
	workDir=args.output

	print "Starting...\n"
	while(n>0):
		try:
			print "Iteration " + str(n)
			start_shuffle(output)
			print "finished iteration " + str(n) +"\n" 
			n-=1
		except:
			time.sleep(1)
			print "exception raised" 
			start_shuffle(output)

if __name__ == "__main__":
	main()


# add chromosome count of shuffling
# add count of results (number of overlaps)
# add count of chromosome with overlaps
# add a way to visualize distribution of shuffles per chromosome, with a manhattan plot or similar. maybe make bins
