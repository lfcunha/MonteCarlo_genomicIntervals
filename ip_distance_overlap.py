import time
import cPickle as pickle
import os
import random
import sys
import collections
from multiprocessing import Pool
import scipy
sys.path.append(os.path.abspath("bin"))
import chr_position_tree as tree
import instersection_bp as Ibp
import json
from max_overlap import max_overlap

arrayLength=74494878-1
distance=-1500000   #minimum overlap


path_to_bed="data/positionArrayRefSeq_chrsizeX.bed"
dic=tree.IntervalTreeDict(path_to_bed)


overall_overlaps={0:0, 1:0, 2:0, 3:0, 4:0,5:0,6:0,7:0, 8:0,9:0,10:0,11:0,12:0,13:0,14:0,15:0,16:0,17:0,18:0,19:0,20:0,21:0,22:0,23:0}

count_shuffle={}
count={}
total=0
numsim=0





class Initialize:
	def __init__(self):
		pass	
	def load_intervals_all_individuals():
        	if not os.path.isfile("data/allintervals.p"):
                	#ind = ['04764', '05375','04860','03552','03476', '02318','04661','00688','01970','05118','04911','03121','05476','04924','04923','03107','00287', '03765', '04710']
			ind = ['01498', '01322','05059','01687','04334', '00161','02315']
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
	


	positionArray = pickle.load(open('data/positionArrayRefSeqX.p', 'rb'))
	chrArray = pickle.load(open('data/positionArrayRefSeq_chrsizeX.p', 'rb'))
	individuals = load_intervals_all_individuals()


initialize = Initialize()
allIndividuals = initialize.individuals
positionArray = initialize.positionArray
chrArray = initialize.chrArray




def find_chr_from_array_position(pos):
	result= dic.find("1", pos, pos)
	#print "start: ", pos, result
	try:
		#print int(result[0][3].strip())
		return int(result[0][3].strip())
	except:
		print "error in find_chr_from_array_position function; pos: ", pos
		return 0

def test_non_overlap(intervals):
    #print intervals
    #print intervals.keys()
    #print intervals
    #print "**\n"

    res={}
    for x in intervals.values()[0]:
	try:
            chr=x[0]
            if chr not in res:
                res[a]={}
            res[a][x[1]]=x[2]
        except:
            pass

    for f,y in res:
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
	random.seed()
	time.sleep(random.randint(0,9))
	intervals ={}
	intervals[individual_intervals[0]]=[]
	
	for x in individual_intervals[1]:
		randomPosition=random.randint(0,arrayLength)
		start = positionArray[randomPosition]
		end = start+int(x.values()[0])
		chr = find_chr_from_array_position(randomPosition)
		intervals[individual_intervals[0]].append((chr, start, end))
		if chr not in count_shuffle: count_shuffle[chr]=0
		count_shuffle[chr]+=1

                with open('tmp/random_positions.txt', "a") as outfile:
                        outfile.write(str(randomPosition)+"\t"+str(chr)+"\t"+str(start)+"\t"+str(end)+"\n")

        if test_non_overlap(intervals):

                #write intervals to file to later validate result
                with open('tmp/shuffled_intervals.json', "a") as outfile:
                        json.dump(intervals, outfile)
			outfile.write("\n")
                return intervals
        else:
                shuffle(individual_intervals)


	

def start_shuffle(n):
	#g.write(str(n))
	#g.write("\n")
	i=0
	startTime  =  time.time()
	individualIntervals = allIndividuals.items()
	try:
		pool = Pool(8)  #or os.nice(100)
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
		elapsedTime  =  time.time() - startTime
		f=open("tmp/shuffle.bed", "w")
		reads={}
		persons_reads={}
		for x in result:
			for y in x.values()[0]:
				s=""
				id=str(next(inc(i)))
				s+="chr"+str(y[0])+"\t"+str(y[1])+"\t"+str(y[2])+"\t"+id+"\t"+str(x.keys()[0])
				f.write(s)
				f.write("\n")
				reads[id]=[str(y[0]),str(y[1]),str(y[2]),str(x.keys()[0])]
				if str(x.keys()[0]) not in persons_reads: persons_reads[str(x.keys()[0])]=[]
				persons_reads[str(x.keys()[0])].append(id)
				i+=1

		f.close()			

		cluster_trees=Ibp.build_cluster_trees(Ibp.parse_alignment("tmp/shuffle.bed"),distance)
		clusters={}
		local_overall_overlaps={0:0, 1:0, 2:0, 3:0, 4:0,5:0,6:0,7:0, 8:0,9:0,10:0,11:0,12:0,13:0,14:0,15:0,16:0,17:0,18:0,19:0,20:0,21:0,22:0,23:0}
		for chrom, cluster_tree in cluster_trees.items():
			
			if chrom not in count: count[chrom]={}
    			for start, end, read_ids in cluster_tree.getregions():
        			read_info=[]			
				clus={}
				s=[]
				e=[]
				for id in set(read_ids):
					s.append(reads[str(id)][1])
					e.append(reads[str(id)][2])
					read_info.append((reads[str(id)][1], reads[str(id)][2]))
				cluster_name=str(chrom)+"_"+str(start)
				with open("clusters.txt", "a") as out:
					
					clus[cluster_name]=read_info		

					json.dump(clus, out) 
					out.write("\n")
				
				ol=max_overlap(s,e, -distance)
				for x in ol:
					if x[0] not in count[chrom]: count[chrom][x[0]]=0
					if x[0] not in overall_overlaps: overall_overlaps[x[0]]=0
					if local_overall_overlaps[x[0]]==0: local_overall_overlaps[x[0]]=1
					count[chrom][x[0]]+=1

				for item in ol:
					it=map(str,item)
					line= ",".join(str(it))

		#accumulate overlaps
		#for x in local_overall_overlaps: overall_overlaps[x]+= local_overall_overlaps[x]
                #total+=sum([overall_overlaps[x] for x in overall_overlaps])
                #numsim+=1
		with open('tmp/result.txt', 'a') as outfile:
			json.dump(local_overall_overlaps, outfile)
			outfile.write("\n")
	
                with open('tmp/result_per_chrm.txt', 'a') as outfile2:
                        json.dump(count, outfile2)
                        outfile2.write("\n")	
	
		print numsim
		#print total
		#print overall_overlaps
		print elapsedTime
		
#overall_overlaps={0:0, 1:0, 2:0, 3:0, 4:0,5:0,6:0,7:0, 8:0,9:0,10:0,11:0,12:0,13:0,14:0,15:0,16:0,17:0,18:0,19:0,20:0,21:0,22:0,23:0}

#count_shuffle={}
#count={}
#total=0
#numsim=0

n=4030

while(n<10000):
	try:
		start_shuffle(n)
		n+=1
#		with open('tmp/shuffled_intervals.json', "a") as outfile:
#			outfile.write("sim: " + n +"\n")
	except:
		time.sleep(1)
		print "exception raised"
		start_shuffle(n)

with open('tmp/data2.json', "w") as outfile:
	json.dump(count, outfile)
print count_shuffle


#g.close()


# add chrmosome count of shuffling
# add count of results (number of overlaps)
# add count of chromosome with overlaps
# add a way to visualize distribution of shuffles per chromosome, with a manhattan plot or similar. maybe make bins
