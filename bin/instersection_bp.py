import collections
import pylab
from bx.intervals.cluster import ClusterTree
 
 
num_overlaps=2 #minimal number of individuals in the same cluster
distance=-500000   #minimal distance for two reads to be considered in the same cluster. Set a negative value for requiring minimal overlaps
 
 
 
#attention: all parameters must be int, including read_id and match_id.
#match_id: chrm
#read_id: use a seperate dictionary to match int value with gene name
 
def your_parser(align_handle):
    while True:
        data = align_handle.readline()
        if not data:
            break
        try:
            yield int(data.split("\t")[3].rstrip()), int(data.split("\t")[0][3:]), int(data.split("\t")[1]), int(data.split("\t")[2]), int(0)
        except:
            pass
       
def parse_alignment(align_file, errors_allowed=0):
    with open(align_file) as align_handle:
        for (read_id, match_id, match_start, match_end, errors) in \
                your_parser(align_handle):
            if (errors <= errors_allowed):
                yield read_id, match_id, match_start, match_end
 
 
 
def build_cluster_trees(align_generator, cluster_distance):
    # arguments to ClusterTree are:
    # - Distance in basepairs for two reads to be in the same cluster;
    #   for instance 20 would group all reads with 20bp of each other
    # - Number of reads necessary for a group to be considered a cluster;
    #   2 returns all groups with 2 or more overlapping reads
    cluster_trees = collections.defaultdict(lambda:
            ClusterTree(cluster_distance, num_overlaps))
    for read_id, match_id, start, end in align_generator:
        cluster_trees[match_id].insert(start, end, read_id)
    
    return dict(cluster_trees)
 
 
 
 
if __name__ == "__main__": 
	path_to_bed="/Volumes/mac/intervalpermutation/tmp/shuffle2B.bed"

	#path_to_bed="/Volumes/mac/all_t.bed"
	cluster_trees=build_cluster_trees(parse_alignment(path_to_bed), distance)
 
#	print cluster_trees
	'''
	for chrm, cluster_tree in cluster_trees.items():
   		print chrm, cluster_tree
	'''    
 
 
	for chrom, cluster_tree in cluster_trees.items():
		for start, end, read_ids in cluster_tree.getregions():
			print chrom, start, end, read_ids
 
   
