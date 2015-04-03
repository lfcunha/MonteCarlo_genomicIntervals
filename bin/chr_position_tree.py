from bx.intervals.intersection import IntervalTree
import sys

 
class IntervalTreeDict(object):
    def __init__(self, bed_file_path):
 
        self.interval_tree_dict = dict()
        error_message = "Skipping line {0} - too short, only {1} column(s):\n{2}"
 
        with open(bed_file_path, 'r') as bed_file:
            for count, line in enumerate(bed_file):
                split_line = line.split("\t")
                number_of_columns = len(split_line)
 
                try:
                    chromosome, start, end, name = split_line[:4]
                except ValueError:
                    print error_message.format(count+1, number_of_columns, line.strip())
                    continue
                   
                start, end = int(start), int(end)
                tree = None
 
                if chromosome in self.interval_tree_dict:
                    tree = self.interval_tree_dict[chromosome]
                else:
                    tree = IntervalTree()
                    self.interval_tree_dict[chromosome] = tree
 
                tree.add(start, end, tuple(split_line[:4]))
 
    def find(self, chromosome, start, end):
 
        return self.interval_tree_dict[chromosome].find(start, end)
 

if __name__ == "__main__":
	path_to_bed="data/positionArrayRefSeq_chrsizeXB.bed"
	dic=IntervalTreeDict(path_to_bed)
	print dic.find("1", int(sys.argv[1]), int(sys.argv[1]))
