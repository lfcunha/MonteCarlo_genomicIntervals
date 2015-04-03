import os, cPickle as pickle


class readIntervals(object):
	def __init__(self):
		self.individualIds=[]
	def readIntervals(self):
		f=open("inputIntervals.bed", "r")
		lines=f.readlines()

		oldId=""
		intervals=[]
		if not os.path.exists("data"):
			os.makedirs("data")

		def wr(intervals, id):
			name=id+".id.bed"
			g=open("data/" +name, "a")
			for x in intervals:
				g.write(str(x))
	       			g.write("\n")
    			g.close()
		
		for line in lines:
        		l=line.split()
        		id=l[-1]
        		if oldId==id:
        			if id not in self.individualIds:
        				intervals.append(line.rstrip())
				self.individualIds.append(id)
        		else:	
        			wr(intervals, id)
        			oldId=id
        			intervals=[]

			wr(intervals, id)
		
		pickle.dump(self.individualIds, open("data/IDs.p", "wb" ))
		return self.individualIds



if __name__ == "__main__":
        a=readIntervals()
	b=a.readIntervals()
	print b
