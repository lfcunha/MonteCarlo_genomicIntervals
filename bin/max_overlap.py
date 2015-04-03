"""
given two arrays of start and end position of a cluster of segments, find number of discrete overlaps
returns a list of lists of (number of segments in overlap, start position, end position)

>>> max_overlap([10,20,25,30,50,80],[35,40,60,70,90,100], 1)
[[4, 30, 35], [3, 50, 60], [2, 80, 90]]

"""


def max_overlap(s1,e1, overlap_size=50000):
    

    #convert start and end strings to ints
    e_=map(int, e1)
    s_=map(int, s1)
    
    #sort start and end arrays
    e=sorted(e_)
    s=sorted(s_)

    #make a deep deep copy of the array
    e2=e[:]
    s2=s[:]



    start=[]
    ends=[]
    lv=0        # local variable to count overlaps
    flag=1	# flag to decide when to quit
    maxi=0	# max overlaps
    maxis=[]	# array of mx overlaps:  count multiple seperate overlaps within the same cluster
    lv2=[]	# array of localvariable count
    ol=[]	# array of overlaps. will contain: (number of segments overlapping, start, end)
    while(flag==1):
        if len(s)>0:
            if s[0]<e[0]:
                lv+=1
                # if lv>maxi: maxi=lv
                st=s.pop(0)
                if lv>maxi:
                    maxi=lv
            else:
                lv-=1
                end=e[0]
                e.pop(0)
        else:
            if len(e)>0:
                lv-=1
                if lv<maxi: 
                    flag=0
                e.pop(0)
        lv2.append(lv)
    #print maxi , lv2

    fl=0
    n=0
    n2=0
    for x in range(len(lv2)-1):
        if lv2[x+1]>lv2[x]:
            n+=1
            fl=1
        else:
            if fl==1:
                start.append(n)
                ends.append(n2)
                maxis.append(lv2[x])
                if  int(e2[n2])-int(s2[n])>overlap_size:
			ol.append([lv2[x], s2[n], e2[n2]])
                else: print "too short"
            n2+=1
            fl=0

    #print start
    #print ends
    #print maxis
    return ol




if __name__=="__main__":
	import doctest
    	#print doctest.testmod()
	o=max_overlap([0,10,20,25,30,50,80],[5,35,40,60,70,90,100], 4)
	print o	
