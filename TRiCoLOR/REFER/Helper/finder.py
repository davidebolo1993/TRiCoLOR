#!/usr/bin/python env

#python 3 standard library

import re
from bisect import bisect_left, bisect_right
from collections import defaultdict, Counter
from operator import itemgetter

#additional libraries

import pyfaidx
import pysam
import editdistance




def nestover(SortedIntervals, string, coords):

	extended=[]

	i=0

	while i < len(SortedIntervals):

		if extended==[]:

			s_,e_=SortedIntervals[i][1], SortedIntervals[i][2]
			string_s,string_e = bisect_left(coords,s_), bisect_right(coords, e_)
			num=string[string_s:string_e].count(SortedIntervals[i][0])
			extended.append((SortedIntervals[i][0], s_, e_, num))
		

		else:

			if extended[-1][2] >= SortedIntervals[i][1] and extended[-1][2] < SortedIntervals[i][2]: #the two intervals overlap in end(1)/start(2)

				if abs(extended[-1][2]-SortedIntervals[i][1]) >= round((extended[-1][2]-extended[-1][1])/2): #big overlap, check if we want to keep both or just one

					#count for the old interval
					c1=extended[-1][2]

					#count for the new interval
					s_,e_=SortedIntervals[i][1], SortedIntervals[i][2]
					ss__, se__=bisect_left(coords,s_), bisect_right(coords, e_)
					m2=SortedIntervals[i][0]
					c2=string[ss__:se__].count(SortedIntervals[i][0])

					#motifs
					m1=extended[-1][0]

					#extend and rank
					new_s,new_e=min(extended[-1][1],SortedIntervals[i][1]), max(extended[-1][2],SortedIntervals[i][2])
					string_s,string_e = bisect_left(coords,new_s), bisect_right(coords, new_e)
					st_=string[string_s:string_e]
					rank1,count1=Markovchain(m1,st_),st_.count(m1)
					rank2,count2=Markovchain(m2,st_),st_.count(m2)

					if rank2 > rank1: #rank2 higher

						if count2 > c2:

							extended.remove(extended[-1])
							extended.append((m2,new_s, new_e, count2))

						else: #extending did not increase the number of repetitions, keep original for second

							extended.remove(extended[-1])
							extended.append((m2,s_,e_, c2))


					elif rank2 < rank1: #rank1 higher

						if count1 > c1:

							extended.remove(extended[-1])
							extended.append((m1,new_s, new_e, count1))

						else: #extending did not increase the number of repetitions, keep original interval for first

							pass


					else: #two ranks are the same

						if m2 in m1*2: #if rep in second interval is the same of rep in first interval but simply rotated keep interval with higher count

							if count2 > count1:

								if count2 > c2:

									extended.remove(extended[-1])
									extended.append((m2,new_s, new_e, count2))

								else:

									extended.remove(extended[-1])
									extended.append((m2,s_,e_, count2))


							else:

								if count1 > c1:

									extended.remove(extended[-1])
									extended.append((m1,new_s,new_e,count1))

								else:

									pass

						else: #different reps with same probability, keep previous and new but separate. They will be merged in VCF but motifs will be kept separated

							extended.append((m2,s_, e_, c2))

				else: #overlapping region is shorter than half the length of previous. Keep both reps, they will be merged in VCF but motifs will be kept separated

					s_,e_=SortedIntervals[i][1], SortedIntervals[i][2]
					string_s,string_e = bisect_left(coords,s_), bisect_right(coords, e_)
					num = string[string_s:string_e].count(SortedIntervals[i][0])
					extended.append((SortedIntervals[i][0], s_, e_, num))
					
			
			elif extended[-1][1] <=  SortedIntervals[i][1] and extended[-1][2] >= SortedIntervals[i][2]:

				pass
			
			elif extended[-1][2] < SortedIntervals[i][1]: #following does not overlap and is not nested

				s_,e_=SortedIntervals[i][1], SortedIntervals[i][2]
				string_s,string_e = bisect_left(coords,s_), bisect_right(coords, e_)
				num = string[string_s:string_e].count(SortedIntervals[i][0])
				extended.append((SortedIntervals[i][0], s_, e_, num))

		i+=1


	return extended



def RepeatsFinder(string,kmer,times, maxkmerlength, overlapping): 

	seen=set()

	if kmer != 0 and times == 0:

		my_regex = r'(.{'  + str(kmer) + r',})\1+' if not overlapping else r'(?=(.{'  + str(kmer) + r',})\1+)'

	elif kmer == 0 and times != 0:

		my_regex = r'(.+?)\1{' + str(times-1) + r',}' if not overlapping else r'(?=(.+?)\1{' + str(times-1) + r',})'

	elif kmer != 0 and times != 0:

		my_regex = r'(.{'  + str(kmer) + r',})\1{' + str(times-1) + r',}' if not overlapping else r'(?=(.{'  + str(kmer) + r',})\1{' + str(times-1) + r',})'

	else:

		my_regex = r'(.+?)\1+' if not overlapping else r'(?=(.+?)\1+)'

	r=re.compile(my_regex)

	for match in r.finditer(string):

		motif=match.group(1)

		if motif not in seen and len(motif) <= maxkmerlength:

			seen.add(motif)

	return seen


def Get_Alignment_Positions(bamfilein): #as the consensus sequence is supposed to generate just one sequence aligned to the reference, secondary alignments are removed
  
	coords=[]
	seq=[]

	bamfile=pysam.AlignmentFile(bamfilein,'rb')

	for read in bamfile.fetch():

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			coords = read.get_reference_positions(full_length=True)
			seq=read.seq

		else:

			bamfile.close()
			return coords,seq

	bamfile.close()

	return coords,seq


def modifier(coordinates): #fast way to remove None and substitute with closest number in list

	
	coordinates=[el+1 if el is not None else el for el in coordinates] #get true coordinates
	start=next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:

			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates



def Get_Alignment_Coordinates(coord_list,repetitions):

	rep_coord_list=[(a,b,c) for a,b,c in zip([repetition[0] for repetition in repetitions],[coord_list[repetition[1]] for repetition in repetitions],[coord_list[repetition[2]] for repetition in repetitions])]

	return rep_coord_list


def look_for_self(rep,sequence): #build another regular expression, looking for the core of the repetition

	for match in re.finditer(rep,sequence):

		yield (match.group(), match.start(), match.start()+len(match.group())*int(len(match.group())/len(match.group()))-1,int(len(match.group())/len(match.group())))


def dfs(adj_list, visited, vertex, result, key):

	visited.add(vertex)
	result[key].append(vertex)

	for neighbor in adj_list[vertex]:

		if neighbor not in visited:

			dfs(adj_list, visited, neighbor, result, key)


def check_edit(string1, string2, allowed): #cython-based way to check for edit distance, faster than check every time for neighbors

	if len(string2) <= allowed: #string 2 is empty or is motif of length leq allowed

		return True

	elif editdistance.eval(string1,string2) <= allowed:

		return True

	return False



def Markovchain(motif,string):

	STATE_LEN = len(motif)

	model = defaultdict(Counter)

	for i in range(len(string) - STATE_LEN):
		state = string[i:i + STATE_LEN]
		next_ = string[i + STATE_LEN :i + STATE_LEN*2]
		model[state][next_] += 1


	motif_probability=(model[motif][motif]/len(string))*len(motif) + len(motif)/len(string) #reward based on motif length

	return motif_probability



def check_ref(string1, string2): #check if ref_string and test_string are different

	return string1 != string2



def corrector(reference, string, repetitions, coordinates, size, allowed): # correct for n-nucleotides insertions/deletions/substitutions inside a repetitive range

	corr_=[]

	coords=modifier(coordinates)

	for reps in repetitions:

		self_=list(look_for_self(reps,string))
		self__= Get_Alignment_Coordinates(coords,self_)

		ranges=[]

		for i in range(len(self_)-1):

			if self__[i+1][1]-self__[i][1] <= len(reps): #coords are subsequent or in inserted region

				if check_edit(reps,string[self_[i][1]+len(reps):self_[i+1][1]], allowed): #reps are subsequent or separated by n-allowed chars or by a string with edit distance n from reps

					ranges.append((self_[i][1],self_[i+1][1])) #accept as putative interval with reps

				else:

					continue

			else: #coordinates are more distant than len reps

				if check_edit(reps,string[self_[i][1]+len(reps):self_[i+1][1]], allowed) and check_ref(reference[self__[i][1]-1:self__[i+1][1]],string[self_[i][1]: self_[i+1][1]+1]):

					ranges.append((self_[i][1],self_[i+1][1])) #allows same motifs separated by a string with edit distance n from the rep

				else:

					continue


		collapsed_ranges= defaultdict(list)
	  
		for x, y in ranges:

			collapsed_ranges[x].append(y)
			collapsed_ranges[y].append(x)

		result = defaultdict(list)
		visited = set()
	  
		for vertex in collapsed_ranges:

			if vertex not in visited:

				dfs(collapsed_ranges, visited, vertex, result, vertex) #collapse ranges that have common bounds


		if len(result) != 0:

			corr_.extend(Get_Alignment_Coordinates(coords,[(reps, interv[0], interv[-1]+(len(reps)-1)) for interv in result.values() if interv[-1]+len(reps)-interv[0] >= size]))

		else:

			continue

	if corr_ == []: #nearly impossible

		return corr_

	s_corr_=sorted(corr_, key=itemgetter(1,2))

	mod_int=nestover(s_corr_, string, coords)

	return mod_int
