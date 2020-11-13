#!/usr/bin/python3 env


#python 3 standard library

import re
from collections import defaultdict, Counter
from operator import itemgetter

#additional modules

import editdistance


def possible_rotations(word):

	'''
	Build possible rotations of a motif
	'''

	p = []

	for i in range(len(word)):

		p.append(word[i:]+word[:i])
		
	return p


def RegexBuilder(c):

	'''
	Build a RegEx to find repeated patterns in a consensus string
	'''

	if (c.motif == 0 or c.motif==1) and (c.times == 0 or c.times==1):

		if not c.precisemotif:

			my_regex = r'(.+?)\1+' if not c.overlapping else r'(?=(.+?)\1+)'

		else:

			if c.motif == 1:

				my_regex = r'(.{1})\1+' if not c.overlapping else r'(?=(.{1})\1+)'

	elif (c.motif == 0 or c.motif == 1) and (c.times !=0 and c.times!=1):

		if not c.precisemotif:

			if not c.precisetimes:

				my_regex = r'(.+?)\1{' + str(c.times-1) + r',}' if not c.overlapping else r'(?=(.+?)\1{' + str(c.times-1) + r',})'

			else:

				my_regex = r'(.+?)\1{' + str(c.times-1) + r'}' if not c.overlapping else r'(?=(.+?)\1{' + str(c.times-1) + r'})'

		else:

			if c.motif==1:

				if not c.precisetimes:

					my_regex = r'(.{1})\1{' + str(c.times-1) + r',}' if not c.overlapping else r'(?=(.{1})\1{' + str(c.times-1) + r',})'

				else:

					my_regex = r'(.{1})\1{' + str(c.times-1) + r'}' if not c.overlapping else r'(?=(.{1})\1{' + str(c.times-1) + r'})'                    

	elif (c.motif != 0 and c.motif !=1) and (c.times==0 or c.times ==1):

		if not c.precisemotif:

			my_regex = r'(.{'  + str(c.motif) + r',})\1+' if not c.overlapping else r'(?=(.{'  + str(c.motif) + r',})\1+)'

		else:

			my_regex = r'(.{'  + str(c.motif) + r'})\1+' if not c.overlapping else r'(?=(.{'  + str(c.motif) + r'})\1+)'

	elif (c.motif != 0 and c.motif !=1) and (c.times !=0 and c.times!=1):

		if not c.precisemotif and not c.precisetimes:

			my_regex = r'(.{'  + str(c.motif) + r',})\1{' + str(c.times-1) + r',}' if not c.overlapping else r'(?=(.{'  + str(c.motif) + r',})\1{' + str(c.times-1) + r',})'

		elif c.precisemotif and not c.precisetimes:

			my_regex = r'(.{'  + str(c.motif) + r'})\1{' + str(c.times-1) + r',}' if not c.overlapping else r'(?=(.{'  + str(c.motif) + r'})\1{' + str(c.times-1) + r',})'

		elif not c.precisemotif and c.precisetimes:

			my_regex = r'(.{'  + str(c.motif) + r',})\1{' + str(c.times-1) + r'}' if not c.overlapping else r'(?=(.{'  + str(c.motif) + r',})\1{' + str(c.times-1) + r'})'

		else:

			my_regex = r'(.{'  + str(c.motif) + r'})\1{' + str(c.times-1) + r'}' if not c.overlapping else r'(?=(.{'  + str(c.motif) + r'})\1{' + str(c.times-1) + r'})'

	return my_regex


def RepeatsFinder(string,c):

	'''
	Find repeated motifs in consensus string, having a RegEx built
	'''

	seen=set()

	r=re.compile(c.regex)

	for match in r.finditer(string):

		motif=match.group(1)

		if len(motif) <= c.maxmotif:

			seen.add(motif)

	return seen


def SolveNestedR(SortedIntervals):

	'''
	Solve nested repeats in reference sequence (no errors)
	'''

	extended=[]

	i=0

	while i < len(SortedIntervals):

		if extended==[]:

			extended.append(SortedIntervals[i])

		else:

			if extended[-1][2] >= SortedIntervals[i][1]:

				if SortedIntervals[i][2] - SortedIntervals[i][1] > extended[-1][2] - extended[-1][1]:

					extended.remove(extended[-1])
					extended.append(SortedIntervals[i])
			else:

				extended.append(SortedIntervals[i])
		
		i+=1
	
	return extended


def dfs(adj_list, visited, vertex, result, key):

	'''
	Collapse overlapping edges recursively
	'''

	visited.add(vertex)
	result[key].append(vertex)

	for neighbor in adj_list[vertex]:

		if neighbor not in visited:

			dfs(adj_list, visited, neighbor, result, key)


def ReferenceFilter(chromosome,reference_reps,wanted,c,start):

	'''
	Identify repeated substring/s
	'''

	corr_=[]

	for reps in reference_reps:
	
		self_=list(look_for_self(reps,wanted))
		ranges=[]

		for i in range(len(self_)-1):

			if self_[i+1][1]-self_[i][1] == len(reps):

				ranges.append((self_[i][1],self_[i+1][1]))

		collapsed_ranges= defaultdict(list)
		
		for x, y in ranges:

			collapsed_ranges[x].append(y)
			collapsed_ranges[y].append(x)

		result = defaultdict(list)
		visited = set()
		
		for vertex in collapsed_ranges:

			if vertex not in visited:

				dfs(collapsed_ranges, visited, vertex, result, vertex)

		if len(result) !=0:

			new_reps=[(reps, start+val[0], start+val[-1]+len(reps)-1, len(val)) for val in list(result.values()) if val[-1]+len(reps)-val[0] >= c.size]
			corr_.extend(new_reps)

	if corr_ == []:

		return corr_

	else:

		s_corr_=sorted(corr_, key=itemgetter(1,2))
		mod_int=SolveNestedR(s_corr_)

		return [(chromosome,a,b,c,d) for a,b,c,d in zip([x[0] for x in mod_int],[x[1]-1 for x in mod_int],[x[2]-1 for x in mod_int],[x[3] for x in mod_int])] #convert to 0-based


def look_for_self(motif,sequence):

	'''
	Find all occurence of a motif in string

	'''

	r=re.compile(motif)

	for match in r.finditer(sequence):

		yield (match.group(), match.start(), match.start()+len(match.group())*int(len(match.group())/len(match.group()))-1,int(len(match.group())/len(match.group())))


def check_edit(string1, string2, allowed):

	'''
	Check if the edit-distance between motif and substring is accepted
	'''


	if len(string2) <= allowed:

		return True

	elif editdistance.eval(string1,string2) <= allowed:

		return True

	return False


def Markovchain(motif,string):

	'''
	If 2 motifs overlap, choose the most-likely based on N-gram occurences (Markov chain)
	'''

	state_len = len(motif)
	model = defaultdict(Counter)

	for i in range(len(string) - state_len):
		state = string[i:i + state_len]
		next_ = string[i + state_len :i + state_len*2]
		model[state][next_] += 1

	motif_probability=((model[motif][motif])/len(string))*state_len + state_len/len(string)

	return motif_probability, string.count(motif)


def SolveNestedH(SortedIntervals, string, size):

	'''
	If 2 motifs overlap, choose the most likely
	'''

	extended=[]
	
	i=0

	while i < len(SortedIntervals):

		if extended==[]:

			s_,e_=SortedIntervals[i][1],SortedIntervals[i][2]
			num=string[s_:e_+1].count(SortedIntervals[i][0])
			extended.append((SortedIntervals[i][0], s_, e_, num))
		
		else:

			if extended[-1][2] >= SortedIntervals[i][1] and extended[-1][2] < SortedIntervals[i][2]:

				if abs(extended[-1][2]-SortedIntervals[i][1]) >= round((extended[-1][2]-extended[-1][1])/2):

					c1=string[extended[-1][1]:extended[-1][2]+1].count(extended[-1][0])
					c2=string[SortedIntervals[i][1]:SortedIntervals[i][2]+1].count(SortedIntervals[i][0])
					new_s,new_e=min(extended[-1][1],SortedIntervals[i][1]), max(extended[-1][2],SortedIntervals[i][2])
					st_=string[new_s:new_e+1]
					m1=extended[-1][0]
					m2=SortedIntervals[i][0]
					rank1,count1=Markovchain(m1,st_) #?  OTHER WAY TO DECIDE WHICH REPETITIVE MOTIF IS MORE SIGNIFICANT. NOT PRIORITY. IT WORKS WELL
					rank2,count2=Markovchain(m2,st_)

					if rank2 > rank1:

						if count2 > c2:

							extended.remove(extended[-1])
							extended.append((m2,new_s, new_e, count2))

						else:

							extended.remove(extended[-1])
							extended.append((m2,SortedIntervals[i][1],SortedIntervals[i][2], c2))

					elif rank2 < rank1:

						if count1 > c1:

							extended.remove(extended[-1])
							extended.append((m1,new_s, new_e, count1))

					else: 

						if m2 not in possible_rotations(m1):

							if count2 > count1: 

								extended.append((m2,new_s, new_e, count2))

				else:
					
					new_s,new_e=extended[-1][2]+1,SortedIntervals[i][2]
					st_=string[new_s:new_e+1]

					if len(st_) >=size:

						extended.append((SortedIntervals[i][0],new_s, new_e, st_.count(SortedIntervals[i][0])))
						
			elif extended[-1][2] < SortedIntervals[i][1]:

				s_,e_=SortedIntervals[i][1],SortedIntervals[i][2]
				num=string[s_:e_+1].count(SortedIntervals[i][0])
				extended.append((SortedIntervals[i][0], s_, e_, num))

		i+=1

	return extended


def corrector(chromosome,c, string, repetitions, coordinates):

	'''
	Solve imperfect repetitions. Only perfect occurences are counted but fuzzy occurences are merged into the repeated string
	'''

	corrected=[]

	for reps in repetitions:

		self_=list(look_for_self(reps,string)) #indices in consensus
		self__=[(a,b,c,d) for a,b,c,d in zip([x[0] for x in self_],[coordinates[x[1]] for x in self_],[coordinates[x[2]] for x in self_],[x[3] for x in self_])] #corresponding BAM coordinates

		ranges=[]

		for i in range(len(self_)-1):

			if self_[i+1][1]-self_[i][1] == len(reps): #subsequent motif

				ranges.append((self_[i][1],self_[i+1][1]))

			else:

				if len(reps) > 1: 

					if check_edit(reps,string[self_[i][1]+len(reps):self_[i+1][1]], c.editdistance): #if interruption of perfect reps below treshold

						ranges.append((self_[i][1],self_[i+1][1]))

		collapsed_ranges= defaultdict(list)
	  
		for x, y in ranges:

			collapsed_ranges[x].append(y)
			collapsed_ranges[y].append(x)

		result = defaultdict(list)
		visited = set()

		for vertex in collapsed_ranges:

			if vertex not in visited:

				dfs(collapsed_ranges, visited, vertex, result, vertex)

		if len(result) != 0:

			corrected.extend([(reps, interv[0], interv[-1]+(len(reps)-1)) for interv in result.values() if interv[-1]+len(reps)-interv[0] >= c.size])

	corrected_srtd=sorted(corrected, key=itemgetter(1,2))
	corrected_srtd_fltrd=SolveNestedH(corrected_srtd,string,c.size)

	return [(chromosome,a,b,c,d) for a,b,c,d in zip([x[0] for x in corrected_srtd_fltrd],[coordinates[x[1]] for x in corrected_srtd_fltrd],[coordinates[x[2]] for x in corrected_srtd_fltrd],[x[3] for x in corrected_srtd_fltrd])] #corresponding BAM coordinates
