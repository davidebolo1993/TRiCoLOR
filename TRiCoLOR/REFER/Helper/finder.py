#!/usr/bin/python env


#python 3 standard library

import re
from collections import defaultdict, Counter
from operator import itemgetter


#additional modules

import pyfaidx
import pysam
import editdistance


## FUNCTIONS


def SolveNestedR(SortedIntervals):


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


def ReferenceFilter(reference_reps,wanted,size,start): 


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

            new_reps=[(reps, start+val[0], start+val[-1]+len(reps)-1, len(val)) for val in list(result.values()) if val[-1]+len(reps)-val[0] >= size]
            corr_.extend(new_reps)

    if corr_ == []:

        return corr_

    else:

        s_corr_=sorted(corr_, key=itemgetter(1,2))
        mod_int=SolveNestedR(s_corr_)

        return mod_int


def SolveNestedH(SortedIntervals, string, size):


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


def possible_rotations(word):

    p = []

    for i in range(len(word)):

        p.append(word[i:]+word[:i])
        
    return p


def RegexBuilder(kmer,times,overlapping,strictmotif,strictimes):

    if (kmer == 0 or kmer==1) and (times == 0 or times==1):

        if not strictmotif:

            my_regex = r'(.+?)\1+' if not overlapping else r'(?=(.+?)\1+)'

        else:

            if kmer == 1:

                my_regex = r'(.{1})\1+' if not overlapping else r'(?=(.{1})\1+)'

    elif (kmer == 0 or kmer == 1) and (times !=0 and times!=1):

        if not strictmotif:

            if not strictimes:

                my_regex = r'(.+?)\1{' + str(times-1) + r',}' if not overlapping else r'(?=(.+?)\1{' + str(times-1) + r',})'

            else:

                my_regex = r'(.+?)\1{' + str(times-1) + r'}' if not overlapping else r'(?=(.+?)\1{' + str(times-1) + r'})'

        else:

            if kmer==1:

                if not strictimes:

                    my_regex = r'(.{1})\1{' + str(times-1) + r',}' if not overlapping else r'(?=(.{1})\1{' + str(times-1) + r',})'

                else:

                    my_regex = r'(.{1})\1{' + str(times-1) + r'}' if not overlapping else r'(?=(.{1})\1{' + str(times-1) + r'})'                    

    elif (kmer != 0 and kmer !=1) and (times==0 or times ==1):

        if not strictmotif:

            my_regex = r'(.{'  + str(kmer) + r',})\1+' if not overlapping else r'(?=(.{'  + str(kmer) + r',})\1+)'

        else:

            my_regex = r'(.{'  + str(kmer) + r'})\1+' if not overlapping else r'(?=(.{'  + str(kmer) + r'})\1+)'

    elif (kmer != 0 and kmer !=1) and (times !=0 and times!=1):

        if not strictmotif and not strictimes:

            my_regex = r'(.{'  + str(kmer) + r',})\1{' + str(times-1) + r',}' if not overlapping else r'(?=(.{'  + str(kmer) + r',})\1{' + str(times-1) + r',})'

        elif strictmotif and not strictimes:

            my_regex = r'(.{'  + str(kmer) + r'})\1{' + str(times-1) + r',}' if not overlapping else r'(?=(.{'  + str(kmer) + r'})\1{' + str(times-1) + r',})'

        elif not strictmotif and strictimes:

            my_regex = r'(.{'  + str(kmer) + r',})\1{' + str(times-1) + r'}' if not overlapping else r'(?=(.{'  + str(kmer) + r',})\1{' + str(times-1) + r'})'

        else:

            my_regex = r'(.{'  + str(kmer) + r'})\1{' + str(times-1) + r'}' if not overlapping else r'(?=(.{'  + str(kmer) + r'})\1{' + str(times-1) + r'})'

    return my_regex


def RepeatsFinder(string, regex, maxkmerlength): 

    seen=set()

    r=re.compile(regex)

    for match in r.finditer(string):

        motif=match.group(1)

        if len(motif) <= maxkmerlength:

            seen.add(motif)

    return seen


def Get_Alignment_Positions(bamfilein):

  
    coords=[]
    seq=[]

    bamfile=pysam.AlignmentFile(bamfilein,'rb')

    for read in bamfile.fetch():

        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

            coords = read.get_reference_positions(full_length=True)
            seq=read.seq

    bamfile.close()

    return coords,seq


def modifier(coordinates):

    
    coordinates=[el+1 if el is not None else el for el in coordinates]
    start=next(ele for ele in coordinates if ele is not None)

    for ind, ele in enumerate(coordinates):
        
        if ele is None:

            coordinates[ind] = start
        
        else:

            start = ele

    return coordinates


def Get_Alignment_Coordinates(coord_list,repetitions):


    rep_coord_list=[(a,b,c,d) for a,b,c,d in zip([repetition[0] for repetition in repetitions],[coord_list[repetition[1]] for repetition in repetitions],[coord_list[repetition[2]] for repetition in repetitions],[repetition[3] for repetition in repetitions])]

    return rep_coord_list


def look_for_self(rep,sequence):


    r=re.compile(rep)

    for match in r.finditer(sequence):

        yield (match.group(), match.start(), match.start()+len(match.group())*int(len(match.group())/len(match.group()))-1,int(len(match.group())/len(match.group())))


def dfs(adj_list, visited, vertex, result, key):


    visited.add(vertex)
    result[key].append(vertex)

    for neighbor in adj_list[vertex]:

        if neighbor not in visited:

            dfs(adj_list, visited, neighbor, result, key)


def check_edit(string1, string2, allowed):


    if len(string2) <= allowed:

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

    motif_probability=(model[motif][motif]/len(string))*len(motif) + len(motif)/len(string)

    return motif_probability, model[motif][motif]+1


def check_ref(string1, string2): 


    return string1 != string2


def corrector(reference, string, repetitions, coordinates, size, allowed):


    corr_=[]
    coords=modifier(coordinates)

    for reps in repetitions:

        self_=list(look_for_self(reps,string))
        self__= Get_Alignment_Coordinates(coords,self_)

        ranges=[]

        for i in range(len(self_)-1):

            if self_[i+1][1]-self_[i][1] == len(reps): 

                ranges.append((self_[i][1],self_[i+1][1]))

            else:

                if len(reps) > 1:

                    if check_edit(reps,string[self_[i][1]+len(reps):self_[i+1][1]], allowed) and check_ref(reference[self__[i][1]-1:self__[i+1][1]],string[self_[i][1]: self_[i+1][1]+1]):

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

            corr_.extend([(reps, interv[0], interv[-1]+(len(reps)-1)) for interv in result.values() if interv[-1]+len(reps)-interv[0] >= size])

        else:

            continue

    if corr_ == []:

        return corr_

    s_corr_=sorted(corr_, key=itemgetter(1,2))

    mod_int=SolveNestedH(s_corr_, string,size)

    return Get_Alignment_Coordinates(coords,mod_int)
