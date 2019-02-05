#!/usr/bin/python env


import sys
import os
import glob
import pysam
import math
from bisect import bisect_left
import subprocess


def Sub_None(list_of_coord):

	return [-999999 if v is None else v for v in list_of_coord] #return a very small number so that it won't be take into account


def Start_Index(list_of_coord,start):

	if start in list_of_coord:

		index=bisect_left(list_of_coord, start)
		return index

	else:

		if any(x < start and x != -999999 for x in list_of_coord): #if there is something lower than start, take closest

			alternative=min(list_of_coord, key=lambda x:abs(x-start))
			index= bisect_left(list_of_coord, alternative)
			return index

		else:

			return #no overlapping start


def End_Index(list_of_coord,end):

	if end in list_of_coord:

		index=bisect_left(list_of_coord, end)
		return index

	else:

		if any(x > end for x in list_of_coord): #check if any coord is higher in read than end, else return None (sequence has not an end in interval)

			alternative = min(list_of_coord, key=lambda x:abs(x-end))
			index=bisect_left(list_of_coord, alternative)

		else:

			return


def Basic_Infos(pysam_AlignedSegment):

	identifier=pysam_AlignedSegment.query_name	
	seq=pysam_AlignedSegment.seq #give the correct result even for reverse strand, no need to reverse-complement
	
	return (identifier,seq)



def Get_Coords(pysam_AlignedSegment,start,end):


	a_coord=pysam_AlignedSegment.get_reference_positions(full_length=True) #list that indicates the reference position on which each base alignes (None included for insertion/soft clipped bases)
	identifier=pysam_AlignedSegment.query_name
	
	a_coord=Sub_None(a_coord)	

	ind_start=Start_Index(a_coord,start)
	ind_end=End_Index(a_coord,end)

	if ind_start is None or ind_end is None:

		return

	elif ind_end - ind_start <= 10: #set a cut-off. If start and end indexes are approximatley the same, it means that wanted coordinates are not in the sequence (closest are taken)

		return

	else:
		
		return (identifier,ind_start,ind_end)


def split_equal(value, parts):

	value = float(value)
	
	return int(value/parts), len([i*value/parts for i in range(1,parts+1)])


def SizeChecker(start, end, treshold=2): #adjust size of interval so that we don't have sequences too short and they can me mapped with more accuracy

	if end-start <= treshold*1000:

		return start, end, end-start

	else:

		size, len_=split_equal(end-start, math.ceil((end-start)/(treshold*1000)))

		return start, start+(size*len_), size 
	


def check_coverage(pysam_AlignmentFile, chromosome, start, end, coverage):

	counter = 0

	for reads in pysam_AlignmentFile.fetch(chromosome, start, end):

		counter +=1

	if counter >= coverage:

		return True

	else:

		return False


def Bamfile_Analyzer(bamfilein,chromosome,start,end, coverage): #double sliding window function to retrive informations

	bamfile=pysam.AlignmentFile(bamfilein,'rb')
	start,end,size=SizeChecker(start,end) #split the period in sub-period of max length 2000 bp

	i=0

	pointer=start+1 #slicing window

	inf=[] #list for appending results
	coord=[] #list for appending coordinates

	while end-1 >= start:

		inf.append([])
		coord.append([])

		if check_coverage(bamfile, chromosome, start, pointer, coverage):

			for read in bamfile.fetch(chromosome,start,pointer):

				if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

					inf[i].append(Basic_Infos(read))

					if end > start+size:

						coord[i].append(Get_Coords(read,start,start+size))

					else:

						coord[i].append(Get_Coords(read,start,end))

		if end-start > size:

			start += size
			pointer += size

		elif end-start <= size and start < end-1:

			start = end-1
			pointer = end

		else: 

			break

		i+=1 #update counter in order to store results progressively

	bamfile.close()

	return inf,coord[:len(coord)-1] 


def Reorder(l1,l2): #reorder and filtering

	l2_2=list(filter(None.__ne__,l2)) #filtering None for reads that do not have an end in interval
	l2_el=[el[0] for el in l2_2]
	mapping = dict((x[0], x[1:]) for x in l1)
	l1[:] = [(x,) + mapping[x] for x in l2_el]

	return l1,l2_2


def InCommon(info_list,coord_list):

	filtered_seq=[]
	filtered_coord=[]

	for i in range(len(info_list)-1):

		fip=info_list[i]
		sip=info_list[i+1]


		if fip != [] and sip != []:

			id_fip=[a for a,b in fip]
			id_sip=[a for a,b in sip]
			int_id=list(set(id_fip).intersection(id_sip))
			com_seq=[(a,b) for a,b in fip if a in int_id]
			l1,l2=Reorder(com_seq,coord_list[i])
			filtered_seq.append(l1)
			filtered_coord.append(l2)

		else:

			continue

	return [el for el in filtered_seq if el != []],[el for el in filtered_coord if el !=[]]


def Fasta_Generator(filtered_seq,filtered_coord,out):

	for i in range(len(filtered_seq)):

		dirname=os.path.abspath(out)

		if not os.path.exists(dirname):

			os.makedirs(dirname)

		fname=os.path.abspath(dirname+'/Fasta'+str(i+1)+'.unaligned.fa')

		f=open(fname,'a')

		for j in range(len(filtered_seq[i])):

			fasta_in_window='>'+filtered_seq[i][j][0]+'\n'+filtered_seq[i][j][1][filtered_coord[i][j][1]:filtered_coord[i][j][2]+1]+'\n' #get sequence from start to end (end included)
			f.write(fasta_in_window)

		f.close() 



def MA(out,mmi_ref): #MSA function that uses Alfred consensus as it is much faster than a MSA python implementation.


	fafile=glob.glob(os.path.abspath(out)+'/*.unaligned.fa')

	outfa=[]
	outall=[]

	for fa in fafile:

		outfa.append(fa.replace('unaligned.fa','consensus.fa.gz'))
		outall.append(fa.replace('unaligned.fa','alignment.txt.gz'))

	for fastain,consensusout,alignmentout in zip(fafile,outfa,outall):

		subprocess.call(['/home/bolognin/Downloads/alfred/src/alfred', 'consensus', '-a', alignmentout, '-c', consensusout, '-f', 'fasta', '-t', 'ont', fastain]) #slowest part

		with open(consensusout.replace('.fa.gz','.sam'), 'w') as f:

			subprocess.call(['minimap2', '-ax', 'map-ont', mmi_ref, consensusout], stdout=f)

		with open(consensusout.replace('.fa.gz','.bam'), 'w') as f:

				subprocess.call(['htsbox', 'samview', '-bS', consensusout.replace('.fa.gz','.sam')], stdout=f)

		with open(consensusout.replace('.fa.gz','.srt.bam'), 'w') as f:

			subprocess.call(['samtools', 'sort', consensusout.replace('.fa.gz','.bam')], stdout=f)

		subprocess.call(['samtools', 'index', consensusout.replace('.fa.gz','.srt.bam')])
		
	types = ['/*.bed', '/*srt.bam', '/*srt.bam.bai']
	tokeep = []

	for files in types:

		tokeep.extend(glob.glob(os.path.abspath(out)+files))


	for dirpath,_,filenames in os.walk(os.path.abspath(out)):

		for f in filenames:

			if os.path.abspath(os.path.join(dirpath, f)) not in tokeep:

				os.remove(os.path.abspath(os.path.join(dirpath, f)))
