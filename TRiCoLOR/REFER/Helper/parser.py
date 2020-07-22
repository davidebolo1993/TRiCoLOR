#!/usr/bin/python env


#python 3 standard library

import os
import glob
import math
import statistics
import subprocess


# additional modules

import pysam


## FUNCTIONS


def sub_none(list_of_coord):


	return [-999999 if v is None else v for v in list_of_coord]


def Bamfile_Analyzer(bamfilein,chromosome,start,end, coverage, out, processor):


	cov=0
	
	headers=[]
	sequences=[]
	lengths=[]

	bamfile=pysam.AlignmentFile(bamfilein,'rb')	

	for read in bamfile.fetch(chromosome,start,end):

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			read_start=read.reference_start
			read_end=read.reference_end

			if read_start <= start and read_end >= end:

				sequence=read.seq
				coord=sub_none(read.get_reference_positions(full_length=True))
				header=read.query_name

				s_,e_=min(coord, key=lambda x:abs(x-start)), min(coord, key=lambda x:abs(x-end)) 
				s_i,e_i = [i for i,e in enumerate(coord) if e == s_][0], [i for i,e in enumerate(coord) if e == e_][0]
				finalseq=sequence[s_i:e_i+1]

				headers.append(header)
				sequences.append(finalseq)
				lengths.append(len(finalseq))


	bamfile.close()

	if len(sequences) >= coverage:

		averagelen=statistics.mean(lengths)
		stdevlen=statistics.stdev(lengths)

		minbound=averagelen-3*stdevlen
		maxbound=averagelen+3*stdevlen

		fasta=''

		for head,seq,leng in zip(headers,sequences,lengths):

			if leng <= minbound or leng >= maxbound:

				continue

			else:

				cov+=1
				
				if cov > 100:
					
					break #do not add other sequences, as it slows down the subsequent consensus computation
					
				else:
				
					fasta+='>' + head + '\n' + seq + '\n'

		if cov >= coverage:

			with open(os.path.abspath(out+'/' + processor + '.unaligned.fa'),'w') as fastaout:

				fastaout.write(fasta.rstrip())

	else:

		cov=len(sequences)

	return cov
