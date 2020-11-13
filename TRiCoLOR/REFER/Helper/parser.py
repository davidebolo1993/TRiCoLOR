#!/usr/bin/python env


#python 3 standard library

import os
import math
import statistics

# additional modules

import pysam
import numpy as np


def sub_none(coordinates):

	'''
	Substitute None (soft-clipped/inserted) coordinates with a negative value
	'''

	return [-9999999 if v is None else v for v in coordinates]


def find_nearest(array, value):

	'''
	Find closest match in array, given value
	'''

	return (np.abs(array - value)).argmin()


def Bamfile_Analyzer(bamfilein,chromosome,start,end,c,out,processor):

	'''
	Parse haplotype-resolved BAM and store sequences in FASTA for consensus calculation
	'''

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

				sequence=read.query_sequence #changed seq with query sequence
				coord=np.asarray(sub_none(read.get_reference_positions(full_length=True)))
				header=read.query_name

				#s_,e_=min(coord, key=lambda x:abs(x-start)), min(coord, key=lambda x:abs(x-end)) 
				#s_i,e_i = [i for i,e in enumerate(coord) if e == s_][0], [i for i,e in enumerate(coord) if e == e_][0]
				#finalseq=sequence[s_i:e_i+1]
				si,ei=find_nearest(coord,start),find_nearest(coord,end)
				finalseq=sequence[si:ei]

				headers.append(header)
				sequences.append(finalseq)
				lengths.append(len(finalseq))


	bamfile.close()

	if len(sequences) >= c.coverage:

		averagelen=statistics.mean(lengths)
		stdevlen=statistics.stdev(lengths)

		minbound=averagelen-3*stdevlen
		maxbound=averagelen+3*stdevlen

		fasta=''

		for head,seq,leng in zip(headers,sequences,lengths):

			if leng <= minbound or leng >= maxbound:

				continue #exclude outliers

			else:

				cov+=1
				
				if cov > 50: #arbitrary treshold to avoid problems with coverage too high
					
					break #do not add other sequences, as it slows down the subsequent consensus computation
					
				else:
				
					fasta+='>' + head + '\n' + seq + '\n'

		if cov >= c.coverage:

			with open(os.path.abspath(out+'/' + processor + '.unaligned.fa'),'w') as fastaout:

				fastaout.write(fasta)

	else:

		cov=len(sequences)

	return cov


def Bamfile_Analyzer_Single(bamfilein,chromosome,start,end,c,processor):

	'''
	Parse haplotype-tagged BAM and store sequences in FASTA for consensus calculation
	'''

	cov1=0
	cov2=0
	headers1=[]
	sequences1=[]
	lengths1=[]
	headers2=[]
	sequences2=[]
	lengths2=[]

	bamfile=pysam.AlignmentFile(bamfilein,'rb')	

	for read in bamfile.fetch(chromosome,start,end):

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			if not read.has_tag('HP'):

				continue

			else:


				read_start=read.reference_start
				read_end=read.reference_end

				if read_start <= start and read_end >= end:

					sequence=read.query_sequence #changed seq with query sequence
					coord=np.asarray(sub_none(read.get_reference_positions(full_length=True)))
					header=read.query_name

					#s_,e_=min(coord, key=lambda x:abs(x-start)), min(coord, key=lambda x:abs(x-end)) 
					#s_i,e_i = [i for i,e in enumerate(coord) if e == s_][0], [i for i,e in enumerate(coord) if e == e_][0]
					#finalseq=sequence[s_i:e_i+1]
					si,ei=find_nearest(coord,start),find_nearest(coord,end)
					finalseq=sequence[si:ei]

					if read.get_tag('HP') == 1: 

						headers1.append(header)
						sequences1.append(finalseq)
						lengths1.append(len(finalseq))

					else:

						headers2.append(header)
						sequences2.append(finalseq)
						lengths2.append(len(finalseq))

	bamfile.close()

	#hap1

	if len(sequences1) >= c.coverage:

		averagelen=statistics.mean(lengths1)
		stdevlen=statistics.stdev(lengths1)

		minbound=averagelen-3*stdevlen
		maxbound=averagelen+3*stdevlen

		fasta=''

		for head,seq,leng in zip(headers1,sequences1,lengths1):

			if leng <= minbound or leng >= maxbound:

				continue #exclude outliers

			else:

				cov1+=1
				
				if cov1 > 50: #arbitrary treshold to avoid problems with coverage too high
					
					break #do not add other sequences, as it slows down the subsequent consensus computation
					
				else:
				
					fasta+='>' + head + '\n' + seq + '\n'

		if cov1 >= c.coverage:

			with open(os.path.abspath(c.OUT+'/haplotype1/' + processor + '.unaligned.fa'),'w') as fastaout:

				fastaout.write(fasta)

	else:

		cov1=len(sequences1)


	#hap2

	if len(sequences2) >= c.coverage:

		averagelen=statistics.mean(lengths2)
		stdevlen=statistics.stdev(lengths2)

		minbound=averagelen-3*stdevlen
		maxbound=averagelen+3*stdevlen

		fasta=''

		for head,seq,leng in zip(headers2,sequences2,lengths2):

			if leng <= minbound or leng >= maxbound:

				continue #exclude outliers

			else:

				cov2+=1
				
				if cov2 > 50: #arbitrary treshold to avoid problems with coverage too high
					
					break #do not add other sequences, as it slows down the subsequent consensus computation
					
				else:
				
					fasta+='>' + head + '\n' + seq + '\n'

		if cov2 >= c.coverage:

			with open(os.path.abspath(c.OUT+'/haplotype2/' + processor + '.unaligned.fa'),'w') as fastaout:

				fastaout.write(fasta)

	else:

		cov2=len(sequences2)


	return cov1,cov2