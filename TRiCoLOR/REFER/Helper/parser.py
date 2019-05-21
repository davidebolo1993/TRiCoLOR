#!/usr/bin/python env


#python 3 standard library

import os
import glob
import math
import subprocess


# additional modules

import pysam


def sub_none(list_of_coord):


	return [-999999 if v is None else v for v in list_of_coord] #return a very small number so that it won't be take into account


def check_coverage(pysam_AlignmentFile, chromosome, start, end, coverage): #check if we have at least the wanted coverage for the interval


	counter = 0

	for read in pysam_AlignmentFile.fetch(chromosome, start, end):
		
		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			if read.reference_start <= start and read.reference_end >= end: #make sure that the same read cover the entire region

				counter +=1

	if counter >= coverage:

		return True

	else:

		return False


def fastagen(header,sequence,iteration,out):


	filename=os.path.abspath(out+'/Fasta'+str(iteration+1)+'.unaligned.fa')

	with open(os.path.abspath(out+'/Fasta'+str(iteration+1)+'.unaligned.fa'), 'a') as fastaout:

		fastaout.write('>' + header + '\n' + sequence + '\n')


def split_equal(value, parts):


	value = float(value)
	
	return int(value/parts), len([i*value/parts for i in range(1,parts+1)])


def sizechecker(start, end): #adjust size of interval so that we don't have sequences too short and they can me mapped with more accuracy


	if end-start <= 2000:

		return start, end, end-start

	else:

		size, len_=split_equal(end-start, math.ceil((end-start)/2000))

		return start, start+(size*len_), size 


def Bamfile_Analyzer(bamfilein,chromosome,start,end, coverage, fastaout):


	bamfile=pysam.AlignmentFile(bamfilein,'rb')	
	start,end,size=sizechecker(start,end) #split the period in sub-period of equal length
	next_=start+size
	iteration=0
	final=False

	while not final:

		if check_coverage(bamfile, chromosome, start, next_, coverage): # coverage check for the region

			for read in bamfile.fetch(chromosome,start,next_):

				if not read.is_unmapped and not read.is_secondary and not read.is_supplementary: #skip everything but primary, as in coverage check

					read_start=read.reference_start
					read_end=read.reference_end

					if read_start <= start and read_end >= read_end: #make sure that the same read cover the entire region

						sequence=read.seq
						coord=sub_none(read.get_reference_positions(full_length=True))
						header=read.query_name

						s_,e_=min(coord, key=lambda x:abs(x-start)), min(coord, key=lambda x:abs(x-(start+size))) #must exist
						s_i,e_i = [i for i,e in enumerate(coord) if e == s_][0], [i for i,e in enumerate(coord) if e == e_][0] #must exist

						fastagen(header,sequence[s_i:e_i+1],iteration,fastaout)

		iteration+=1

		if end-next_ >= size:

			start += size
			next_ += size

		else:

			final=True
	
	bamfile.close()

