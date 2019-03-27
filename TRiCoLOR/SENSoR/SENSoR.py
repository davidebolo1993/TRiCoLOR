#!/usr/bin/python env

#python 3 standard library

import os
import sys
import math
import itertools
from operator import itemgetter
from multiprocessing import Process
from shutil import which
import timeit
import logging
import subprocess


# additional libraries

import pysam
import numpy as np



def run(parser, args):


	if not os.path.exists(os.path.abspath(args.output)):

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('It was not possible to create the results folder. Specify a path for which you have write permissions')
			sys.exit(1)

	else: #path already exists

		if not os.access(os.path.dirname(os.path.abspath(args.output)),os.W_OK): #path exists but no write permissions on that folder

			print('You do not have write permissions on the directory in which results will be stored. Specify a folder for which you have write permissions')
			sys.exit(1)

	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR_SENSoR.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


	external_tools=['samtools', 'bedops']

	for tools in external_tools:

		if which(tools) is None:

			logging.error(tools + ' was not found as an executable command. Install ' + tools + ' and re-run TRiCoLOR SENSoR')
			sys.exit(1)

	#check inputs validity

	try:

		subprocess.check_call(['samtools','quickcheck', os.path.abspath(args.bamfile1)],stderr=open(os.devnull, 'wb'))

	except:

		logging.error('.bam file 1 does not exist, is not readable or is not a valid .bam file')
		sys.exit(1)


	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bamfile2)],stderr=open(os.devnull, 'wb'))
		

	except:

		logging.error('.bam file 2 does not exist, is not readable or is not a valid .bam file')
		sys.exit(1)


	if not os.path.exists(os.path.abspath(args.bamfile1 + '.bai')):

		logging.info('Creating index for .bam file 1, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile1)])

		except:

			logging.error('.bam file 1 could not be indexed')
			sys.exit(1)


	if not os.path.exists(os.path.abspath(args.bamfile2 + '.bai')):

		logging.info('Creating index for .bam file 2, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile2)])

		except:

			logging.error('.bam file 2 could not be indexed')
			sys.exit(1)

	start=timeit.default_timer()
	logging.info('Analysis starts now')


	try:

		runInParallel(BScanner, (args.bamfile1, os.path.abspath(args.output + '/' + args.label + '_hap1.bed'), args.scansize, args.entropytreshold, args.calltreshold, args.lengthtreshold),(args.bamfile2, os.path.abspath(args.output + '/' + args.label + '_hap2.bed'), args.scansize, args.entropytreshold, args.calltreshold,args.lengthtreshold))

	except:

		logging.exception('Something went wrong during the scanning process. Log is attached below')
		sys.exit(1)


	end=timeit.default_timer()
	elapsed=(end-start)/60

	logging.info('.bam files scanned in ' + str(elapsed) + ' minutes')


	with open(os.path.abspath(args.output + '/' + args.label + '_hap1.srt.bed'), 'w') as srtbed1:
	
		subprocess.call(['sort-bed', os.path.abspath(args.output + '/' + args.label + '_hap1.bed')],stdout=srtbed1)


	with open(os.path.abspath(args.output + '/' + args.label + '_hap2.srt.bed'), 'w') as srtbed2:
	
		subprocess.call(['sort-bed', os.path.abspath(args.output + '/' + args.label + '_hap2.bed')],stdout=srtbed2)


	os.remove(os.path.abspath(args.output + '/' + args.label + '_hap1.bed')) #remove unsorted
	os.remove(os.path.abspath(args.output + '/' + args.label + '_hap2.bed')) #remove unsorted


	with open(os.path.abspath(args.output + '/' + args.label + '.merged.bed'), 'w') as bedout:

		subprocess.call(['bedops', '-m', os.path.abspath(args.output + '/' + args.label + '_hap1.srt.bed'), os.path.abspath(args.output + '/' + args.label + '_hap2.srt.bed')], stdout=bedout)


	os.remove(os.path.abspath(args.output + '/' + args.label + '_hap1.srt.bed')) #remove un-merged
	os.remove(os.path.abspath(args.output + '/' + args.label + '_hap2.srt.bed')) #remove un-merged


	if args.genometype is None:

		logging.info('.bed file ready')
		logging.info('Done')

	else:

		data_path=os.path.abspath(os.path.dirname(__file__) + '/data') #path to data for filtering


		with open(os.path.abspath(args.output + '/' + args.label + '.filtered.merged.bed'), 'w') as bedout:

			subprocess.call(['bedops', '-d', os.path.abspath(args.output + '/' + args.label + '.merged.bed'), os.path.abspath(data_path + '/' + args.genometype + '.srt.bed')], stdout=bedout)

		os.remove(os.path.abspath(args.output + '/' + args.label + '.merged.bed'))

		logging.info('Filtered .bed file ready')
		logging.info('Done')



def runInParallel(function, *arguments):

	proc = []

	for args in arguments:

		p = Process(target=function, args=args)
		p.start()
		proc.append(p)

	for p in proc:

		p.join()



def entropy(string): #Shannon entropy scanner

	prob = [float(string.count(c)) / len(string) for c in dict.fromkeys(list(string))]
	entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])

	return entropy


def modifier(coordinates): #fast way to remove None and substitute with closest number in list. Mantain 0-based coordinates

	start = next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:

			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates



def entropy_finder(sequence,coordinates,scansize,entropy_treshold): # get coordinates for intervals of certain scansize in a sequence in which entropy is lower than treshold 

	ind_start=0
	ind_end=scansize

	terminal_ind=len(sequence)-1

	hit=[]


	while terminal_ind > ind_end:

		if terminal_ind-ind_end >= scansize:

			if entropy(sequence[ind_start:ind_end]) < entropy_treshold:

				hit.append((coordinates[ind_start],coordinates[ind_end]))

			ind_start+=scansize
			ind_end+=scansize

		else:

			break

	return hit



def merge_intervals(intervals): #merge overlapping tuples in the same list

	sorted_by_lower_bound = sorted(intervals, key=itemgetter(0))
	merged = []

	for higher in sorted_by_lower_bound:
	
		if not merged:
			
			merged.append(higher)
		
		else:

			lower = merged[-1]

			if higher[0] <= lower[1]:
			
				upper_bound = max(lower[1], higher[1])
				merged[-1] = (lower[0], upper_bound)  
		
			else:
			
				merged.append(higher)

	return merged




def BScanner(bamfilein,bedfileout,scansize,entropy_treshold,call_treshold, dist_treshold): #entropy treshold was trained for a scansize of 20 bp: call treshold is the number of supporting entropy drop to consider this drop as true

	bamfile=pysam.AlignmentFile(bamfilein,'rb') #open the bamfile
	header=bamfile.header
	chromosomes_info=list(header.items())[1][1]

	with open (bedfileout, 'w') as fin:

		fin.write('#Bed file with regions in which an entropy lower than ' + str(entropy_treshold) + ' was found at least ' + str(call_treshold) + ' times. Minimum size of repetitive region is ' + str(dist_treshold) + '. Intervals overlapping in +- 500 bp range are merged' + '\n')

	classic_chrs = ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y']]

	chrom_dict=dict()

	for infos in chromosomes_info:

		if infos['SN'] in classic_chrs:

			chrom_dict[infos['SN']]=infos['LN']

	for chromosome in chrom_dict.keys():

		chr_array=np.zeros(chrom_dict[chromosome])

		for reads in bamfile.fetch(chromosome):

			if not reads.is_unmapped and not reads.is_secondary and not reads.is_supplementary and reads.mapping_quality >= 10: #exclude reads with mapping quality really low

				coordinates=modifier(reads.get_reference_positions(full_length=True))
				hits=entropy_finder(reads.seq,coordinates,scansize,entropy_treshold)

				for hit in hits:

					chr_array[hit[0]:hit[1]+1]+=1


		to_get= np.concatenate(np.where(chr_array>=call_treshold)).tolist() #at least the number of entropy drop specified in call_treshold must support the entropy drop in that point
		intervals=[]

		for k, g in itertools.groupby(enumerate(to_get),lambda x:x[0]-x[1]):

			group = (map(itemgetter(1),g))
			group = list(map(int,group))

			if len(group) >= dist_treshold: #exclude intervals that are shorter than dist_treshold

				intervals.append((group[0]-500,group[-1]+500)) # extend interval to the left and to the right, so that intervals in this range can be merged

		#merge overlapping intervals

		intervals=merge_intervals(intervals)

		with open (bedfileout, 'a') as fin:

			for inter in intervals:

					fin.write(chromosome + '\t' +  str(inter[0]) + '\t' + str(inter[1]) + '\n') 


	bamfile.close()
