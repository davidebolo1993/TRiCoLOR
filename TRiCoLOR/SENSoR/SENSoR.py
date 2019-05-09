#!/usr/bin/python env

#python 3 standard library

import os
import sys
import math
import itertools
from operator import itemgetter
from multiprocessing import Process
from shutil import which
import logging
import subprocess


# additional libraries

import pysam
import numpy as np



def run(parser, args):


	if not os.path.exists(os.path.abspath(args.output)): #folder does not exist

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('Cannot create the output folder') #no write permission, probably			
			sys.exit(1)

	else: #path already exists

		if not os.access(os.path.abspath(args.output),os.W_OK): #folder exists but no write permissions

			print('Missing write permissions on the output folder')			
			sys.exit(1)
			
		elif os.listdir(os.path.abspath(args.output)): #folder exists but isn't empty. Do not overwrite results.

			print('The output folder is not empry. Specify another output folder or clean the previsouly chosen')
			sys.exit(1)
			

	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR_SENSoR.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

	external_tools=['samtools', 'bedops']

	for tools in external_tools: #this tools are required

		if which(tools) is None:

			logging.error(tools + ' cannot be executed. Install ' + tools + ' and re-run TRiCoLOR SENSoR')
			sys.exit(1)

	#check all inputs

	bams=args.bamfile[0]

	if len(bams) > 2:

		logging.error('TRiCoLOR supports haploid and diploid genomes only')


	for bam in bams:

		try:

			subprocess.check_call(['samtools','quickcheck', os.path.abspath(bam)],stderr=open(os.devnull, 'wb'))

		except:

			logging.error('BAM ' + bam + ' does not exist, is not readable or is not a valid BAM')
			sys.exit(1)


		if not os.path.exists(os.path.abspath(bam + '.bai')):

			logging.warning('Missing index for BAM ' + bam + '. Creating ...')

			try:

				subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile1)], stderr=open(os.devnull, 'wb'))

			except:

				logging.error('BAM ' + bam + ' could not be indexed')
				sys.exit(1)


	logging.info('Scansize: ' + str(args.scansize))
	logging.info('Entropy treshold: ' + str(args.entropy))
	logging.info('Call treshold: ' + str(args.call))
	logging.info('Entropy drops treshold: ' + str(args.length))
	logging.info('Ploidy: ' + str(len(bams)))
	logging.info('Cores: ' + str(len(bams)))

	if args.chromosomes is None:

		logging.info('Chromosomes: all')

	else:

		logging.info('Chromosomes: ' + '-'.join(x for x in args.chromosomes[0]))




	logging.info('Scanning ...')


	if len(bams) == 1:

		try:

			BScanner(bams[0], args.chromosomes, os.path.abspath(args.output + '/' + args.label + '.H1.bed'), args.scansize, args.entropy, args.call, args.length)

		except:

			logging.exception('Unexpected error while scanning. Log is below')
			sys.exit(1)

		logging.info('Done')

		logging.info('Writing final BED to output folder')


		with open(os.path.abspath(args.output + '/' + args.label + '.merged.srt.bed'), 'w') as srtbed1:
		
			subprocess.call(['sort-bed', os.path.abspath(args.output + '/' + args.label + '.H1.bed')],stdout=srtbed1, stderr=open(os.devnull, 'wb'))

		os.remove(os.path.abspath(args.output + '/' + args.label + '.H1.bed')) #remove unsorted


	else:


		try:

			runInParallel(BScanner, (bams[0], args.chromosomes, os.path.abspath(args.output + '/' + args.label + '.H1.bed'), args.scansize, args.entropy, args.call, args.length),(bams[1], args.chromosomes, os.path.abspath(args.output + '/' + args.label + '.H2.bed'), args.scansize, args.entropy, args.call,args.length))

		except:

			logging.exception('Unexpected error while scanning. Log is below')
			sys.exit(1)


		logging.info('Done')

		logging.info('Writing final BED to output folder')


		with open(os.path.abspath(args.output + '/' + args.label + '.H1.srt.bed'), 'w') as srtbed1:
		
			subprocess.call(['sort-bed', os.path.abspath(args.output + '/' + args.label + '.H1.bed')],stdout=srtbed1, stderr=open(os.devnull, 'wb'))


		with open(os.path.abspath(args.output + '/' + args.label + '.H2.srt.bed'), 'w') as srtbed2:
		
			subprocess.call(['sort-bed', os.path.abspath(args.output + '/' + args.label + '.H2.bed')],stdout=srtbed2, stderr=open(os.devnull, 'wb'))


		os.remove(os.path.abspath(args.output + '/' + args.label + '.H1.bed')) #remove unsorted
		os.remove(os.path.abspath(args.output + '/' + args.label + '.H2.bed')) #remove unsorted


		with open(os.path.abspath(args.output + '/' + args.label + '.merged.bed'), 'w') as bedout:

			subprocess.call(['bedops', '-m', os.path.abspath(args.output + '/' + args.label + '.H1.srt.bed'), os.path.abspath(args.output + '/' + args.label + '.H2.srt.bed')], stderr=open(os.devnull, 'wb'), stdout=bedout)


		os.remove(os.path.abspath(args.output + '/' + args.label + '.H1.srt.bed')) #remove unmerged
		os.remove(os.path.abspath(args.output + '/' + args.label + '.H2.srt.bed')) #remove unmerged


	if args.exclude is None:

		logging.info('No region to exclude')

	else:

		if not os.path.exists(os.path.abspath(args.exclude)):

			logging.warning('BED to -x/--exclude does not exists. No region excluded')

		else:

			try:

				with open(os.path.abspath(args.output + '/exclude.srt.bed'), 'w') as ebedout:

					subprocess.check_call(['sort-bed', os.path.abspath(args.exclude)],stdout=ebedout, stderr=open(os.devnull, 'wb'))

				with open(os.path.abspath(args.output + '/' + args.label + '.merged.bed.tmp'), 'w') as excludeout:

					subprocess.check_call(['bedops', '-d', os.path.abspath(args.output + '/' + args.label + '.merged.bed'), os.path.abspath(args.output + '/exclude.srt.bed')],stdout=excludeout, stderr=open(os.devnull, 'wb'))

				os.remove(os.path.abspath(args.output + '/' + args.label + '.merged.bed'))
				os.rename(os.path.abspath(args.output + '/' + args.label + '.merged.bed.tmp'),os.path.abspath(args.output + '/' + args.label + '.merged.bed'))


			except:

				logging.exception('Unexpected error while excluding regions from final BED. Log is below')
				sys.exit(1)


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



def modifier(coordinates): #fast way to remove None (soft-clipped coordinates) and substitute with closest number in list. Mantain 0-based coordinates

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




def BScanner(bamfilein, chromosomes, bedfileout,scansize,entropy_treshold,call_treshold, dist_treshold):


	bamfile=pysam.AlignmentFile(bamfilein,'rb')
	header=bamfile.header
	chromosomes_info=list(header.items())[1][1]
	chrom_dict=dict()

	for infos in chromosomes_info:

		if not chromosomes is None:

			if infos['SN'] in chromosomes[0]:

				chrom_dict[infos['SN']]=infos['LN']

		else:

			chrom_dict[infos['SN']]=infos['LN']


	for chromosome in chrom_dict.keys():

		chr_array=np.zeros(chrom_dict[chromosome])

		for reads in bamfile.fetch(chromosome):

			if not reads.is_unmapped and not reads.is_secondary and not reads.is_supplementary:

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
