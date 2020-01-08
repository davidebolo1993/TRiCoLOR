#!/usr/bin/python env

#python 3 standard library

import os
import sys
import math
import itertools
import multiprocessing
import logging
import subprocess
from operator import itemgetter
from shutil import which

# additional modules

import pysam
import numpy as np


def run(parser, args):


	if not os.path.exists(os.path.abspath(args.output)):

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('Cannot create the output folder')		
			sys.exit(1)

	else:

		if not os.access(os.path.abspath(args.output),os.W_OK):

			print('Missing write permissions on the output folder')			
			sys.exit(1)
			
		elif os.listdir(os.path.abspath(args.output)):

			print('The output folder is not empty. Specify another output folder or clean the one previsouly chosen')
			sys.exit(1)
			
	command_dict= vars(args)
	
	notkey=['func']
	command_string= ' '.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey)
	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR.SENSoR.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	
	print('Initialized .log file ' + os.path.abspath(args.output + '/TRiCoLOR.SENSoR.log'))

	logging.info('main=TRiCoLOR ' + command_string)
	external_tools=['samtools', 'bedtools', 'bedops']

	for tools in external_tools:

		if which(tools) is None:

			logging.error(tools + ' cannot be executed. Install ' + tools + ' and re-run TRiCoLOR SENSoR')
			exitonerror()

	bams=args.bamfile[0]

	if len(bams) != 2:

		logging.error('TRiCoLOR supports only diploid individuals')
		exitonerror()

	for bam in bams:

		try:

			subprocess.call(['samtools','quickcheck', os.path.abspath(bam)],stderr=open(os.devnull, 'wb'))

		except:

			logging.error('BAM ' + bam + ' does not exist, is not readable or is not a valid BAM')
			exitonerror()

		if not os.path.exists(os.path.abspath(bam + '.bai')):

			logging.warning('Missing index for BAM ' + bam + '. Creating ...')

			try:

				subprocess.call(['samtools', 'index', os.path.abspath(bam)], stderr=open(os.devnull, 'wb'))

			except:

				logging.error('BAM ' + bam + ' could not be indexed')
				exitonerror()

	logging.info('Scansize: ' + str(args.scansize))
	logging.info('Entropy treshold: ' + str(args.entropy))
	logging.info('Minimum suppporting calls: ' + str(args.call))
	logging.info('Minimum entropy drops length: ' + str(args.length))
	logging.info('Haplotypes: ' + str(len(bams)))
	logging.info('Cores: ' + str(len(bams)))

	if args.chromosomes is None:

		logging.info('Chromosomes: all')

	else:

		logging.info('Chromosomes: ' + '-'.join(x for x in args.chromosomes[0]))

	header='#CHROM\tSTART\tEND\tMEANCOV\tSTDEVCOV\tCOLLAPSEDCOV\n'

	with open(os.path.abspath(args.output + '/TRiCoLOR.srt.bed'), 'w') as finalbedout:

		finalbedout.write(header)

	logging.info('Scanning ...')
	print('Scanning ...')

	try:

		runInParallel(BScanner, (bams[0], args.chromosomes, os.path.abspath(args.output + '/H1.bed'), args.scansize, args.entropy, args.call, args.length),(bams[1], args.chromosomes, os.path.abspath(args.output + '/H2.bed'), args.scansize, args.entropy, args.call,args.length))

	except:

		logging.exception('Unexpected error while scanning. Log is below')
		exitonerror()

	logging.info('Done')
	logging.info('Writing final BED to output folder ...')

	with open(os.path.abspath(args.output + '/H1.srt.bed'), 'w') as srtbed1:
		
		subprocess.call(['bedtools', 'sort', '-i', os.path.abspath(args.output + '/H1.bed')],stdout=srtbed1, stderr=open(os.devnull, 'wb'))

	with open(os.path.abspath(args.output + '/H2.srt.bed'), 'w') as srtbed2:
		
		subprocess.call(['bedtools', 'sort', '-i', os.path.abspath(args.output + '/H2.bed')],stdout=srtbed2, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(args.output + '/H1.bed'))
	os.remove(os.path.abspath(args.output + '/H2.bed'))

	with open(os.path.abspath(args.output + '/H1.srt.merged.bed'), 'w') as mergedbed1:
		
		subprocess.call(['bedtools', 'merge', '-i', os.path.abspath(args.output + '/H1.srt.bed'), '-c', '4,4,4', '-o', 'mean,stdev,collapse', '-d', str(args.innerdistance)],stdout=mergedbed1, stderr=open(os.devnull, 'wb'))

	with open(os.path.abspath(args.output + '/H2.srt.merged.bed'), 'w') as mergedbed2:
		
		subprocess.call(['bedtools', 'merge', '-i', os.path.abspath(args.output + '/H2.srt.bed'), '-c', '4,4,4', '-o', 'mean,stdev,collapse', '-d', str(args.innerdistance)],stdout=mergedbed2, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(args.output + '/H1.srt.bed'))
	os.remove(os.path.abspath(args.output + '/H2.srt.bed'))

	with open(os.path.abspath(args.output + '/TRiCoLOR.srt.bed.tmp'), 'w') as bedout:

		subprocess.call(['bedops', '-u', os.path.abspath(args.output + '/H1.srt.merged.bed'), os.path.abspath(args.output + '/H2.srt.merged.bed')], stderr=open(os.devnull, 'wb'), stdout=bedout)

	os.remove(os.path.abspath(args.output + '/H1.srt.merged.bed'))
	os.remove(os.path.abspath(args.output + '/H2.srt.merged.bed'))

	with open(os.path.abspath(args.output + '/TRiCoLOR.srt.bed'), 'a') as finalbedout:

		subprocess.call(['bedtools', 'merge', '-i', os.path.abspath(args.output + '/TRiCoLOR.srt.bed.tmp'), '-c', '4,4,4', '-o', 'mean,stdev,collapse', '-d', str(args.outerdistance)],stderr=open(os.devnull, 'wb'), stdout=finalbedout)

	os.remove(os.path.abspath(args.output + '/TRiCoLOR.srt.bed.tmp'))

	logging.info('Done')
	print('Done')


##FUNCTIONS


def exitonerror():

	print('An error occured. Check the log file for more details')
	sys.exit(1)


def runInParallel(function, *arguments): #? AS THIS DOES NOT REQUIRE TOO MUCH TIME EVEN TO SCAN ENTIRE GENOMES AT HIGH COVERAGES (~10 HRS FOR INITIAL 56X BAM), SIMPLY RUN IN PARALLEL THE 2 HAPLOTYPES WHEN 2 HAPLOTYPES ARE PROVIDED. DO WE NEED THIS TO BE EVEN FASTER? NOT A PRIORITY.


	proc = []

	for args in arguments:

		p = multiprocessing.Process(target=function, args=args)
		p.start()
		proc.append(p)

	for p in proc:

		p.join()


def entropy(string):


	prob = [float(string.count(c)) / len(string) for c in dict.fromkeys(list(string))]
	entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])

	return entropy


def modifier(coordinates):


	start = next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:

			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates


def entropy_finder(sequence,coordinates,scansize,entropy_treshold):


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

		to_get= np.concatenate(np.where(chr_array>=call_treshold)).tolist()
		intervals=[]

		for k, g in itertools.groupby(enumerate(to_get),lambda x:x[0]-x[1]):

			group = (map(itemgetter(1),g))
			group = list(map(int,group))

			if len(group) >= dist_treshold:

				value=np.median(chr_array[group])
				intervals.append((group[0]-350,group[-1]+350, value))

		with open (bedfileout, 'a') as bedout:

			for inter in intervals:

					bedout.write(chromosome + '\t' +  str(inter[0]) + '\t' + str(inter[1]) + '\t' + str(inter[2]) + '\n') 

	bamfile.close()
