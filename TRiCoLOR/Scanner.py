#!/usr/bin/python env


import os
import math
import sys
import pysam
import numpy as np
from operator import itemgetter
import itertools
import argparse
from multiprocessing import Process
from shutil import which
import subprocess
import timeit
import logging


def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Shannon entropy-based scanner to identify putative repetitions in whole bam files''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''')
	parser.add_argument('-b1', '--bam1', help='haplotype-resolved .bam 1 file to be scanned.', metavar='',required=True)
	parser.add_argument('-b2', '--bam2', help='haplotype-resolved .bam 2 file to be scanned.', metavar='',required=True)
	parser.add_argument('-s', '--scansize', type=int, help='scansize to use when scanning .bam files. Lower the scansize, higher the resolution,more it takes to scan', metavar='',default=20)
	parser.add_argument('-et', '--entropytreshold', type=int, help='entropy treshold to call for repetition. Must be trained before changing', metavar='',default=1.3)
	parser.add_argument('-ct', '--calltreshold', type=float, help='number of entropy drop needed to call for interval in .bed', metavar='',default=5)
	parser.add_argument('-l', '--label', type=str, help='label to identify the outputted .bed files', metavar='',required=True) 
	parser.add_argument('-o', '--output', help='path to where the resulting .bed files will be saved', metavar='',required=True)
	parser.add_argument('-rt', '--reftype', help='Are you using Hg38 or Hg19?', metavar='',default='Hg38')  
	args = parser.parse_args()

	start=timeit.default_timer()

	if not os.path.exists(os.path.abspath(args.output)):

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('It was not possible to create the results folder. Specify a path for which you have write permission')
			sys.exit(1)


	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR_scanner.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


	logging.info('Analysis starts now')

	external_tools=['samtools', 'bedops']

	for tools in external_tools:

		if which(tools) is None:

			logging.error(tools + ' was not found as an executable command. Install ' + tools + ' and re-run TRiCoLOR')
			sys.exit(1)

	#check inputs validity

	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bam1)],stderr=open(os.devnull, 'wb'))

	except:

		logging.error('.bam file 1 does not exist, is not readable or is not a valid .bam file')
		sys.exit(1)


	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bam2)],stderr=open(os.devnull, 'wb'))
		

	except:

		logging.error('.bam file 2 does not exist, is not readable or is not a valid .bam file')
		sys.exit(1)


	if not os.path.exists(os.path.abspath(args.bam1 + '.bai')):

		logging.info('Creating index for .bam file 1, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bam1)])

		except:

			logging.error('.bam1 file could not be indexed')
			sys.exit(1)


	if not os.path.exists(os.path.abspath(args.bam2 + '.bai')):

		logging.info('Creating index for .bam file 2, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bam2)])

		except:

			logging.error('.bam2 file could not be indexed')
			sys.exit(1)


	#check permissions on output folder

	if not os.access(os.path.dirname(os.path.abspath(args.output)),os.W_OK):

		logging.error('You do not have write permissions on the directory in which results will be stored: you must specify a folder for which you have write permissions')
		sys.exit(1)

	try:

		runInParallel(BScanner, (args.bam1, os.path.abspath(args.output + '/' + args.label + '_hap1.bed'), args.scansize, args.entropytreshold, args.calltreshold),(args.bam2, os.path.abspath(args.output + '/' + args.label + '_hap2.bed'), args.scansize, args.entropytreshold, args.calltreshold))

	except:

		logging.exception('Something went wrong during the scanning process.')
		sys.exit(1)


	end=timeit.default_timer()
	elapsed=end-start

	logging.info('.bam files scanned in', elapsed, 'seconds')



	with open(os.path.abspath(args.output + '/' + args.label + '_hap1.srt.bed'), 'w') as srtbed1:
	
		subprocess.call(['sort-bed', os.path.abspath(args.output + '/' + args.label + '_hap1.bed')],stdout=srtbed1)


	with open(os.path.abspath(args.output + '/' + args.label + '_hap2.srt.bed'), 'w') as srtbed2:
	
		subprocess.call(['sort-bed', os.path.abspath(args.output + '/' + args.label + '_hap2.bed')],stdout=srtbed2)


	os.remove(os.path.abspath(args.output + '/' + args.label + '_hap1.bed')) #remove unsorted
	os.remove(os.path.abspath(args.output + '/' + args.label + '_hap2.bed')) #remove unsorted


	with open(os.path.abspath(args.output + '/' + args.label + '.merged.bed'), 'w') as bedout:

		subprocess.call(['bedops', '-m', os.path.abspath(args.output + '/' + args.label + '_hap1.srt.bed'), os.path.abspath(args.output + '/' + args.label + '_hap2.srt.bed')], stdout=bedout)


	if args.reftype != 'Hg38' and args.reftype != 'Hg19':

		logging.warning('You are not using Hg38 or Hg19 as reference genomes. Merged .bed file is ready but centromeres and teleomeres regions were not filtered out.')
		sys.exit(1)


	else:

		with open(os.path.abspath(args.output + '/' + args.label + '.filtered.tmp.bed'), 'w') as bedout:

			subprocess.call(['bedops', '-d', os.path.abspath(args.output + '/' + args.label + '.bed'), os.path.abspath('/home/bolognin/TRiCoLOR_py/Data/'+ args.reftype + 'Telomeres.sorted.bed')], stdout=bedout)

		with open(os.path.abspath(args.output + '/' + args.label + '.filtered.merged.bed'), 'w') as bedout:

			subprocess.call(['bedops', '-d', os.path.abspath(args.output + '/' + args.label + '.tmp.bed'), os.path.abspath('/home/bolognin/TRiCoLOR_py/Data/'+ args.reftype + 'Centromeres.sorted.bed')], stdout=bedout)

		os.remove(os.path.abspath(args.output + '/' + args.label + '.filtered.tmp.bed'))

		logging.info('Merged and filtered .bed file ready. Done.')



def runInParallel(function, *arguments):

	proc = []

	for args in arguments:

		p = Process(target=function, args=args)
		p.start()
		proc.append(p)

	for p in proc:

		p.join()



def entropy(string): #Shannon Entropy

	prob = [float(string.count(c)) / len(string) for c in dict.fromkeys(list(string))]
	entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])

	return entropy


def modifier(coordinates): #fast way to remove None and substitute with closest number in list

	start = next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:

			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates


def entropy_finder(sequence,coordinates,scansize,entropy_treshold): # get coordinates for intervals of scansize bp in a sequence in which entropy is lower than treshold 

	ind_start=0
	ind_end=scansize

	hit=[]


	while len(sequence)-1 >= ind_start:


		if len(sequence) -1 >= ind_end:

			if entropy(sequence[ind_start:ind_end]) < entropy_treshold:

				hit.append((coordinates[ind_start]+1,coordinates[ind_end]+1)) 

				ind_start+=scansize
				ind_end+=scansize

			else:

				ind_start+=scansize
				ind_end+=scansize

		else:

			if len(sequence)-ind_end > scansize-int(scansize/4):

				if entropy(sequence[ind_start:len(sequence)]) < entropy_treshold:

					hit.append((coordinates[ind_start]+1, coordinates[-1]+1))

					break #reached the end

				else:

					break #reached the end
			
			else: 

				break #for short intervals the treshold for a sequence with low entropy can be different
	
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




def BScanner(bamfilein,bedfileout,scansize,entropy_treshold,call_treshold): #entropy treshold was trained for a scansize of 20/50 bp: call treshold is the number of supporting entropy drop to consider this as true

	bamfile=pysam.AlignmentFile(bamfilein,'rb') #open the bamfile
	header=bamfile.header
	chromosomes_info=list(header.items())[1][1]

	with open (bedfileout, 'w') as fin:

		fin.write('#Bed file with regions in which an entropy lower than ' + str(entropy_treshold) + ' was found at least ' + str(call_treshold) + ' times. Intervals overlapping in +- 500 bp range are merged' + '\n')

	classic_chrs = ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y']]

	chrom_dict=dict()

	for infos in chromosomes_info:

		if infos['SN'] in classic_chrs:

			chrom_dict[infos['SN']]=infos['LN']

	for chromosome in chrom_dict.keys():

		chr_array=np.zeros(chrom_dict[chromosome])

		for reads in bamfile.fetch(chromosome):

			if not reads.is_unmapped and not reads.is_secondary and reads.mapping_quality >= 10:

				coordinates=modifier(reads.get_reference_positions(full_length=True))
				hits=entropy_finder(reads.seq,coordinates,scansize,entropy_treshold)

				for hit in hits:

					chr_array[hit[0]:hit[1]+1]+=1

			#print(reads.query_name)


		to_get= np.concatenate(np.where(chr_array>=call_treshold)).tolist() #at least the number of entropy drop specified in call_treshold must support the entropy drop in that point
		intervals=[]

		for k, g in itertools.groupby(enumerate(to_get),lambda x:x[0]-x[1]):

			group = (map(itemgetter(1),g))
			group = list(map(int,group))
			intervals.append((group[0]-500,group[-1]+500)) # extend interval to the left and to the right, so that close intervals overlap and can be merged 

		#merge overlapping intervals

		intervals=merge_intervals(intervals)

		with open (bedfileout, 'a') as fin:

			for inter in intervals:

				if (inter[1] - inter[0])%2 == 0: # make sure to have only even intervals that can be simply splitted in 2 if too large

					fin.write(chromosome + '\t' +  str(inter[0]) + '\t' + str(inter[1]) + '\n') 

				else:

					fin.write(chromosome + '\t' +  str(inter[0]) + '\t' + str(inter[1]+1) + '\n')


	bamfile.close()




if __name__ == main():


	main()





