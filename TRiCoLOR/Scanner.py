#!/usr/bin/python env


import os
import math
import sys
import pysam
import numpy as np
from operator import itemgetter
import itertools
import argparse
from argparse import HelpFormatter
from multiprocessing import Process
from shutil import which
import subprocess
import timeit
import logging


def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Shannon entropy-based scanner to identify putative repetitions in whole bam files''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat)

	required = parser.add_argument_group('Required arguments')


	required.add_argument('-bam1', '--bamfile1', help='haplotype-resolved .bam 1 file to be scanned.',metavar='.bam',required=True)
	required.add_argument('-bam2', '--bamfile2', help='haplotype-resolved .bam 2 file to be scanned.', metavar='.bam',required=True)
	required.add_argument('-O', '--output', metavar='folder', help='where the resulting .bed files will be saved',required=True)
	
	algorithm = parser.add_argument_group('Entropy-scanning algorithm')

	algorithm.add_argument('-s', '--scansize', type=int, help='scansize (bp) to use when scanning .bamfiles. Lower the scansize, higher the resolution,more it takes to scan. Default to 20', metavar='',default=20)
	algorithm.add_argument('-et', '--entropytreshold', type=int, help='entropy treshold to call for repetition. Default to 1.3, trained with 20 bp scansize', metavar='',default=1.3)
	algorithm.add_argument('-ct', '--calltreshold', type=float, help='number of entropy drops needed to call for putative repetitions in .bed. Default to 5', metavar='',default=5)

	genome = parser.add_argument_group('Filter .bed file')

	parser.add_argument('-rt', '--reftype', help='Are you using Hg38 or Hg19? Default to Hg38', metavar='',default='Hg38')  

	args = parser.parse_args()

	start=timeit.default_timer()

	if not os.path.exists(os.path.abspath(args.output)):

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('It was not possible to create the results folder. Specify a path for which you have write permissions')
			sys.exit(1)


	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR_scanner.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

	logging.info('Analysis starts now')

	external_tools=['samtools', 'bedops']

	for tools in external_tools:

		if which(tools) is None:

			logging.error(tools + ' was not found as an executable command. Install ' + tools + ' and re-run TRiCoLOR')
			sys.exit(1)

	#check inputs validity

	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bamfile1)],stderr=open(os.devnull, 'wb'))

	except:

		logging.error('.bamfile 1 does not exist, is not readable or is not a valid .bam file')
		sys.exit(1)


	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bamfile2)],stderr=open(os.devnull, 'wb'))
		

	except:

		logging.error('.bamfile 2 does not exist, is not readable or is not a valid .bam file')
		sys.exit(1)


	if not os.path.exists(os.path.abspath(args.bamfile1 + '.bai')):

		logging.info('Creating index for .bamfile 1, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile1)])

		except:

			logging.error('.bamfile 1 could not be indexed')
			sys.exit(1)


	if not os.path.exists(os.path.abspath(args.bamfile2 + '.bai')):

		logging.info('Creating index for .bamfile 2, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile2)])

		except:

			logging.error('.bamfile 2 could not be indexed')
			sys.exit(1)


	#check permissions on output folder

	if not os.access(os.path.dirname(os.path.abspath(args.output)),os.W_OK):

		logging.error('You do not have write permissions on the directory in which results will be stored: you must specify a folder for which you have write permissions')
		sys.exit(1)

	try:

		runInParallel(BScanner, (args.bamfile1, os.path.abspath(args.output + '/' + args.label + '_hap1.bed'), args.scansize, args.entropytreshold, args.calltreshold),(args.bamfile2, os.path.abspath(args.output + '/' + args.label + '_hap2.bed'), args.scansize, args.entropytreshold, args.calltreshold))

	except:

		logging.exception('Something went wrong during the scanning process.')
		sys.exit(1)


	end=timeit.default_timer()
	elapsed=end-start

	logging.info('.bamfiles scanned in ' + str(elapsed) + ' seconds')


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

		with open(os.path.abspath(args.output + '/' + args.label + '.filtered.merged.tmp.bed'), 'w') as bedout:

			subprocess.call(['bedops', '-d', os.path.abspath(args.output + '/' + args.label + '.merged.bed'), os.path.abspath('/home/bolognin/TRiCoLOR_py/Data/'+ args.reftype + 'Telomeres.sorted.bed')], stdout=bedout)

		with open(os.path.abspath(args.output + '/' + args.label + '.filtered.merged.bed'), 'w') as bedout:

			subprocess.call(['bedops', '-d', os.path.abspath(args.output + '/' + args.label + '.filtered.merged.tmp.bed'), os.path.abspath('/home/bolognin/TRiCoLOR_py/Data/'+ args.reftype + 'Centromeres.sorted.bed')], stdout=bedout)

		os.remove(os.path.abspath(args.output + '/' + args.label + '.filtered.merged.tmp.bed'))

		logging.info('Merged and filtered .bed file ready. Done.')



class CustomFormat(HelpFormatter):

	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()




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

	terminal_ind=len(sequence)-1

	hit=[]


	while terminal_ind > ind_end:

		if terminal_ind-ind_end >= scansize:

			if entropy(sequence[ind_start:ind_end]) < entropy_treshold:

				hit.append((coordinates[ind_start],coordinates[ind_end])) #use 0-based coordinates

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





