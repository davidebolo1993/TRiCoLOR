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



def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Shannon entropy-based scanner to identify putative repetitions in whole bam files''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''')
	parser.add_argument('-b1', '--bam1', help='haplotype-resolved .bam 1 file to be scanned.', metavar='',required=True)
	parser.add_argument('-b2', '--bam2', help='haplotype-resolved .bam 2 file to be scanned.', metavar='',required=True)
	parser.add_argument('-s', '--scansize', help='scansize to use when scanning .bam files. Lower the scansize, higher the resolution,more it takes to scan', metavar='',default=20)
	parser.add_argument('-et', '--entropytreshold', help='entropy treshold to call for repetition. Must be trained before changing', metavar='',default=1.3)
	parser.add_argument('-ct', '--calltreshold', help='number of entropy drop needed to call for interval in .bed', metavar='',default=5)
	parser.add_argument('-l', '--label', help='label to identify the outputted .bed files', metavar='',required=True)	
	parser.add_argument('-o', '--output', help='path to where the resulting bed files will be saved', metavar='',required=True)
	parser.add_argument('-rt', '--reftype', help='Are you using Hg38 or Hg19?', metavar='',default='Hg38')	
	args = parser.parse_args()

	import timeit

	start=timeit.default_timer()

	runInParallel(BScanner, (args.bam1, os.path.abspath(args.output + '/' + args.label + '_hap1.bed'), args.scansize, args.entropytreshold, args.calltreshold),(args.bam2, os.path.abspath(args.output + '/' + args.label + '_hap2.bed'), args.scansize, args.entropytreshold, args.calltreshold))

	end=timeit.default_timer()
	elapsed=(end-start)/60

	print('Bam files scanned in', elapsed, 'minutes')

	try:

		assert(which('bedops') is not None)

		with open(os.path.abspath(args.output + '/' + args.label + '.bed'), 'w') as bedout:

			subprocess.call(['bedops', '--header', '-m', os.path.abspath(args.output + '/' + args.label + '_hap1.bed'), os.path.abspath(args.output + '/' + args.label + '_hap2.bed')], stdout=bedout)

		with open(os.path.abspath(args.output + '/' + args.label + '.filtered.tmp.bed'), 'w') as bedout:

			subprocess.call(['bedops', '-d', os.path.abspath(args.output + '/' + args.label + '.bed'), os.path.abspath('/home/bolognin/TRiCoLOR_py/Data/'+ args.reftype + 'Telomeres.sorted.bed')], stdout=bedout)

		with open(os.path.abspath(args.output + '/' + args.label + '.filtered.bed'), 'w') as bedout:

			subprocess.call(['bedops', '-d', os.path.abspath(args.output + '/' + args.label + '.tmp.bed'), os.path.abspath('/home/bolognin/TRiCoLOR_py/Data/'+ args.reftype + 'Centromeres.sorted.bed')], stdout=bedout)

		os.remove(os.path.abspath(args.output + '/' + args.label + '.filtered.tmp.bed'))


	except:

		sys.exit('bedops was not found as an executable command. Install ' + tools + ' and run, on the two .bed files, bedops --header -m/--merge. The resulting .bed file can be given, as input, to TRiCoLOR')



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





