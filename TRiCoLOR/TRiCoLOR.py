#!/usr/bin/python env

import sys
import os
import re
import pyfaidx
import pysam
import glob
import subprocess
import itertools
import csv
import math
from collections import defaultdict
import editdistance
from operator import itemgetter
from multiprocessing import Process,Manager
from shutil import which
import pandas as pd
from bisect import bisect_left, bisect_right
from RepFinder import *
from BamParser import *
from VCFwriter import *
import argparse
from argparse import HelpFormatter
import logging
import datetime
import timeit


def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Tandem Repeats Caller fOr LOng Reads''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 
	
	required = parser.add_argument_group('Required arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='.fa',required=True)
	required.add_argument('-bed','--bedfile', help='.merged.srt.bed file generated during the sensing step or proprietary .bed file in the same format to use for STR searching', metavar='.bed',required=True)
	required.add_argument('-bam1','--bamfile1', help='haplotype-resolved .bam file, first haplotype',metavar='.bam',required=True)
	required.add_argument('-bam2','--bamfile2', help='haplotype-resolved .bam file, second haplotype',metavar='.bam',required=True)
	required.add_argument('-O','--output', help='path to where results will be saved',metavar='folder',required=True)
	
	algorithm = parser.add_argument_group('Regex-search algorithm')

	algorithm.add_argument('-m','--motif', type=int, help='size of the motif of the repetition: use 0 for any. Default to 0',metavar='',default=0)
	algorithm.add_argument('-t','--times', type=int, help='consencutive times a repetition must occur at least to be detected: use 0 for any. Default to 3',metavar='',default=3)
	algorithm.add_argument('-s','--size', type=int, help='what size repetitions - even fuzzy ones - must have at least to be called. Default to 10',metavar='',default=10)
	algorithm.add_argument('-o','--overlapping', type=str2bool, help='check for overlapping repetitions. Default to False',metavar='', default='False')
	algorithm.add_argument('-mm','--maxmotif', type=int, help='exclude motifs which length is greater than value. Default to 6',metavar='', default=6)
	algorithm.add_argument('-edit','--editdistance', type=int, help='allow specified number of insertions, deletions or substitutions in repetitive patterns. Default to 1',metavar='', default=1)

	utilities = parser.add_argument_group('Coverage')

	utilities.add_argument('-c', '--coverage', type=int, help='minimum number of reads to call consensus sequence', metavar='', default=3)

	parser.add_argument('-mmi','--mmiref', help='path to minimap2 .mmi chromosome-resolved references folder',metavar='', default=None)
	parser.add_argument('-sname','--samplename', help='sample name to use in .vcf header',metavar='',default='Sample')

	args = parser.parse_args()


	if not os.path.exists(os.path.abspath(args.output)):

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('It was not possible to create the results folder. Specify a path for which you have write permission')
			sys.exit(1)

	
	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

	command_dict= vars(args)
	command_string=",".join(("{}={}".format(*i) for i in command_dict.items())) 

	if os.path.exists(os.path.abspath(args.output + '/TRiCoLOR.vcf')):

		print('Specified output folder already contains TRiCoLOR results. Clean that folder and run TRiCoLOR again or specify a different folder')
		sys.exit(1)

	else:

		VCF_headerwriter(args.bamfile1, args.bamfile2, args.samplename, command_string, os.path.abspath(args.output))


	start_t=timeit.default_timer()

	logging.info('Analysis starts now')

	#check the presence of needed external tools

	external_tools=['samtools', 'minimap2', 'htsbox', 'bgzip', 'tabix', 'vcftools', 'bcftools']

	for tools in external_tools:

		if which(tools) is None:

			logging.error(tools + ' was not found as an executable command. Install ' + tools + ' and re-run TRiCoLOR')
			sys.exit(1)


	#check inputs validity

	#check if the genome file exists, is readable and is in .fasta format

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>'))

	except:

		logging.error('Reference file does not exist, is not readable or is not a valid .fasta file')
		sys.exit(1)


	#check if the two .bam files exist, are readable and are valid .bam files


	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bamfile1)],stderr=open(os.devnull, 'wb'))

	except:

		logging.error('.bam file 1 does not exist, is not readable or is not a valid .bam file')
		sys.exit(1)


	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bamfile2)],stderr=open(os.devnull, 'wb'))
		
	except:

		logging.error('.bam file 2 does not exist, is not readable or is not a valid .bam file')
		sys.exit(1)

	
	#look for the index of the  two .bam files

	if not os.path.exists(os.path.abspath(args.bamfile1 + '.bai')):

		logging.info('Creating index for .bam file 1, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile1)])

		except:

			logging.error('.bam1 file could not be indexed')
			sys.exit(1)


	if not os.path.exists(os.path.abspath(args.bamfile2 + '.bai')):

		logging.info('Creating index for .bam file 2, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile2)])

		except:

			logging.error('.bam2 file could not be indexed')
			sys.exit(1)

	#check write permissions on needed folders

	if not os.access(os.path.dirname(os.path.abspath(args.output)),os.W_OK):

		logging.error('You do not have write permissions on the directory in which results will be stored: you must specify a folder for which you have write permissions')
		sys.exit(1)

	if args.mmiref is not None:

		mmi_abspath=os.path.abspath(args.mmiref)

	else:

		mmi_abspath=os.path.dirname(os.path.abspath(args.genome))

	ref=pyfaidx.Fasta(args.genome)
	b_in=Bed_Reader(args.bedfile)
	chromosomes_seen=set()
	it_ = iter(b_in)

	skipped=0
	ambiguous=0

	for i in range(b_in.length()):

		chromosome, start, end=next(it_)

		if chromosome not in chromosomes_seen:

			if chromosomes_seen: #set is not empty

				res=CleanResults(list(chromosomes_seen)[-1], args.output, os.path.abspath(args.bamfile1), os.path.abspath(args.bamfile2)) #clean results for previous chromosome

				if type(res) == str:

					logging.error(res)
				
			chrom=ref[chromosome]
			ref_seq=chrom[:len(chrom)].seq
			chromosomes_seen.add(chromosome)


			#check if the .mmi reference exists, otherwise create it for the wanted chromosome

			if not os.path.exists(os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi')):

			#check if we have write permissions on reference directory

				if not os.access(mmi_abspath, os.W_OK):

					logging.error('You do not have write permissions on the reference folder: create a new folder with a copy of the reference and use that folder - or change permissions -')
					sys.exit(1)

				logging.info('Creating a .mmi index for ' + chromosome)

				try:

					subprocess.check_call(['samtools','faidx',os.path.abspath(args.genome),chromosome],stderr=open(os.devnull, 'wb'),stdout=os.path.abspath(mmi_abspath + '/' + chromosome + '.fa'))

				except:

					logging.error('Something went wrong with the creation of the .mmi chromosome index. Most likely your genome is in .fasta format but does not contain the information for ' + chromosome)
					sys.exit(1)

				subprocess.call(['minimap2','-d', os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi'),os.path.abspath(mmi_abspath + '/' + chromosome + '.fa')]) #must work if the previous step was done correctly

			mmi_ref=os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi')

		
		else:

			pass

		try:

			further = Ref_Repeats(ref_seq, chromosome, start, end, args.motif, args.times, args.maxmotif, args.overlapping, args.size, args.output)

			if further:

				manager = Manager()

				repetitions_h1 = manager.list()
				repetitions_h2 = manager.list()

				p1=Process(target=Haplo1_Repeats, args=(os.path.abspath(args.bamfile1), chromosome, start, end, args.coverage, args.motif, args.times,args.maxmotif, args.overlapping, args.size, args.editdistance, ref_seq, mmi_ref, args.output, i,repetitions_h1))
				p2=Process(target=Haplo2_Repeats, args=(os.path.abspath(args.bamfile2), chromosome, start, end, args.coverage, args.motif, args.times, args.maxmotif, args.overlapping, args.size, args.editdistance, ref_seq, mmi_ref, args.output, i,repetitions_h2))

				p1.start()
				p2.start()

				p1.join()
				p2.join()

				VCF_writer(chromosome, further, ref_seq, repetitions_h1, os.path.abspath(args.output + '/haplotype1/' + str(i+1) + '.srt.bam'), repetitions_h2, os.path.abspath(args.output + '/haplotype2/' + str(i+1) + '.srt.bam'), os.path.abspath(args.output))

				if os.stat(os.path.abspath(args.output + '/haplotype1/' + str(i+1) + '.srt.bam')).st_size == 0:

					os.remove(os.path.abspath(args.output + '/haplotype1/' + str(i+1) + '.srt.bam'))
					os.remove(os.path.abspath(args.output + '/haplotype1/' + str(i+1) + '.srt.bam.bai'))


				if os.stat(os.path.abspath(args.output + '/haplotype2/' + str(i+1) + '.srt.bam')).st_size == 0:

					os.remove(os.path.abspath(args.output + '/haplotype2/' + str(i+1) + '.srt.bam'))
					os.remove(os.path.abspath(args.output + '/haplotype2/' + str(i+1) + '.srt.bam.bai'))

									
			else:

				logging.info('Skipped ambiguous region ' + chromosome + ':' + str(start) + '-' + str(end))
				ambiguous +=1
				
				continue

		except:

			logging.exception('Something went wrong for ' + chromosome + ':' + str(start) + '-' +str(end))
			skipped +=1



	res=CleanResults(list(chromosomes_seen)[-1], args.output, os.path.abspath(args.bamfile1), os.path.abspath(args.bamfile2)) #clean results at the end of the process

	if type(res) == str:

		logging.error(res)


	#.vcf file post-processing

	try:

		with open(os.path.abspath(args.output + '/TRiCoLOR.srt.vcf'), 'w') as srtvcfout:

			subprocess.check_call(['vcf-sort', os.path.abspath(args.output + '/TRiCoLOR.vcf')],stderr=open(os.devnull, 'wb'),stdout=srtvcfout) #outputted .vcf file is supposed to be sorted

		os.remove(os.path.abspath(args.output + '/TRiCoLOR.vcf'))

		subprocess.check_call(['bgzip', os.path.abspath(args.output + '/TRiCoLOR.srt.vcf')])
		subprocess.check_call(['tabix', os.path.abspath(args.output + '/TRiCoLOR.srt.vcf.gz')])

		subprocess.check_call(['bcftools', 'norm', '-f', os.path.abspath(args.genome), '-o', os.path.abspath(args.output + '/TRiCoLOR.norm.srt.vcf.gz'), '-O', 'z', os.path.abspath(args.output + '/TRiCoLOR.srt.vcf.gz')],stderr=open(os.devnull, 'wb'))
		subprocess.check_call(['tabix', os.path.abspath(args.output + '/TRiCoLOR.norm.srt.vcf.gz')])

	except:

		logging.exception('Something went wrong during .vcf processing. Log is below:')

	end_t=timeit.default_timer()
	elapsed=end_t-start_t

	logging.info('Analysis completed in ' + str(elapsed) + ' seconds')
	logging.info('Number of regions: ' + str(b_in.length()))
	logging.info('Ambiguous regions: ' + str(ambiguous))
	logging.info('Other regions skipped: ' + str(skipped))



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


class Bed_Reader():

	def __init__(self,bedfile):

		self.bedfile=bedfile

	def __iter__(self):

		with open (self.bedfile, 'r') as bedin:

			for line in csv.reader(bedin, delimiter='\t'):

				if not line[0].startswith('#') and line !=[]: 

					if len(line) < 3:

						logging.error('.bed input must be a .bed file with at least 3 fields: chromosome, start and end')
						sys.exit(1)

					else:

						yield (line[0], int(line[1]), int(line[2])) # exclude other fields

	def length(self):

		with open(self.bedfile, 'r') as bedin:

			size=sum(1 for _ in bedin if not _.startswith('#') and not _.strip()=='')

			return size


class EmptyTable():

	def __init__ (self, tablepath):

		self.tablepath=tablepath


	def write(self):

		if os.path.exists(os.path.abspath(self.tablepath)):

			return

		else:

			with open(self.tablepath, 'w') as refout:

				Empty=pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Repeated Motif','Repetitions Number'])
				Empty.to_csv(refout, index=False, sep='\t')



def str2bool(v):

	if v.lower() == 'true':

		return True

	elif v.lower() == 'false':

		return False
	
	else:

		raise argparse.ArgumentTypeError('Boolean value expected for argument --overlapping. Use True or False.')


def isEmpty(list_obj): #recursive function to check for empty list

	if list_obj == []:

		return True

	else:

		return all((isinstance(sli, list) and isEmpty(sli)) for sli in list_obj)



	
def TableWriter(chromosome,repetitions_with_coord, out): #Table is in .bed (chromosome, start, end) format with header. Append to the existing table, so that we don't have to merge after

	seq=[el[0] for el in repetitions_with_coord]
	start=[el[1] for el in repetitions_with_coord]
	end=[el[2] for el in repetitions_with_coord]
	rep=[el[3] for el in repetitions_with_coord]
	chrom=[chromosome]*len(start)
	Table=pd.DataFrame({'Chromosome':chrom, 'Start':start,'End':end, 'Repeated Motif':seq,'Repetitions Number':rep},columns=['Chromosome', 'Start', 'End', 'Repeated Motif', 'Repetitions Number'])
	
	if os.path.exists(os.path.abspath(out + '/' + chromosome + '.repetitions.bed')):

		with open(os.path.abspath(out + '/' + chromosome + '.repetitions.bed'), 'a') as refout:

			Table.to_csv(refout ,sep='\t',index=False, header=False)

	else:

		with open(os.path.abspath(out + '/' + chromosome + '.repetitions.bed'), 'a') as refout:

			Table.to_csv(refout ,sep='\t',index=False)



def ref_nestover(SortedIntervals, string):

	extended=[]
	extended.append(SortedIntervals[0])

	i=1

	while i < len(SortedIntervals):

		if extended[-1][2] > SortedIntervals[i][1]: #the two intervals overlap

			if SortedIntervals[i][2] - SortedIntervals[i][1] > extended[-1][2] - extended[-1][1]:

				extended.remove(extended[-1])
				extended.append(SortedIntervals[i])
				i+=1

			else:

				i+=1

		else:

			extended.append(SortedIntervals[i])
			i+=1
	
	return extended



def Reference_Filter(reference_reps,wanted,size,start): 

	corr_=[]

	for reps in reference_reps:
	
		self_=list(look_for_self(reps,wanted))

		ranges=[]

		for i in range(len(self_)-1):

			if self_[i+1][1]-self_[i][1] == len(reps):

				ranges.append((self_[i][1],self_[i+1][1]))

		collapsed_ranges= defaultdict(list)
		
		for x, y in ranges:

			collapsed_ranges[x].append(y)
			collapsed_ranges[y].append(x)

		result = defaultdict(list)
		visited = set()
		
		for vertex in collapsed_ranges:

			if vertex not in visited:

				dfs(collapsed_ranges, visited, vertex, result, vertex)

		if len(result) !=0:

			new_reps=[(reps,val[0], val[-1]+len(reps)-1, len(val)) for val in list(result.values())]
			corr_.extend(new_reps)


	s_corr_=sorted(corr_, key=itemgetter(1,2))
	mod_int=ref_nestover(s_corr_, wanted)

	return [(a,start+b, start+c, d) for (a,b,c,d) in mod_int if len(a)*d >= size]


def Ref_Repeats(reference_seq, chromosome, start, end, kmer, times, maxmotif, overlapping, size, out):

	out_=os.path.abspath(out+'/reference')

	if not os.path.exists(out_):

		os.makedirs(out_)

	wanted=reference_seq[start-1:end] #pyfaidx way to get start-end


	if 'N' in wanted: #region with ambiguous bases

		Table=EmptyTable(os.path.abspath(out_ + '/' + chromosome + '.repetitions.bed'))
		Table.write()

		return False

	else:

		repetitions=list(RepeatsFinder(wanted,kmer,times, maxmotif, overlapping))

		filtered=Reference_Filter(repetitions,wanted,size,start)

		if isEmpty(filtered):

			Table=EmptyTable(os.path.abspath(out_ + '/' + chromosome + '.repetitions.bed'))
			Table.write()

		else:

			TableWriter(chromosome,filtered,out_)

		
		return filtered



def Haplo1_Repeats(bamfile1, chromosome, start, end, coverage, kmer, times, maxmotif, overlapping, size, allowed, ref_seq, mmi_ref, out, iteration,repetitions_h1):

	out_=os.path.abspath(out+'/haplotype1')

	if not os.path.exists(out_):

		os.makedirs(out_)

	seq,coord = Bamfile_Analyzer(bamfile1,chromosome,start,end, coverage)
	filseq,filcoord=InCommon(seq,coord)

	if isEmpty(filseq): #no informations for that region in this haplotype

		Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
		Table.write()

		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam'), 'w').close() #create empty
		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam.bai'), 'w').close() #create empty

		return

	else:

		Fasta_Generator(filseq,filcoord,out_)
		MA(out_,mmi_ref)
		consensus_bams=glob.glob(os.path.abspath(out_+'/*consensus.srt.bam'))

		for bam in consensus_bams:

			coords,seq=Get_Alignment_Positions(bam)

			if isEmpty(seq):

				Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
				Table.write()

				continue

			else:

				repetitions=list(RepeatsFinder(seq,kmer,times, maxmotif, overlapping))

				if isEmpty(repetitions): #no repetitions found in the region

					if bam != consensus_bams[-1]: #not the last one

						continue

					else:

						Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
						Table.write()


				else:

					cor_coord_reps=corrector(ref_seq, seq, repetitions, coords, size, allowed) #probably an exception here is needed

					if isEmpty(cor_coord_reps): #size dimension exclude repetitions previously found

						Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
						Table.write()

					else:

						TableWriter(chromosome, cor_coord_reps,out_)
						repetitions_h1.extend(cor_coord_reps)

		#merge and clean

		if len(consensus_bams) == 1:

			os.rename(consensus_bams[0], consensus_bams[0].replace('.'.join(consensus_bams[0].split('.',2)[:2]),os.path.abspath(out_+"/" +str(iteration+1))))
			os.remove(consensus_bams[0].replace('.bam', '.bam.bai'))

		else:

			with open(os.path.abspath(out_ + '/FileToMerge.txt'), 'a') as fin:

				for file in consensus_bams:

					fin.write(file + '\n')

			subprocess.call(['samtools', 'merge', '-b', os.path.abspath(out_ + '/FileToMerge.txt'), os.path.abspath(out_+'/'+str(iteration+1) + '.srt.bam')])
			os.remove(os.path.abspath(out_ + '/FileToMerge.txt'))

			for bams in consensus_bams:

				os.remove(bams)
				os.remove(bams.replace('.bam', '.bam.bai'))



def Haplo2_Repeats(bamfile2, chromosome, start, end, coverage, kmer, times, maxmotif, overlapping, size,allowed, ref_seq, mmi_ref, out, iteration,repetitions_h2):


	out_=os.path.abspath(out+'/haplotype2')

	if not os.path.exists(out_):

		os.makedirs(out_)

	seq,coord = Bamfile_Analyzer(bamfile2,chromosome,start,end, coverage)
	filseq,filcoord=InCommon(seq,coord)

	if isEmpty(filseq): #no informations for that region in this haplotype

		Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
		Table.write()

		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam'), 'w').close() #create empty
		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam.bai'), 'w').close() #create empty


		return

	else:

		Fasta_Generator(filseq,filcoord,out_)
		MA(out_,mmi_ref)
		consensus_bams=glob.glob(os.path.abspath(out_+'/*consensus.srt.bam'))

		for bam in consensus_bams:

			coords,seq=Get_Alignment_Positions(bam)

			if isEmpty(seq):

				Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
				Table.write()

				continue

			else:

				repetitions=list(RepeatsFinder(seq,kmer,times, maxmotif, overlapping))

				if isEmpty(repetitions): #no repetitions found in the region

					if bam != consensus_bams[-1]: #not the last one

						continue

					else:

						Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
						Table.write()

				else:

					cor_coord_reps=corrector(ref_seq, seq, repetitions, coords, size, allowed) #probably an exception here is needed

					if isEmpty(cor_coord_reps): #size dimension exclude repetitions previously found

						Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
						Table.write()

					else:

						TableWriter(chromosome, cor_coord_reps,out_)
						repetitions_h2.extend(cor_coord_reps)

		#merge and clean


		if len(consensus_bams) == 1:

			os.rename(consensus_bams[0], consensus_bams[0].replace('.'.join(consensus_bams[0].split('.',2)[:2]),os.path.abspath(out_+"/" +str(iteration+1))))
			os.remove(consensus_bams[0].replace('.bam', '.bam.bai'))

		else:

			with open(os.path.abspath(out_ + '/FileToMerge.txt'), 'a') as fin:

				for file in consensus_bams:

					fin.write(file + '\n')

			subprocess.call(['samtools', 'merge', '-b', os.path.abspath(out_ + '/FileToMerge.txt'), os.path.abspath(out_+'/'+str(iteration+1) + '.srt.bam')])
			os.remove(os.path.abspath(out_ + '/FileToMerge.txt'))

			for bams in consensus_bams:

				os.remove(bams)
				os.remove(bams.replace('.bam', '.bam.bai'))



def CleanResults(chromosome, out, bam1, bam2):

	#merge all the .srt.bam files for the two haplotypes

	out_=[os.path.abspath(out+j) for j in ['/reference', '/haplotype1', '/haplotype2']]


	if not os.listdir(out_[1]): #directory is empty

		return

	else:


		try:

			subprocess.check_call(['sh', '/home/bolognin/TRiCoLOR_py/Merging.sh', bam1, os.path.abspath(out_[1]),chromosome],stderr=open(os.devnull, 'wb'))

		except:

			message= 'Something wrong in merging for haplotype 1, ' + chromosome
			return message


	if not os.listdir(out_[2]): #directory is empty

		return

	else:

		try:

			subprocess.check_call(['sh', '/home/bolognin/TRiCoLOR_py/Merging.sh', bam2, os.path.abspath(out_[2]),chromosome],stderr=open(os.devnull, 'wb'))

		except:

			message= 'Something wrong in merging for haplotype 2, ' + chromosome
			return message




if __name__ == main():

	main()
