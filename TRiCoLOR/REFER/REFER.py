#!/usr/bin/python env

#python 3 standard library

import sys
import os
import re
import glob
import subprocess
import itertools
import csv
import math
from collections import defaultdict
from operator import itemgetter
from multiprocessing import Process,Manager
from shutil import which
import logging
import datetime
import timeit


# additional libraries

import pyfaidx
import pysam
import editdistance
import pandas as pd


# import submodules

from .Helper import parser
from .Helper import finder
from .Helper import writer


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

		elif os.listdir(os.path.abspath(args.output)):

			print('Output folder is not empty. Clean the folder and re-run TRiCoLOR REFER')
			sys.exit(1)
	
	
	command_dict= vars(args)
	command_string= ','.join(("{}={}".format(*i) for i in command_dict.items())) 

	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR_REFER.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	

	#check the presence of needed external tools

	external_tools=['minimap2', 'samtools', 'bgzip', 'tabix', 'vcftools', 'bcftools']

	for tools in external_tools:

		if which(tools) is None:

			logging.error(tools + ' was not found as an executable command. Install ' + tools + ' and re-run TRiCoLOR REFER')
			sys.exit(1)


	## validate all inputs

	#check if the genome file exists, is readable and is in .fasta format

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>')) #genome .file starts with '>'

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

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile1)],stderr=open(os.devnull, 'wb'))

		except:

			logging.error('.bam file 1 could not be indexed')
			sys.exit(1)


	if not os.path.exists(os.path.abspath(args.bamfile2 + '.bai')):

		logging.info('Creating index for .bam file 2, as it was not found in folder...')

		try:

			subprocess.check_call(['samtools', 'index', os.path.abspath(args.bamfile2)],stderr=open(os.devnull, 'wb'))

		except:

			logging.error('.bam file 2 could not be indexed')
			sys.exit(1)



	mmi_abspath=os.path.abspath(os.path.dirname(args.genome))

	alfred_path=os.path.abspath(os.path.dirname(__file__) + '/alfred/bin/alfred')

	merging_path=os.path.abspath(os.path.dirname(__file__) + '/merging.sh')

	start_t=timeit.default_timer()

	logging.info('Analysis starts now')


	#write VCF header

	writer.VCF_headerwriter(args.bamfile1, args.bamfile2, args.samplename, command_string, os.path.abspath(args.output))


	ref=pyfaidx.Fasta(args.genome)
	b_in=Bed_Reader(args.bedfile)
	chromosomes_seen=set()
	it_ = iter(b_in)

	skipped=0
	ambiguous=0
	analyzed=0

	for i in range(b_in.length()):

		chromosome,start,end=next(it_)

		if chromosome not in chromosomes_seen:

			if chromosomes_seen: #set is not empty

				res=CleanResults(merging_path, list(chromosomes_seen)[-1], args.output, os.path.abspath(args.bamfile1), os.path.abspath(args.bamfile2)) #clean results for previous chromosome

				if type(res) == str:

					logging.error(res)
				
			chrom=ref[chromosome]
			ref_seq=chrom[:len(chrom)].seq
			chromosomes_seen.add(chromosome)

			#check if the .mmi reference exists, otherwise create it for the wanted chromosome

			if not os.path.exists(os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi')):

				if not os.access(mmi_abspath, os.W_OK): #check if we have write permissions on the directory

					logging.error('You do not have write permissions on the referene genome folder. This is necessary to create .mmi indexes')
					sys.exit(1)


				logging.info('Creating a .mmi index for ' + chromosome)

				try:

					subprocess.check_call(['samtools','faidx',os.path.abspath(args.genome),chromosome],stderr=open(os.devnull, 'wb'),stdout=os.path.abspath(mmi_abspath + '/' + chromosome + '.fa'))

				except:

					logging.error('Something went wrong with the creation of the .mmi chromosome index. Does your reference genome contain informations for ' + chromosome + ' ?')
					sys.exit(1)

				subprocess.call(['minimap2','-d', os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi'),os.path.abspath(mmi_abspath + '/' + chromosome + '.fa')]) #must work as everything has been checked
			
			mmi_ref=os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi')


		try:

			further = Ref_Repeats(ref_seq, chromosome, start, end, args.motif, args.times, args.maxmotif, args.overlapping, args.size, args.output)

			if further:

				manager = Manager()

				repetitions_h1 = manager.list()
				repetitions_h2 = manager.list()

				p1=Process(target=Haplo1_Repeats, args=(alfred_path, os.path.abspath(args.bamfile1), chromosome, start, end, args.coverage, args.motif, args.times,args.maxmotif, args.overlapping, args.size, args.editdistance, ref_seq, mmi_ref, args.output, i,repetitions_h1))
				p2=Process(target=Haplo2_Repeats, args=(alfred_path, os.path.abspath(args.bamfile2), chromosome, start, end, args.coverage, args.motif, args.times, args.maxmotif, args.overlapping, args.size, args.editdistance, ref_seq, mmi_ref, args.output, i,repetitions_h2))

				p1.start()
				p2.start()

				p1.join()
				p2.join()

				writer.VCF_writer(chromosome, further, ref_seq, repetitions_h1, os.path.abspath(args.output + '/haplotype1/' + str(i+1) + '.srt.bam'), repetitions_h2, os.path.abspath(args.output + '/haplotype2/' + str(i+1) + '.srt.bam'), os.path.abspath(args.output))

				if os.stat(os.path.abspath(args.output + '/haplotype1/' + str(i+1) + '.srt.bam')).st_size == 0:

					os.remove(os.path.abspath(args.output + '/haplotype1/' + str(i+1) + '.srt.bam'))
					os.remove(os.path.abspath(args.output + '/haplotype1/' + str(i+1) + '.srt.bam.bai'))


				if os.stat(os.path.abspath(args.output + '/haplotype2/' + str(i+1) + '.srt.bam')).st_size == 0:

					os.remove(os.path.abspath(args.output + '/haplotype2/' + str(i+1) + '.srt.bam'))
					os.remove(os.path.abspath(args.output + '/haplotype2/' + str(i+1) + '.srt.bam.bai'))

				analyzed+=1

									
			else:

				logging.info('Skipped ambiguous region ' + chromosome + ':' + str(start) + '-' + str(end))
				ambiguous +=1
				
				continue

		except:

			logging.exception('Something went wrong for region ' + chromosome + ':' + str(start) + '-' +str(end) + ". Log is attached below:")
			skipped +=1


	res=CleanResults(merging_path, list(chromosomes_seen)[-1], args.output, os.path.abspath(args.bamfile1), os.path.abspath(args.bamfile2)) #clean results at the end of the process

	if type(res) == str: #function is not supposed to return a message unless something unexpected happened

		logging.error(res)

	#.vcf file post-processing

	try:

		with open(os.path.abspath(args.output + '/TRiCoLOR.srt.vcf'), 'w') as srtvcfout:

			subprocess.check_call(['vcf-sort', os.path.abspath(args.output + '/TRiCoLOR.vcf')],stderr=open(os.devnull, 'wb'),stdout=srtvcfout) #outputted .vcf file is supposed to be sorted but sort, just in case

		os.remove(os.path.abspath(args.output + '/TRiCoLOR.vcf'))

		subprocess.check_call(['bgzip', os.path.abspath(args.output + '/TRiCoLOR.srt.vcf')])
		subprocess.check_call(['tabix', os.path.abspath(args.output + '/TRiCoLOR.srt.vcf.gz')])

		#also normalize

		subprocess.check_call(['bcftools', 'norm', '-f', os.path.abspath(args.genome), '-o', os.path.abspath(args.output + '/TRiCoLOR.norm.bcf'), '-O', 'b', os.path.abspath(args.output + '/TRiCoLOR.srt.vcf.gz')],stderr=open(os.devnull, 'wb'))
		subprocess.check_call(['tabix', os.path.abspath(args.output + '/TRiCoLOR.norm.bcf')])

	except:

		logging.exception('Something went wrong during .vcf processing. Log is attached below:')

	end_t=timeit.default_timer()
	elapsed=(end_t-start_t)/60 #convert time to minutes

	logging.info('Analysis completed in ' + str(elapsed) + ' minutes')
	logging.info('Total regions: ' + str(b_in.length()))
	logging.info('Regions analyzed: ' + str(analyzed))
	logging.info('Ambiguous regions skipped: ' + str(ambiguous))
	logging.info('Regions not analyzed: ' + str(skipped))




class Bed_Reader():

	def __init__(self,bedfile):

		self.bedfile=bedfile

	def __iter__(self):

		with open (self.bedfile, 'r') as bedin:

			for line in csv.reader(bedin, delimiter='\t'):

				if not line[0].startswith('#') and line !=[]: 

					if len(line) < 3:

						logging.error('.bed given to TRiCoLOR REFER -bed/--bedfile must be a .bed file with at least 3 fields: chromosome, start and end')
						sys.exit(1)

					else:

						yield (line[0], int(line[1]), int(line[2])) # exclude other fields if they are present

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



def isEmpty(list_obj): #recursive function to check for empty list

	if list_obj == []:

		return True

	else:

		return all((isinstance(sli, list) and isEmpty(sli)) for sli in list_obj)


	
def TableWriter(chromosome,repetitions_with_coord, out): #Table is in .bed (chromosome, start, end) format with header. Append to the existing table, so that we don't have to merge at the end, saving time

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



def ref_nestover(SortedIntervals, string): #for reference, consider largest between nested repetitions

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
	
		self_=list(finder.look_for_self(reps,wanted))

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

				finder.dfs(collapsed_ranges, visited, vertex, result, vertex)

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

		repetitions=list(finder.RepeatsFinder(wanted,kmer,times, maxmotif, overlapping))

		filtered=Reference_Filter(repetitions,wanted,size,start)

		if isEmpty(filtered):

			Table=EmptyTable(os.path.abspath(out_ + '/' + chromosome + '.repetitions.bed'))
			Table.write()

		else:

			TableWriter(chromosome,filtered,out_)
		
		return filtered


def Haplo1_Repeats(alfred_path, bamfile1, chromosome, start, end, coverage, kmer, times, maxmotif, overlapping, size, allowed, ref_seq, mmi_ref, out, iteration,repetitions_h1):

	out_=os.path.abspath(out+'/haplotype1')

	if not os.path.exists(out_):

		os.makedirs(out_)

	parser.Bamfile_Analyzer(bamfile1,chromosome,start,end, coverage, out_)

	if isEmpty(glob.glob(os.path.abspath(out_)+'/*.unaligned.fa')): #no informations for that region in this haplotype

		Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
		Table.write()

		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam'), 'w').close() #create empty
		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam.bai'), 'w').close() #create empty

		return

	else:
		parser.MSA(alfred_path, out_,mmi_ref)
		consensus_bams=glob.glob(os.path.abspath(out_+'/*consensus.srt.bam'))

		for bam in consensus_bams:

			coords,seq=finder.Get_Alignment_Positions(bam)

			if isEmpty(seq):

				Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
				Table.write()

				continue

			else:

				repetitions=list(finder.RepeatsFinder(seq,kmer,times, maxmotif, overlapping))

				if isEmpty(repetitions): #no repetitions found in the region

					if bam != consensus_bams[-1]: #not the last one

						continue

					else:

						Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
						Table.write()


				else:

					cor_coord_reps=finder.corrector(ref_seq, seq, repetitions, coords, size, allowed) #probably an exception here is needed

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



def Haplo2_Repeats(alfred_path, bamfile2, chromosome, start, end, coverage, kmer, times, maxmotif, overlapping, size,allowed, ref_seq, mmi_ref, out, iteration,repetitions_h2):


	out_=os.path.abspath(out+'/haplotype2')

	if not os.path.exists(out_):

		os.makedirs(out_)

	parser.Bamfile_Analyzer(bamfile2,chromosome,start,end, coverage, out_)

	if isEmpty(glob.glob(os.path.abspath(out_)+'/*.unaligned.fa')): #no informations for that region in this haplotype

		Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
		Table.write()

		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam'), 'w').close() #create empty
		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam.bai'), 'w').close() #create empty

		return

	else:

		parser.MSA(alfred_path, out_,mmi_ref)
		consensus_bams=glob.glob(os.path.abspath(out_+'/*consensus.srt.bam'))

		for bam in consensus_bams:

			coords,seq=finder.Get_Alignment_Positions(bam)

			if isEmpty(seq):

				Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
				Table.write()

				continue

			else:

				repetitions=list(finder.RepeatsFinder(seq,kmer,times, maxmotif, overlapping))

				if isEmpty(repetitions): #no repetitions found in the region

					if bam != consensus_bams[-1]: #not the last one

						continue

					else:

						Table=EmptyTable(os.path.abspath(out_ +'/' + chromosome + '.repetitions.bed'))
						Table.write()


				else:

					cor_coord_reps=finder.corrector(ref_seq, seq, repetitions, coords, size, allowed) #probably an exception here is needed

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



def CleanResults(merging_path, chromosome, out, bam1, bam2):

	#merge all the .srt.bam files for the two haplotypes

	out_=[os.path.abspath(out+j) for j in ['/haplotype1', '/haplotype2']]


	if not os.listdir(out_[0]): #directory is empty

		return

	else:

		try:

			subprocess.check_call(['bash', merging_path, bam1, os.path.abspath(out_[0]),chromosome],stderr=open(os.devnull, 'wb'))

		except:

			message= 'Something wrong while mergning .bam files for haplotype 1, ' + chromosome
			return message


	if not os.listdir(out_[0]): #directory is empty

		return

	else:

		try:

			subprocess.check_call(['bash', merging_path, bam2, os.path.abspath(out_[0]),chromosome],stderr=open(os.devnull, 'wb'))

		except:

			message= 'Something wrong while mergning .bam files for haplotype 2, ' + chromosome
			return message
