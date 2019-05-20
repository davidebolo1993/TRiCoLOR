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
from collections import defaultdict, Counter
from operator import itemgetter
from multiprocessing import Process,Manager
from bisect import bisect_left,bisect_right
from shutil import which
import logging
import datetime


# additional modules

import pyfaidx
import pysam
import editdistance
import pandas as pd


# import submodules

from .Helper import parser
from .Helper import finder
from .Helper import writer


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

			print('The output folder is not empty. Specify another output folder or clean the previsouly chosen')
			sys.exit(1)
	

	command_dict= vars(args)
	
	notkey=['func']
	command_string= ' '.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey)
	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR_REFER.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')	
	logging.info('main=TRiCoLOR ' + command_string)

	#check the presence of needed external tools

	external_tools=['minimap2', 'samtools', 'bcftools']

	for tools in external_tools: 

		if which(tools) is None:

			logging.error(tools + ' cannot be executed. Install ' + tools + ' and re-run TRiCoLOR REFER')
			sys.exit(1)

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>'))

	except:

		logging.error('Reference file does not exist, is not readable or is not a valid FASTA')
		sys.exit(1)

	mmi_abspath=os.path.abspath(os.path.dirname(args.genome))
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

				subprocess.check_call(['samtools', 'index', os.path.abspath(bam)], stderr=open(os.devnull, 'wb'))

			except:

				logging.error('BAM ' + bam + ' could not be indexed')
				sys.exit(1)

	#external tools/scripts included

	alfred_path=os.path.abspath(os.path.dirname(__file__) + '/alfred/bin/alfred') #ask for a stand-alone consensus script
	merging_path=os.path.abspath(os.path.dirname(__file__) + '/merging.sh')

	#write infos for debugging

	if not args.precisemotif:

		if args.motif == 0 or args.motif == 1:

			logging.info('Repetition motif length: any')

		else:

			logging.info('Repetition motif length: at least ' + str(args.motif))

	else:

		if args.motif == 0:

			logging.info('Repetition motif length: any')

		else:

			logging.info('Repetition motif length: ' + str(args.motif))

	if not args.precisetimes:

		if args.times == 0 or args.times == 1:

			logging.info('Number of repetitions: any')

		else:

			logging.info('Number of repetitions: at least ' + str(args.times))

	else:

		if args.times == 1:

			logging.warning('--precisetimes is ignored if -t/--times is 1. Number of repetitions: any')

		else:

			logging.info('Number of repetitions: ' + str(args.times))
	
	logging.info('Check for overlapping repeated motif: ' + str(args.overlapping))

	regex=finder.RegexBuilder(args.motif,args.times,args.overlapping, args.precisemotif, args.precisetimes)		

	logging.info('Regex built: ' + regex)
	logging.info('Maximum repeated motif length: ' + str(args.maxmotif))
	logging.info('Allowed edit distance: ' + str(args.editdistance))
	logging.info('Minimum repeated region size: ' + str(args.size))
	logging.info('Ploidy: ' + str(len(bams)))
	logging.info('Cores: ' + str(len(bams)))

	if len(bams) == 2: #double ploidy

		writer.VCF_headerwriter(os.path.abspath(bams[0]), os.path.abspath(bams[1]), args.samplename, ','.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey), os.path.abspath(args.output))

	else: #single ploidy

		writer.VCF_headerwriter(os.path.abspath(bams[0]), None, args.samplename, ','.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey), os.path.abspath(args.output))

	ref=pyfaidx.Fasta(args.genome)
	b_in=Bed_Reader(args.bedfile)
	chromosomes_seen=set()
	it_ = iter(b_in)
	
	skipped=0
	ambiguous=0
	analyzed=0

	logging.info('Finding repetitions ...')

	for i in range(b_in.length()):

		chromosome,start,end=next(it_)

		if chromosome not in chromosomes_seen:

			if chromosomes_seen:

				if len(bams) == 2:

					res=CleanResults(merging_path, list(chromosomes_seen)[-1], args.output, os.path.abspath(bams[0]), os.path.abspath(bams[1]))

				else:

					res=CleanResults(merging_path, list(chromosomes_seen)[-1], args.output, os.path.abspath(bams[0]), None)

				if type(res) == str:

					logging.error(res)
					sys.exit(1)
							
			chromosomes_seen.add(chromosome)			
			chrom=ref[chromosome]
			ref_seq=chrom[:len(chrom)].seq			
			
			logging.info('Analyzing ' + str(chromosome))

			if not os.path.exists(os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi')):

				if not os.access(mmi_abspath, os.W_OK):

					logging.error('Missing write permissions on the reference genome folder. This is necessary to store chromosomes .mmi indexes')
					sys.exit(1)

				logging.info('Creating .mmi index for ' + chromosome)

				try:

					with open(os.path.abspath(mmi_abspath + '/' + chromosome + '.fa'), 'w') as chromout:

						subprocess.check_call(['samtools','faidx',os.path.abspath(args.genome),chromosome],stderr=open(os.devnull, 'wb'),stdout=chromout)

				except:

					logging.error('Unexpected error while creating .mmi index for current chromosome. Does your reference genome contain informations for this chromosome?')
					sys.exit(1)

				subprocess.call(['minimap2','-d', os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi'),os.path.abspath(mmi_abspath + '/' + chromosome + '.fa')],stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb')) #must work as everything has been checked
			
			mmi_ref=os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi')

		try:

			further = Ref_Repeats(ref_seq, chromosome, start, end, regex, args.maxmotif, args.size, args.output)

			if further is not None:	

				manager = Manager()

				repetitions_h1 = manager.list()
				repetitions_h2 = manager.list()

				p1=Process(target=Haplo1_Repeats, args=(alfred_path, os.path.abspath(bams[0]), chromosome, start, end, args.coverage, regex, args.maxmotif, args.size, args.editdistance, ref_seq, mmi_ref, args.output, i,repetitions_h1))

				if len(bams) == 2:

					p2=Process(target=Haplo2_Repeats, args=(alfred_path, os.path.abspath(bams[1]), chromosome, start, end, args.coverage, regex, args.maxmotif, args.size, args.editdistance, ref_seq, mmi_ref, args.output, i,repetitions_h2))

				else:

					p2=Process(target=Haplo2_Repeats, args=(alfred_path, None, chromosome, start, end, args.coverage, args.maxmotif, args.size, regex, args.editdistance, ref_seq, mmi_ref, args.output, i,repetitions_h2))

				p1.start()
				p2.start()

				p1.join()
				p2.join()

				#VCF

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

			logging.exception('Unexpected error for region ' + chromosome + ':' + str(start) + '-' +str(end) + ". Log is below")
			skipped +=1


	if len(bams) == 2:

		res=CleanResults(merging_path, list(chromosomes_seen)[-1], args.output, os.path.abspath(bams[0]), os.path.abspath(bams[1])) #clean results at the end of the process

	else:

		res=CleanResults(merging_path, list(chromosomes_seen)[-1], args.output, os.path.abspath(bams[0]), None) #clean results at the end of the process
		os.rmdir(os.path.abspath(args.output + '/haplotype2'))

	if type(res) == str:

		logging.error(res)
		sys.exit(1)

	#VCF post-processing

	try:

		subprocess.check_call(['bcftools', 'sort', '-o', os.path.abspath(args.output + '/TRiCoLOR.bcf'), '-O', 'b', os.path.abspath(args.output + '/TRiCoLOR.vcf')],stderr=open(os.devnull, 'wb'))		
		subprocess.check_call(['bcftools', 'index', os.path.abspath(args.output + '/TRiCoLOR.bcf')],stderr=open(os.devnull, 'wb'))

		os.remove(os.path.abspath(args.output + '/TRiCoLOR.vcf'))


	except:

		logging.exception('Unexpected error while processing VCF. Log is below')


	logging.info('Total regions: ' + str(b_in.length()))
	logging.info('Regions analyzed: ' + str(analyzed))
	logging.info('Ambiguous regions skipped: ' + str(ambiguous))
	logging.info('Regions not analyzed: ' + str(skipped))

	logging.info('Done')


class Bed_Reader():


	def __init__(self,bedfile):

		self.bedfile=bedfile

	def __iter__(self):

		with open (self.bedfile, 'r') as bedin:

			for line in csv.reader(bedin, delimiter='\t'):

				if not line[0].startswith('#') and line !=[]: 

					if len(line) < 3:

						logging.error('BED to TRiCoLOR REFER -bed/--bedfile must contain at least: chromosome, start, end (other fields are ignored)')
						sys.exit(1)

					else:

						yield (line[0], int(line[1]), int(line[2])) # exclude other fields if they are present

	def length(self):

		with open(self.bedfile, 'r') as bedin:

			size=sum(1 for _ in bedin if not _.startswith('#') and not _.strip()=='')

			return size


def isEmpty(list_obj): #recursive function to check for empty list


	if list_obj == []:

		return True

	else:

		return all((isinstance(sli, list) and isEmpty(sli)) for sli in list_obj)


def TableWriter(chromosome,repetitions_with_coord, out):


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


def SolveNestedR(SortedIntervals): #for reference, keep simply largest between nested


	extended=[]

	i=0

	while i < len(SortedIntervals):

		if extended==[]:

			extended.append(SortedIntervals[i])

		else:

			if extended[-1][2] >= SortedIntervals[i][1]: #the two intervals overlap

				if SortedIntervals[i][2] - SortedIntervals[i][1] > extended[-1][2] - extended[-1][1]:

					extended.remove(extended[-1])
					extended.append(SortedIntervals[i])
			else:

				extended.append(SortedIntervals[i])
		
		i+=1
	
	return extended


def ReferenceFilter(reference_reps,wanted,size,start): 


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

			new_reps=[(reps, start+val[0], start+val[-1]+len(reps)-1, len(val)) for val in list(result.values()) if val[-1]+len(reps)-val[0] >= size]
			corr_.extend(new_reps)

	s_corr_=sorted(corr_, key=itemgetter(1,2))
	mod_int=SolveNestedR(s_corr_)

	return mod_int



def Ref_Repeats(reference_seq, chromosome, start, end, regex, maxmotif, size, out):


	out_=os.path.abspath(out+'/reference')

	if not os.path.exists(out_):

		os.makedirs(out_)

	wanted=reference_seq[start-1:end]

	if 'N' in wanted: #skip ambiguous

		return

	else:

		repetitions=list(finder.RepeatsFinder(wanted,regex,maxmotif))
		filtered=ReferenceFilter(repetitions,wanted,size,start)

		if not isEmpty(filtered):

			TableWriter(chromosome,filtered,out_)
		
		return filtered


def Haplo1_Repeats(alfred_path, bamfile1, chromosome, start, end, coverage, regex, maxmotif, size, allowed, ref_seq, mmi_ref, out, iteration,repetitions_h1):


	out_=os.path.abspath(out+'/haplotype1')

	if not os.path.exists(out_):

		os.makedirs(out_)

	parser.Bamfile_Analyzer(bamfile1,chromosome,start,end, coverage, out_)

	if isEmpty(glob.glob(os.path.abspath(out_)+'/*.unaligned.fa')): #no informations for that region in this haplotype

		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam'), 'w').close() #create empty
		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam.bai'), 'w').close() #create empty

		return

	else:

		parser.MSA(alfred_path, out_,mmi_ref)
		consensus_bams=glob.glob(os.path.abspath(out_+'/*consensus.srt.bam'))

		for bam in consensus_bams:

			coords,seq=finder.Get_Alignment_Positions(bam)

			if isEmpty(seq):

				continue

			else:

				repetitions=list(finder.RepeatsFinder(seq,regex,maxmotif))

				if not isEmpty(repetitions): #no repetitions found in the region

					cor_coord_reps=finder.corrector(ref_seq, seq, repetitions, coords, size, allowed) #probably an exception here is needed

					if not isEmpty(cor_coord_reps): #size dimension exclude repetitions previously found

						TableWriter(chromosome, cor_coord_reps,out_)
						repetitions_h1.extend(cor_coord_reps)

		#merge and clean


		with open(os.path.abspath(out_ + '/FileToMerge.txt'), 'a') as fin:

			for file in consensus_bams:

				fin.write(file + '\n')

		subprocess.call(['samtools', 'merge', '-b', os.path.abspath(out_ + '/FileToMerge.txt'), os.path.abspath(out_+'/'+str(iteration+1) + '.srt.bam')])
		subprocess.call(['samtools', 'index', os.path.abspath(out_+'/'+str(iteration+1) + '.srt.bam')])
			
		os.remove(os.path.abspath(out_ + '/FileToMerge.txt'))

		for bams in consensus_bams:

			os.remove(bams)
			os.remove(bams + '.bai')


def Haplo2_Repeats(alfred_path, bamfile2, chromosome, start, end, coverage, regex, maxmotif, size, allowed, ref_seq, mmi_ref, out, iteration,repetitions_h2):


	out_=os.path.abspath(out+'/haplotype2')

	if not os.path.exists(out_):

		os.makedirs(out_)

	if bamfile2 is None:

		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam'), 'w').close() #create empty
		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam.bai'), 'w').close() #create empty

		return

	parser.Bamfile_Analyzer(bamfile2,chromosome,start,end, coverage, out_)

	if isEmpty(glob.glob(os.path.abspath(out_)+'/*.unaligned.fa')): #no informations for that region in this haplotype

		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam'), 'w').close() #create empty
		open(os.path.abspath(out_ +'/' + str(iteration +1) + '.srt.bam.bai'), 'w').close() #create empty

		return

	else:

		parser.MSA(alfred_path, out_,mmi_ref)
		consensus_bams=glob.glob(os.path.abspath(out_+'/*consensus.srt.bam'))

		for bam in consensus_bams:

			coords,seq=finder.Get_Alignment_Positions(bam)

			if isEmpty(seq):

				continue

			else:

				repetitions=list(finder.RepeatsFinder(seq,regex,maxmotif))

				if not isEmpty(repetitions): #no repetitions found in the region

					cor_coord_reps=finder.corrector(ref_seq, seq, repetitions, coords, size, allowed) #probably an exception here is needed

					if not isEmpty(cor_coord_reps): #size dimension exclude repetitions previously found

						TableWriter(chromosome, cor_coord_reps,out_)
						repetitions_h2.extend(cor_coord_reps)


		#merge and clean

		with open(os.path.abspath(out_ + '/FileToMerge.txt'), 'a') as fin:

			for file in consensus_bams:

				fin.write(file + '\n')

		subprocess.call(['samtools', 'merge', '-b', os.path.abspath(out_ + '/FileToMerge.txt'), os.path.abspath(out_+'/'+str(iteration+1) + '.srt.bam')])
		subprocess.call(['samtools', 'index', os.path.abspath(out_+'/'+str(iteration+1) + '.srt.bam')])
			
		os.remove(os.path.abspath(out_ + '/FileToMerge.txt'))

		for bams in consensus_bams:

			os.remove(bams)
			os.remove(bams + '.bai')


def CleanResults(merging_path, chromosome, out, bam1, bam2):


	out_=[os.path.abspath(out+j) for j in ['/haplotype1', '/haplotype2', '/reference']]


	if not os.listdir(out_[2]):

		open(os.path.abspath(out_[2] +'/' + chromosome + '.repetitions.empty.bed'), 'w').close()


	if not os.listdir(out_[0]):

		open(os.path.abspath(out_[0] +'/' + chromosome + '.repetitions.empty.bed'), 'w').close() 
		open(os.path.abspath(out_[0] +'/' + chromosome + '.empty.bam'), 'w').close() 
		open(os.path.abspath(out_[0] +'/' + chromosome + '.empty.bam.bai'), 'w').close() 

	else:

		try:

			subprocess.check_call(['bash', merging_path, bam1, os.path.abspath(out_[0]),chromosome],stderr=open(os.devnull, 'wb'))

		except BaseException as be:

			message= 'Unexpected error while merging BAM for haplotype 1:' + '\n' + be
			
			return message


	if bam2 is not None:

		if not os.listdir(out_[1]): #directory is empty

			open(os.path.abspath(out_[1] +'/' + chromosome + '.repetitions.empty.bed'), 'w').close()
			open(os.path.abspath(out_[1] +'/' + chromosome + '.empty.bam'), 'w').close() 
			open(os.path.abspath(out_[1] +'/' + chromosome + '.empty.bam.bai'), 'w').close() 

		else:

			try:

				subprocess.check_call(['bash', merging_path, bam2, os.path.abspath(out_[1]),chromosome],stderr=open(os.devnull, 'wb'))

			except BaseException as be:

				message= 'Unexpected error while merging BAM for haplotype 2:' + '\n' + be
				
				return message
