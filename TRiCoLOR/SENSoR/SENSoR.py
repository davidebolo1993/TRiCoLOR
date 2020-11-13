#!/usr/bin/python3 env

#python 3 standard library

import os
import sys
import math
import itertools
import multiprocessing
from operator import itemgetter
from datetime import datetime
from shutil import which

#additional modules

import pysam
import pybedtools
import numpy as np


class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	BAM = list()
	OUT = ''
	scansize=0
	entropy=0
	call=0
	length=0
	innerdistance=0
	outerdistance=0
	chromosomes=None
	exclude=None


def runInParallel(function, *arguments): #? AS THIS DOES NOT REQUIRE TOO MUCH TIME EVEN TO SCAN ENTIRE GENOMES AT HIGH COVERAGES (~10 HRS FOR INITIAL 56X BAM), SIMPLY RUN IN PARALLEL THE 2 HAPLOTYPES WHEN 2 HAPLOTYPES ARE PROVIDED

	'''
	Having a haplotype-resolved (splitted) BAM, this parallelize (one proc for each haplotype) BAM-scanning
	'''

	proc = []

	for args in arguments:

		p = multiprocessing.Process(target=function, args=args)
		p.start()
		proc.append(p)

	for p in proc:

		p.join()


def entropy(substr):

	'''
	Calculate Shannon entropy of a (sub)string

	'''

	prob = [float(substr.count(x)) / len(substr) for x in dict.fromkeys(list(substr))]
	entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])

	return entropy


def modifier(coordinates):

	'''
	Modify list of pysam coordinates. Substitute None with closest not-None
	'''

	start = next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:

			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates


def entropy_finder(sequence,coordinates,c):

	'''
	Generate a list of hits where calculated Shannon entropy is lower than trained treshold (for each sequence in BAM)
	'''

	ind_start=0
	ind_end=c.scansize

	terminal_ind=len(sequence)-1

	hit=[]

	while terminal_ind - ind_end >= c.scansize:

		if entropy(sequence[ind_start:ind_end]) < c.entropy:

			hit.append((coordinates[ind_start],coordinates[ind_end]))

		ind_start+=c.scansize
		ind_end+=c.scansize


	return hit


def BScanner_parallel(bamfilein, bedfileout, c):

	'''
	If a haplotype-resolved BAM, scan the 2 haplotypes in parallel. Else use BScanner

	'''

	bamfile=pysam.AlignmentFile(bamfilein,'rb')
	header=bamfile.header
	chromosomes_info=list(header.items())[1][1]
	chrom_dict=dict()

	for infos in chromosomes_info:

		if not c.chromosomes is None:

			if infos['SN'] in c.chromosomes[0] and infos['SN'] not in c.exclude:

				chrom_dict[infos['SN']]=infos['LN']

		else:

			if infos['SN'] not in c.exclude:

				chrom_dict[infos['SN']]=infos['LN']

	for chromosome in chrom_dict.keys():

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Processing chromosome ' + chromosome)

		chr_array=np.zeros(chrom_dict[chromosome])

		for reads in bamfile.fetch(chromosome):

			if not reads.is_unmapped and not reads.is_secondary and not reads.is_supplementary:

				coordinates=modifier(reads.get_reference_positions(full_length=True))
				hits=entropy_finder(reads.query_sequence,coordinates,c)

				for hit in hits:

					chr_array[hit[0]:hit[1]+1]+=1

		to_get= np.concatenate(np.where(chr_array>=c.call)).tolist()
		
		if len(to_get) != 0:

			intervals=[]

			for k, g in itertools.groupby(enumerate(to_get),lambda x:x[0]-x[1]):

				group = (map(itemgetter(1),g))
				group = list(map(int,group))

				if len(group) >= c.length:

					value=np.median(chr_array[group])
					intervals.append((group[0]-350,group[-1]+350, value))

			with open (bedfileout, 'a') as bedout:

				for inter in intervals:

						bedout.write(chromosome + '\t' +  str(inter[0]) + '\t' + str(inter[1]) + '\t' + str(inter[2]) + '\n') 

	bamfile.close()


def BScanner(c):

	'''
	If a HP-tagged BAM, assign reads from the 2 haplotypes to different groups. Else use BScanner_parallel

	'''

	bamfile=pysam.AlignmentFile(c.BAM[0] ,'rb')
	header=bamfile.header
	chromosomes_info=list(header.items())[1][1]
	chrom_dict=dict()


	for infos in chromosomes_info:

		if not c.chromosomes is None:

			if infos['SN'] in c.chromosomes[0] and infos['SN'] not in c.exclude:

				chrom_dict[infos['SN']]=infos['LN']

		else:

			if infos['SN'] not in c.exclude:

				chrom_dict[infos['SN']]=infos['LN']

	for chromosome in chrom_dict.keys():

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Processing chromosome ' + chromosome)

		chr_array1=np.zeros(chrom_dict[chromosome])
		chr_array2=np.zeros(chrom_dict[chromosome])

		for reads in bamfile.fetch(chromosome):

			if not reads.is_unmapped and not reads.is_secondary and not reads.is_supplementary:

				if not reads.has_tag('HP'):

					continue

				else:

					coordinates=modifier(reads.get_reference_positions(full_length=True))
					hits=entropy_finder(reads.query_sequence,coordinates,c)

					if reads.get_tag('HP') == 1:

						for hit in hits:

							chr_array1[hit[0]:hit[1]+1]+=1

					else:

						for hit in hits:

							chr_array2[hit[0]:hit[1]+1]+=1

		#first BED

		to_get= np.concatenate(np.where(chr_array1>=c.call)).tolist()
		
		if len(to_get) != 0:

			intervals=[]

			for k, g in itertools.groupby(enumerate(to_get),lambda x:x[0]-x[1]):

				group = (map(itemgetter(1),g))
				group = list(map(int,group))

				if len(group) >= c.length:

					value=np.median(chr_array1[group])
					intervals.append((group[0]-350,group[-1]+350, value))

			with open (os.path.abspath(c.OUT + '/H1.bed'), 'a') as bedout:

				for inter in intervals:

						bedout.write(chromosome + '\t' +  str(inter[0]) + '\t' + str(inter[1]) + '\t' + str(inter[2]) + '\n') 

		#second BED

		to_get= np.concatenate(np.where(chr_array2>=c.call)).tolist()

		if len(to_get) != 0:

			intervals=[]

			for k, g in itertools.groupby(enumerate(to_get),lambda x:x[0]-x[1]):

				group = (map(itemgetter(1),g))
				group = list(map(int,group))

				if len(group) >= c.length:

					value=np.median(chr_array2[group])
					intervals.append((group[0]-350,group[-1]+350, value))

			with open (os.path.abspath(c.OUT + '/H2.bed'), 'a') as bedout:

				for inter in intervals:

						bedout.write(chromosome + '\t' +  str(inter[0]) + '\t' + str(inter[1]) + '\t' + str(inter[2]) + '\n') 

	bamfile.close()


def run(parser, args):

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] TRiCoLOR SENSoR v1.1')

	'''
	Check arguments, run functions
	'''

	#fill container

	c.BAM=[os.path.abspath(x) for x in args.bamfile[0]]
	c.OUT=os.path.abspath(args.output)
	c.scansize=args.scansize
	c.entropy=args.entropy
	c.call=args.call
	c.length=args.length
	c.innerdistance=args.inner
	c.outerdistance=args.outer
	c.chromosomes=args.chromosomes
	c.exclude=args.exclude

	#main

	if not os.path.exists(c.OUT):

		try:

			os.makedirs(c.OUT)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Cannot create the output folder')
			sys.exit(1)

	else:

		if not os.access(os.path.abspath(c.OUT),os.W_OK):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Missing write permissions on the output folder')
			sys.exit(1)
			
		elif os.listdir(os.path.abspath(c.OUT)):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] The output folder is not empty: specify another output folder or clean the current one')
			sys.exit(1)

	#check if bedtools is in PATH. As from v1.1, TRiCoLOR does not perform subprocess calls to bedtools but this bedtools is still required in PATH

	if which('bedtools') is None:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] bedtools must be in PATH')
		sys.exit(1)
		
	if len(c.BAM) > 2:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] TRiCoLOR strictly requires a couple of splitted BAM haplotypes or a single HP-tagged BAM')
		sys.exit(1)

	for bam in c.BAM:

		try:

			pysam.quickcheck(bam)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Errror] BAM ' + bam + ' does not exist, is not readable or is not a valid BAM')
			sys.exit(1)

		if not os.path.exists(bam + '.bai'):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Warning] Missing ' + bam + ' index. Creating')

			try:

				pysam.index(bam)

			except:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Errror] BAM ' + bam + ' could not be indexed')
				sys.exit(1)

	skip_c=[]

	if c.exclude is not None:

		if not os.path.exists(os.path.abspath(c.exclude)):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] File containing chromosomes to exclude is provided but does not exist')

		else:

			with open(os.path.abspath(c.exclude)) as fin:

				for line in fin:

					skip_c.append(line.rstrip())

	c.exclude=skip_c

	if len(c.BAM) == 2:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Scanning haplotype-resolved BAM')

		try:

			runInParallel(BScanner_parallel, (c.BAM[0], os.path.abspath(c.OUT + '/H1.bed'), c),(c.BAM[1], os.path.abspath(c.OUT + '/H2.bed'), c))

		except Exception as e:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Unexpected error while scanning: ' + str(e))
			sys.exit(1)

	else:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Scanning haplotype-tagged BAM')

		try:

			BScanner(c)

		except Exception as e:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Unexpected error while scanning: ' + str(e))
			sys.exit(1)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Writing intervals to BED')

	#h1 BED

	bedfile=pybedtools.BedTool(os.path.abspath(c.OUT + '/H1.bed'))
	bedsrtd=bedfile.sort()
	bedmerged1=bedsrtd.merge(d=c.innerdistance, c=[4], o=['mean'])

	#h2 BED

	bedfile=pybedtools.BedTool(os.path.abspath(c.OUT + '/H2.bed'))
	bedsrtd=bedfile.sort()
	bedmerged2=bedsrtd.merge(d=c.innerdistance, c=[4], o=['mean'])

	#cat and sort/merge

	catted=bedmerged1.cat(bedmerged2, postmerge=True, c=[4,4,4], o=['mean', 'stdev', 'collapse'], d=c.outerdistance)
	header='CHROM\tSTART\tEND\tCOVMEAN\tCOVSTD\tCOVCOLLAPSED\n'
	catted.saveas(os.path.abspath(c.OUT + '/TRiCoLOR.srt.bed.gz'), compressed=True, trackline=header) #this are approximate regions. No need to be precise on 0-based start and 1-based end

	os.remove(os.path.abspath(c.OUT + '/H1.bed'))
	os.remove(os.path.abspath(c.OUT + '/H2.bed'))

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Done')
	sys.exit(0)
