#!/usr/bin/python3 env

#python 3 standard library

import sys
import os
import re
import subprocess
import itertools
import math
import statistics
import multiprocessing
from datetime import datetime,date
from collections import defaultdict,Counter
from operator import itemgetter
from bisect import bisect_left,bisect_right
from shutil import which

#additional modules

import pyfaidx
import pysam
import pybedtools
import mappy as mp
import editdistance
import numpy as np

#import submodules

from .Helper import parser
from .Helper import finder
from .Helper import writer


class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	BAM = list()
	OUT = ''
	REF = ''
	BED = ''
	match=0
	mismatch=0
	gapopen=0
	gapextend=0
	motif=0
	maxmotif=0
	times=0
	size=0
	editdistance=0
	overlapping=False
	precisemotif=False
	precisetimes=False
	coverage=0
	softclipping=0
	samplename=''
	threads=0
	mmidir=''
	exclude=''

	#add other stuff that will be used later

	refseq=''
	regex=''
	aligner=None
	Cpath=os.path.abspath(os.path.dirname(__file__) + '/consensus') #this is the default location. Should not be moved


def atoi(text):

	'''
	Convert text to integers
	'''

	return int(text) if text.isdigit() else text


def natural_keys(text):

	'''
	Natural sort, as bcftools wants
	'''
	
	return [ atoi(c) for c in re.split(r'(\d+)', text)]


def Chunks(l,n):

	'''
	Split list in chunks based on number of threads
	'''

	return [l[i:i+n] for i in range(0, len(l), n)]


def ReferenceReps(s,c):

	'''
	Find repeats in reference sequence
	'''

	chromosome,start,end=s.chrom,s.start,s.end
	wanted=c.refseq[start-1:end]

	if 'N' in wanted:

		return None

	else:

		repetitions=finder.RepeatsFinder(wanted,c)
		filtered=finder.ReferenceFilter(chromosome,repetitions,wanted,c,start)
		
		return filtered


def PyCoord(CIGOP,ref_pointer):

	'''
	Extract coords from cigar operations, in the style of pysam (i.e. read.get_reference_positions)
	'''

	s=ref_pointer
	coords=[]

	for i,op in enumerate(list(CIGOP)):

		if i==0 and op == 'M':

			coords.append(s)

		else:

			if op == 'M':
				
				s+=1
				coords.append(s)

			elif op == 'I' or op == 'S':

				coords.append(s) #this also avoids use of modifier (old TRiCoLOR behaviour)

			elif op == 'D':

				s+=1

	return coords


def Map(c,fastain,start,end,cov):

	'''
	Map consensus to reference
	'''

	seq=''
	name=''
	Aldict=dict()

	with open(fastain, 'r') as fin:

		for line in fin:

			if line[0] != '>':

				seq+=line.rstrip()

			else:

				name+=line.rstrip()

	for hit in c.aligner.map(seq,MD=True,cs=True): #store also MD and cs (maybe one wants to re-use)

		if hit.is_primary: #keep only primary alignment (no need to check for forward orientation, consensus are built on forward sequences)

			if hit.r_st <= start and hit.r_en >= end: #keep only if read aligns back properly

				clip = ['' if x == 0 else '{}S'.format(x) for x in (hit.q_st, len(seq) - hit.q_en)] #calculate soft clipped bases, if any
				
				clipped=sum(x for x in (hit.q_st, len(seq) - hit.q_en))
				clippedperc=clipped/len(seq)*100 #clipped percentage

				if clippedperc < c.softclipping:

					cigstr = ''.join((clip[0], hit.cigar_str, clip[1]))
					cig_split = [''.join(x) for _, x in itertools.groupby(cigstr, key=str.isdigit)] #split by operation
					cigop=[cig_split[n:n+2] for n in range(0,len(cig_split),2)] #group by operation
					cigop_conv=''.join([str(x[1])*int(x[0]) for x in cigop]) #convert to string, one char for each operation

					Aldict['QNAME'] = name[1:] + '_' + hit.ctg + '_' + str(hit.r_st)
					Aldict['RNAME'] = hit.ctg
					Aldict['POS'] = hit.r_st
					Aldict['MAPQ'] = hit.mapq
					Aldict['CIGAR'] = cigstr
					Aldict['SEQ'] = seq
					Aldict['MD'] = hit.MD
					Aldict['cs'] = hit.cs
					Aldict['coords'] = PyCoord(cigop_conv,hit.r_st)
					Aldict['coverage'] = cov

	return Aldict


def BamW(header,BAMsegments,bamout):

	'''
	Write aligned segments in BAM format
	'''

	with pysam.AlignmentFile(bamout, mode='wb', header=header) as bout:

		for segments in BAMsegments:

			if segments != {}:

				s = pysam.AlignedSegment(bout.header)
				s.is_unmapped=False #'cause unmapped reads were skipped
				s.is_reverse=False #'cause reads are all translated to forward orientation
				s.is_secondary=False #'cause only primary alignments were retained
				s.query_name=segments['QNAME']
				s.reference_name=segments['RNAME']
				s.reference_start=segments['POS']
				s.mapping_quality=segments['MAPQ']
				s.cigarstring=segments['CIGAR']
				s.query_sequence=segments['SEQ']
				s.set_tags([('MD',segments['MD'], 'Z'), ('cs', segments['cs'], 'Z')])
				bout.write(s)


def HaploReps(s,bamfile,out,c,processor):

	'''
	Find repeats in haplotype-resvolved BAM
	'''

	chromosome,start,end=s.chrom,s.start,s.end
	cov_inbam=parser.Bamfile_Analyzer(bamfile,chromosome,start,end,c,out,processor) #this also generate fasta for consensus
	file=os.path.abspath(out + '/' + processor + '.unaligned.fa')

	Aldict={}
	true_repetitions=[]

	if os.path.exists(file):

		#generate consensus .fa
		consfa=os.path.abspath(out + '/' + processor + '.consensus.fa')
		
		with open(consfa, 'w') as fout:

			subprocess.call([c.Cpath, str(c.match), str(c.mismatch), str(c.gapopen), str(c.gapextend), file], stdout=fout, stderr=open(os.devnull, 'wb'))

		os.remove(file)
		Aldict=Map(c,consfa,start,end,cov_inbam)
		os.remove(consfa)
		
		if Aldict != {}:

			repetitions=finder.RepeatsFinder(Aldict['SEQ'],c)

			if repetitions != set():

				true_repetitions=finder.corrector(chromosome,c,Aldict['SEQ'],repetitions,Aldict['coords'])


	return Aldict,true_repetitions


def HaploReps_Single(s,bamfile,c,processor):

	'''
	Find repeats in haplotype-tagged BAM
	'''

	chromosome,start,end=s.chrom,s.start,s.end
	cov_inbam1,cov_inbam2=parser.Bamfile_Analyzer_Single(bamfile,chromosome,start,end,c,processor)
	
	file1=os.path.abspath(c.OUT + '/haplotype1/' + processor + '.unaligned.fa')
	file2=os.path.abspath(c.OUT + '/haplotype2/' + processor + '.unaligned.fa')

	Aldict1,Aldict2={},{}
	true_repetitions1,true_repetitions2=[],[]

	#H1
	
	if os.path.exists(file1):

		#generate consensus .fa
		consfa1=os.path.abspath(c.OUT + '/haplotype1/' + processor + '.consensus.fa')
		
		with open(consfa1, 'w') as fout:

			subprocess.call([c.Cpath, str(c.match), str(c.mismatch), str(c.gapopen), str(c.gapextend), file1], stdout=fout, stderr=open(os.devnull, 'wb'))

		os.remove(file1)
		Aldict1=Map(c,consfa1,start,end,cov_inbam1)
		os.remove(consfa1)
		
		if Aldict1 != {}:

			repetitions1=finder.RepeatsFinder(Aldict1['SEQ'],c)

			if repetitions1 != set():
			
				true_repetitions1=finder.corrector(chromosome, c,Aldict1['SEQ'],repetitions1,Aldict1['coords'])

	#H2
	
	if os.path.exists(file2):

		#generate consensus .fa
		consfa2=os.path.abspath(c.OUT + '/haplotype2/' + processor + '.consensus.fa')
		
		with open(consfa2, 'w') as fout:

			subprocess.call([c.Cpath, str(c.match), str(c.mismatch), str(c.gapopen), str(c.gapextend), file2], stdout=fout, stderr=open(os.devnull, 'wb'))

		os.remove(file2)
		Aldict2=Map(c,consfa2,start,end,cov_inbam2)
		os.remove(consfa2)

		if Aldict2 != {}:

			repetitions2=finder.RepeatsFinder(Aldict2['SEQ'],c)

			if repetitions2 != set():
			
				true_repetitions2=finder.corrector(chromosome, c,Aldict2['SEQ'],repetitions2,Aldict2['coords'])

	return Aldict1,Aldict2,true_repetitions1,true_repetitions2


def Runner(sli,processor,c,Rrep,H1rep,H2rep,H1bam,H2bam,VCFvariants):


	'''
	Parallelize repeats detection
	'''

	for s in sli:

		try:

			pR=ReferenceReps(s,c)

			if pR is not None: #if reference does not contain Ns, basically
		
					Rrep.append(pR)

					if len(c.BAM) == 2:

						out1=os.path.abspath(c.OUT + '/haplotype1')
						out2=os.path.abspath(c.OUT + '/haplotype2')

						D1,R1=HaploReps(s,c.BAM[0],out1,c,processor)
						D2,R2=HaploReps(s,c.BAM[1],out2,c,processor)
						
					else:

						D1,D2,R1,R2=HaploReps_Single(s,c.BAM[0],c,processor)

					#append results to shared manger lists

					H1rep.append(R1)
					H2rep.append(R2)
					H1bam.append(D1)
					H2bam.append(D2)
			
					if D1 == {}:

						D1['SEQ'],D1['coords'],D1['coverage'],D1['MAPQ'] = [],[],[],[]

					if D2 == {}:

						D2['SEQ'],D2['coords'],D2['coverage'],D2['MAPQ'] = [],[],[],[]
					 
					variants=writer.VCF_writer(s.chrom,[x[1:] for x in pR],c.refseq,[x[1:] for x in R1],D1['SEQ'],D1['coords'],D1['coverage'],D1['MAPQ'],[x[1:] for x in R2],D2['SEQ'],D2['coords'],D2['coverage'],D2['MAPQ'])

					if variants != []:

						VCFvariants.extend(variants)

		except Exception as e:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Unexpected error while processing region ' + s.chrom +':'+str(s.start)+'-'+str(s.end) +': ' + str(e))
			continue #skip current iteration. Better than exiting (?)


def run(parser, args):

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] TRiCoLOR REFER v1.1')

	'''
	Check arguments, run functions
	'''

	#fill container

	c.OUT=os.path.abspath(args.output)
	c.REF=os.path.abspath(args.genome)
	c.BED=os.path.abspath(args.bedfile)

	c.match=args.match
	c.mismatch=args.mismatch
	c.gapopen=args.gapopen
	c.gapextend=args.gapextend
	c.motif=args.motif
	c.maxmotif=args.maxmotif
	c.times=args.times
	c.size=args.size
	c.editdistance=args.editdistance
	c.overlapping=args.overlapping
	c.precisemotif=args.precisemotif
	c.precisetimes=args.precisetimes
	c.coverage=args.coverage
	c.softclipping=args.softclipping
	c.samplename=args.samplename
	c.threads=args.threads
	c.mmidir=args.mmidir
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
			
		elif not args.index_only and os.listdir(os.path.abspath(c.OUT)):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] The output folder is not empty: specify another output folder or clean the current one')
			sys.exit(1)

	#check if bedtools is in PATH. As from v1.1, TRiCoLOR does not perform subprocess calls to bedtools but this bedtools is still required in PATH

	if which('bedtools') is None:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] bedtools must be in PATH')
		sys.exit(1)

	try:

		bedfile=pybedtools.BedTool(c.BED)
		bedsrtd=bedfile.sort()

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Errror] BED does not exist, is not readable or is not a valid BED')
		sys.exit(1)


	try:

		ref=pyfaidx.Fasta(c.REF)

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Errror] Reference file does not exist, is not readable or is not a valid FASTA')
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

	if args.index_only:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Only creating mappy .mmi indexes for chromosomes in BED')
		b_chroms={x.chrom for x in bedsrtd if x.chrom not in skip_c}
	
		for b_chrom in b_chroms:

			chromind=os.path.abspath(c.OUT + '/' + b_chrom + '.mmi')

			if not os.path.exists(chromind):

				try:
					
					chrom=ref[b_chrom] #extract pyfaidx chrom

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Missing ' + b_chrom + ' in reference file')
					sys.exit(1)

				refseq=chrom[:len(chrom)].seq #get chromosome reference sequence
				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Message] Creating .mmi index for chromosome ' + b_chrom)

				#removed subprocess call to samtools
				if not os.path.exists(os.path.abspath(c.OUT + '/' + b_chrom + '.fa')): #create .fa first if this does not exist

					with open(os.path.abspath(c.OUT + '/' + b_chrom + '.fa'), 'w') as chromout:

						chromout.write('>' + b_chrom + '\n' + refseq)

				#removed subprocess call to minimap2
				mp.Aligner(os.path.abspath(c.OUT + '/' + b_chrom + '.fa'),fn_idx_out=chromind, preset='asm10') #standard present for consensus-to-ref alignment
				os.remove(os.path.abspath(c.OUT + '/' + b_chrom + '.fa'))

			else:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Warning] .mmi index for chromosome ' + b_chrom + ' already exists: skipped')

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Done')
		sys.exit(0)

	#this is only if we are calling repetitions
	os.makedirs(os.path.abspath(c.OUT + '/reference'))
	os.makedirs(os.path.abspath(c.OUT + '/haplotype1'))
	os.makedirs(os.path.abspath(c.OUT + '/haplotype2'))

	c.BAM=[os.path.abspath(x) for x in args.bamfile[0]]

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

	c.regex=finder.RegexBuilder(c) #regex builder
	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Using RegEx ' + c.regex)

	if c.threads > multiprocessing.cpu_count():

		c.threads=multiprocessing.cpu_count()-1
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] Specified number of cores exceeds the number of cores available. Using all but one')

	if c.mmidir is None:

		gendir=os.path.dirname(c.REF)

	else:

		gendir=os.path.abspath(c.mmidir)

		if not os.path.exists(gendir):

			os.makedirs(gendir)

		if not os.access(gendir, os.W_OK): #test write permissions one for all

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Missing write permissions on the folder chosen to store .mmi indexes')
			sys.exit(1)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	b_chroms={x.chrom for x in bedsrtd if x.chrom not in skip_c} #all the chromosomes in BED except those to exclude
	#sort b_chroms by the chromosome order bcftools expect (that is, natural order)
	b_chroms=sorted(b_chroms,key=natural_keys)
	
	#initialize multiprocessing instances
	manager = multiprocessing.Manager()
	Rrep=manager.list()
	H1rep=manager.list()
	H2rep=manager.list()
	H1bam=manager.list()
	H2bam=manager.list()
	SortedVCF=list()

	BEDheader='#Chromosome\tStart\tEnd\tRepeated Motif\tRepetitions Number\n'
	rBEDstring=''
	h1BEDstring=''
	h2BEDstring=''

	#write final VCF header, store BAM header as well
	bamh=writer.VCF_headerwriter(c) #write header (unique for all the VCF files, will append variant lines in the end)

	for b_chrom in b_chroms:

		VCFvariants=manager.list()

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Processing chromosome ' + b_chrom)

		try:
			
			chrom=ref[b_chrom] #extract pyfaidx chrom

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Missing ' + b_chrom + ' in reference file')
			sys.exit(1)

		c.refseq=chrom[:len(chrom)].seq #get chromosome reference sequence
		query=pybedtools.Interval(b_chrom,0,len(chrom)) #this is the BED-based interval for the entire chromosome
		chromind=os.path.abspath(gendir + '/' + b_chrom + '.mmi')

		if not os.path.exists(chromind):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Creating .mmi index for current chromosome')

			#removed subprocess call to samtools
			if not os.path.exists(os.path.abspath(gendir + '/' + b_chrom + '.fa')): #create .fa first if this does not exist

				with open(os.path.abspath(gendir + '/' + b_chrom + '.fa'), 'w') as chromout:

					chromout.write('>' + b_chrom + '\n' + c.refseq)

			#removed subprocess call to minimap2
			c.aligner=mp.Aligner(os.path.abspath(gendir + '/' + b_chrom + '.fa'),fn_idx_out=chromind, preset='asm10') #standard present for consensus-to-ref alignment

		else:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Loading .mmi index for current chromosome')
			c.aligner=mp.Aligner(chromind, preset='asm10') #use existent


		#extract all the regions in BED for current chromosome and prepare for multi-processing
		p_reg=bedsrtd.all_hits(query)
		reg=len(p_reg)		
		chunk_size=reg/c.threads
		slices=Chunks(p_reg,math.ceil(chunk_size))
		processes = []

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Calling repetitions')

		for i,sli in enumerate(slices):

			processor='p'+str(i+1)
			p=multiprocessing.Process(target=Runner, args=(sli,processor,c,Rrep,H1rep,H2rep,H1bam,H2bam,VCFvariants))
			p.start()
			processes.append(p)

		for p in processes:
		
			p.join()

		#sort VCF entries for current chromosome

		SortedVCF.extend((sorted(VCFvariants, key=itemgetter(1))))

	#write BED files. These are already sorted
	#write REF BED
	#BED has start 0-based and end 1-based

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Writing reference repetitions to BED')

	for r_rep in Rrep:

		if r_rep != []:

			rBEDstring+=r_rep[0][0] + '\t' + str(r_rep[0][2]) + '\t' + str(r_rep[0][3]+1) + '\t' + r_rep[0][1] + '\t' + str(r_rep[0][4]) + '\n'

	rBEDfile=pybedtools.BedTool(rBEDstring,from_string=True)
	rBEDfile_s=rBEDfile.sort()
	rBEDfile_s.saveas(os.path.abspath(c.OUT + '/reference/TRiCoLOR.srt.bed.gz'), compressed=True, trackline=BEDheader)

	#write H1 BED

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Writing haplotype1 repetitions to BED')

	for h1_rep in H1rep:

		if h1_rep != []:

			h1BEDstring+=h1_rep[0][0] + '\t' + str(h1_rep[0][2]) + '\t' + str(h1_rep[0][3]+1) + '\t' + h1_rep[0][1] + '\t' + str(h1_rep[0][4]) + '\n'

	h1BEDfile=pybedtools.BedTool(h1BEDstring,from_string=True)
	h1BEDfile_s=h1BEDfile.sort()
	h1BEDfile_s.saveas(os.path.abspath(c.OUT + '/haplotype1/TRiCoLOR.srt.bed.gz'), compressed=True, trackline=BEDheader)
 
	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Writing haplotype2 repetitions to BED')

	#write H2 BED

	for h2_rep in H2rep:

		if h2_rep != []:

			h2BEDstring+=h2_rep[0][0] + '\t' + str(h2_rep[0][2]) + '\t' + str(h2_rep[0][3]+1) + '\t' + h2_rep[0][1] + '\t' + str(h2_rep[0][4]) + '\n'

	h2BEDfile=pybedtools.BedTool(h2BEDstring,from_string=True)
	h2BEDfile_s=h2BEDfile.sort()
	h2BEDfile_s.saveas(os.path.abspath(c.OUT + '/haplotype2/TRiCoLOR.srt.bed.gz'), compressed=True, trackline=BEDheader)

	#write H1 bam

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Writing haplotype1 consensus BAM')

	BamW(bamh,H1bam,os.path.abspath(c.OUT + '/haplotype1/TRiCoLOR.bam'))
	pysam.sort('-o', os.path.abspath(c.OUT + '/haplotype1/TRiCoLOR.srt.bam'), '-@', str(c.threads), os.path.abspath(c.OUT + '/haplotype1/TRiCoLOR.bam'))
	pysam.index(os.path.abspath(c.OUT + '/haplotype1/TRiCoLOR.srt.bam'))
	os.remove(os.path.abspath(c.OUT + '/haplotype1/TRiCoLOR.bam'))

	#write H2 bam

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Writing haplotype2 consensus BAM')

	BamW(bamh,H2bam,os.path.abspath(c.OUT + '/haplotype2/TRiCoLOR.bam'))
	pysam.sort('-o', os.path.abspath(c.OUT + '/haplotype2/TRiCoLOR.srt.bam'), '-@', str(c.threads), os.path.abspath(c.OUT + '/haplotype2/TRiCoLOR.bam'))
	pysam.index(os.path.abspath(c.OUT + '/haplotype2/TRiCoLOR.srt.bam'))
	os.remove(os.path.abspath(c.OUT + '/haplotype2/TRiCoLOR.bam'))

	#write VCF variants
	
	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Writing variants to VCF')
	
	for v in SortedVCF: #does this kind of sort works? this is basically pybedtools sort for chromosomes

		with open(os.path.abspath(c.OUT + '/TRiCoLOR.srt.vcf'), 'a') as vcfout:

			vcfout.write(v[0])

	pysam.tabix_index(os.path.abspath(c.OUT + '/TRiCoLOR.srt.vcf'), preset='vcf')

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Done')
	sys.exit(0)