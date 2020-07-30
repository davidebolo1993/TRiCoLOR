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
import statistics
import multiprocessing
import logging
import datetime
from collections import defaultdict, Counter
from operator import itemgetter
from bisect import bisect_left,bisect_right
from shutil import which

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
	

	os.makedirs(os.path.abspath(args.output) + '/reference')
	os.makedirs(os.path.abspath(args.output) + '/haplotype1')
	os.makedirs(os.path.abspath(args.output) + '/haplotype2')
	command_dict= vars(args)	
	notkey=['func']
	command_string= ' '.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey)
	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR.REFER.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')	

	print('Initialized .log file ' + os.path.abspath(args.output + '/TRiCoLOR.REFER.log'))

	logging.info('main=TRiCoLOR ' + command_string)
	external_tools=['minimap2', 'samtools', 'bcftools'] #? ALSO ADD NGMLR. NOT PRIORITY. MINIMAP2 IS FASTER AND MORE ACCURATE

	for tools in external_tools: 

		if which(tools) is None:

			logging.error(tools + ' cannot be executed. Install ' + tools + ' and re-run TRiCoLOR REFER')
			exitonerror()

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>'))

	except:

		logging.error('Reference file does not exist, is not readable or is not a valid FASTA')
		exitonerror()

	bams=args.bamfile[0]

	if len(bams) != 2:

		logging.error('TRiCoLOR supports only diploid individuals')
		exitonerror()

	for bam in bams:

		try:

			subprocess.call(['samtools','quickcheck', os.path.abspath(bam)],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

		except:

			logging.error('BAM ' + bam + ' does not exist, is not readable or is not a valid BAM')
			exitonerror()

		if not os.path.exists(os.path.abspath(bam + '.bai')):

			logging.warning('Missing index for BAM ' + bam + '. Not a BAM file scanned with SENSOR. Creating index ...')

			try:

				subprocess.call(['samtools', 'index', os.path.abspath(bam)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

			except:

				logging.error('BAM ' + bam + ' could not be indexed')
				exitonerror()

	logging.info('Match reward for consensus computation: ' +str(args.match))
	logging.info('Mismatch penalty for consensus computation: '+ str(args.mismatch))
	logging.info('Gap opening penalty for consensus computation: ' + str(args.gapopen))
	logging.info('Gap extending penalty for consensus computation: ' + str(args.gapextend))

	if not args.precisemotif:

		if args.motif == 0 or args.motif == 1:

			logging.info('Repetition motif length: any')

		else:

			logging.info('Repetition motif length: at least ' + str(args.motif))

	else:

		if args.motif == 0:

			logging.error('Cannot search for repetitions with motif length 0')
			exitonerror()

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
	logging.info('Haplotypes: ' + str(len(bams)))
	logging.info('Coverage treshold: ' + str(args.coverage))
	logging.info('Soft-clipping treshold: ' + str(args.softclipping))
	logging.info('Long reads type: ' + str(args.readstype))

	if args.threads > multiprocessing.cpu_count():

		cores=multiprocessing.cpu_count()
		logging.warning(str(args.threads) + ' cores specified but only ' + str(cores) + ' available. Using all cores available')

	else:

		cores=args.threads

	logging.info('Cores: ' + str(cores))

	Cpath=os.path.abspath(os.path.dirname(__file__) + '/consensus') #? FAST ENOUGH
	SHCpath=os.path.abspath(os.path.dirname(__file__) + '/consensus.sh')
	SHMpath=os.path.abspath(os.path.dirname(__file__) + '/merging.sh')	
	ref=pyfaidx.Fasta(os.path.abspath(args.genome))
	
	if args.mmidir is None:

		gendir=os.path.abspath(os.path.dirname(args.genome))

	else:

		gendir=os.path.abspath(args.mmidir)

		if not os.path.exists(gendir):

			os.mkdir(gendir)

	b_in=Bed_Reader(os.path.abspath(args.bedfile))
	b_chroms=sorted(set(iter(b_in)), key=natural_keys)

	for b_chrom in b_chroms:

		logging.info('Parsing input BED for chromosome ' + b_chrom + ' ...')
		print('Parsing input BED for chromosome ' + b_chrom + ' ...')

		b_in2=Bed_Reader(os.path.abspath(args.bedfile), chromosome=b_chrom)
		p_reg=list(iter(b_in2))
	
		logging.info('Finding repetitions on chromosome ' + b_chrom + ' ...')
		print('Finding repetitions on chromosome ' + b_chrom + ' ...')

		chrom=ref[b_chrom]
		refseq=chrom[:len(chrom)].seq

		if args.readstype == 'ONT':

			chromind= os.path.abspath(gendir + '/' + b_chrom + '.ont.mmi')

		else:

			chromind= os.path.abspath(gendir + '/' + b_chrom + '.pb.mmi')

		if not os.path.exists(chromind):

			if not os.access(gendir, os.W_OK):

				logging.error('Missing write permissions on the reference genome folder. This is necessary to store chromosomes indexes for re-alignments')
				exitonerror()

			logging.info('Creating .mmi index for ' + b_chrom)

			if not os.path.exists(os.path.abspath(gendir + '/' + b_chrom + '.fa')):

				try:

					with open(os.path.abspath(gendir + '/' + b_chrom + '.fa'), 'w') as chromout:

						subprocess.call(['samtools','faidx', os.path.abspath(args.genome), b_chrom],stderr=open(os.devnull, 'wb'),stdout=chromout)

				except:

					logging.error('Unexpected error while creating index for current chromosome. Does your reference genome contain informations for this chromosome?')
					exitonerror()

			if args.readstype == 'ONT':

				subprocess.call(['minimap2','-d', os.path.abspath(gendir + '/' + b_chrom + '.ont.mmi'),os.path.abspath(gendir + '/' + b_chrom + '.fa')],stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

			else:

				subprocess.call(['minimap2','-H', '-d', os.path.abspath(gendir + '/' + b_chrom + '.pb.mmi'),os.path.abspath(gendir + '/' + b_chrom + '.fa')],stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

		reg=len(p_reg)
		chunk_size=reg/cores
		slices=Chunks(p_reg,math.ceil(chunk_size))
		manager = multiprocessing.Manager()
		Rrep=manager.dict()
		H1rep=manager.dict()
		H2rep=manager.dict()
		bamfile1=bams[0]
		bamfile2=bams[1]

		processes = []

		for i,sli in enumerate(slices):

			processor='p'+str(i+1)

			if not os.path.exists(os.path.abspath(args.output + '/' + processor + '.TRiCoLOR.vcf')):

				writer.VCF_headerwriter(os.path.abspath(bamfile1), os.path.abspath(bamfile2), args.samplename, ','.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey), os.path.abspath(args.output), processor)

			p=multiprocessing.Process(target=Runner, args=(processor,sli,refseq,regex,args.maxmotif,args.size,bamfile1,bamfile2,args.coverage,args.editdistance,chromind,args.readstype,os.path.abspath(args.output),Rrep,H1rep,H2rep,Cpath,SHCpath,args.match,args.mismatch,args.gapopen,args.gapextend,args.softclipping))
			p.start()
			processes.append(p)
		
		for p in processes:
		
			p.join()

		for key in sorted(Rrep.keys(), key=natural_keys):

			writer.BED_repswriter(b_chrom,Rrep[key],os.path.abspath(args.output + '/reference'))
			writer.BED_repswriter(b_chrom,H1rep[key],os.path.abspath(args.output + '/haplotype1'))
			writer.BED_repswriter(b_chrom,H2rep[key],os.path.abspath(args.output + '/haplotype2'))

		CleanResults(SHMpath, b_chrom, os.path.abspath(args.output), os.path.abspath(bamfile1), os.path.abspath(bamfile2), cores)

		logging.info('Processed chromosome ' + b_chrom)
		print('Processed chromosome ' + b_chrom)

	bcfs=[]

	for i in range(len(glob.glob(os.path.abspath(args.output) + '/*.vcf'))):

		key='p' + str(i+1)
		bcfs.append(os.path.abspath(args.output + '/' + key +  '.TRiCoLOR.bcf'))
		subprocess.call(['bcftools', 'sort', '-o', os.path.abspath(args.output + '/' + key +  '.TRiCoLOR.bcf'), '-O', 'b', os.path.abspath(args.output + '/' + key + '.TRiCoLOR.vcf')],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
		subprocess.call(['bcftools', 'index', os.path.abspath(args.output + '/' + key +  '.TRiCoLOR.bcf')],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
		os.remove(os.path.abspath(args.output + '/' + key + '.TRiCoLOR.vcf'))

	with open(os.path.abspath(args.output + '/TRiCoLOR.bcf.txt'), 'w') as filebcfout:

		for bcf in bcfs:

			filebcfout.write(bcf + '\n')

	subprocess.call(['bcftools', 'concat', '-f', os.path.abspath(args.output + '/TRiCoLOR.bcf.txt'), '-a', '-D', '-O', 'b', '-o', os.path.abspath(args.output + '/TRiCoLOR.bcf'), '--threads', str(cores-1)],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
	os.remove(os.path.abspath(args.output + '/TRiCoLOR.bcf.txt'))
	subprocess.call(['bcftools', 'sort', '-o', os.path.abspath(args.output + '/TRiCoLOR.srt.bcf'), '-O', 'b', os.path.abspath(args.output + '/TRiCoLOR.bcf')],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
	os.remove(os.path.abspath(args.output + '/TRiCoLOR.bcf'))
	subprocess.call(['bcftools', 'index', os.path.abspath(args.output + '/TRiCoLOR.srt.bcf')],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

	for bcf in bcfs:

		os.remove(bcf)
		os.remove(bcf + '.csi')

	logging.info('Done')
	print('Done')


##CLASS


class Bed_Reader():


	def __init__(self,bedfile,chromosome=None):

		self.bedfile=bedfile
		self.chromosome=chromosome

	def __iter__(self):

		with open (self.bedfile, 'r') as bedin:

			for line in csv.reader(bedin, delimiter='\t'): #? FIND ANOTHER, FASTER, WAY TO PARSE BED. NOT PRIORITY

				if not line[0].startswith('#') and line !=[]: 

					if len(line) < 3:

						logging.error('BED to TRiCoLOR REFER -bed/--bedfile must contain at least: chromosome, start, end (other fields are ignored)')
						exitonerror()

					else:

						if self.chromosome is None:

							yield (line[0])

						else:

							if line[0] == self.chromosome:

								yield (line[0], int(line[1]), int(line[2]))


##FUNCTIONS


def exitonerror():

	
	print('An error occured. Check the log file for more details')
	sys.exit(1)

	
def atoi(text):


    return int(text) if text.isdigit() else text


def natural_keys(text):


    return [ atoi(c) for c in re.split(r'(\d+)', text)]


def Chunks(l,n):

	return [l[i:i+n] for i in range(0, len(l), n)]


def Runner(processor,sli,refseq,regex,maxmotif,size,bamfile1,bamfile2,coverage,allowed,chromind,readtype,output,Rrep,H1rep,H2rep,Cpath,SHCpath,match,mismatch,gapopen,gapextend,clipping):

	Ritem = Rrep[processor] = list()
	H1item = H1rep[processor] = list()
	H2item = H2rep[processor] =list()
	
	for i,s in enumerate(sli):

		try:

			pR=ReferenceReps(s,refseq,regex,maxmotif,size)

			if pR is not None:

				if pR != []:
		
					Ritem.extend(pR)

				out1=os.path.abspath(output + '/haplotype1')
				pH1,pS1,pC1,pCOV1,pQ1=HaploReps(SHCpath,Cpath, bamfile1,s, coverage, regex, maxmotif, size, allowed, refseq, chromind, out1, processor,i,readtype,match,mismatch,gapopen, gapextend,clipping)

				if pH1 != []:

					H1item.extend(pH1)

				out2=os.path.abspath(output + '/haplotype2')				
				pH2,pS2,pC2,pCOV2,pQ2=HaploReps(SHCpath,Cpath, bamfile2,s, coverage, regex, maxmotif, size, allowed, refseq, chromind, out2, processor,i,readtype,match,mismatch,gapopen,gapextend, clipping)

				if pH2 != []:

					H2item.extend(pH2)

				writer.VCF_writer(s[0], pR, refseq, pH1, pS1,pC1,pCOV1,pQ1,pH2, pS2,pC2,pCOV2,pQ2,output, processor)

				if os.stat(os.path.abspath(out1 + '/' + processor + '.' + str(i+1) + '.srt.bam')).st_size == 0:

					os.remove(os.path.abspath(out1 + '/' + processor + '.' + str(i+1) + '.srt.bam'))
					os.remove(os.path.abspath(out1 + '/' + processor + '.' + str(i+1) + '.srt.bam.bai'))

				if os.stat(os.path.abspath(out2 + '/' + processor + '.' + str(i+1) + '.srt.bam')).st_size == 0:

					os.remove(os.path.abspath(out2 + '/' + processor + '.' + str(i+1) + '.srt.bam'))
					os.remove(os.path.abspath(out2 + '/' + processor + '.' + str(i+1) + '.srt.bam.bai'))

		except BaseException as BE:

			logging.error('Processor ' + processor + ' generated an error for region ' + s[0] +  ':' + str(s[1]) + '-' + str(s[2]) + '. Check TRiCoLOR.' + processor + '.unexpected.err.log in ' + os.path.abspath(output) + ' for details')
			
			with open (os.path.abspath(output + '/TRiCoLOR.' + processor + '.unexpected.err.log'), 'a') as errorout:

				errorout.write('Unexpected error in region '+ s[0] +  ':' + str(s[1]) + '-' + str(s[2]) + '.' + '\n' + str(BE) + '\n')

			exitonerror()


	Rrep[processor] = Ritem
	H1rep[processor] = H1item
	H2rep[processor] = H2item


def ReferenceReps(s,refseq,regex,maxmotif,size):


	chromosome,start,end=s[0],s[1],s[2]

	wanted=refseq[start-1:end]

	if 'N' in wanted:

		return

	else:

		repetitions=list(finder.RepeatsFinder(wanted,regex,maxmotif))
		filtered=finder.ReferenceFilter(repetitions,wanted,size,start)
		
		return filtered


def HaploReps(SHCpath,Cpath, bamfile,s, coverage, regex, maxmotif, size, allowed, refseq, chromind, out, processor,iteration, readtype,match,mismatch,gapopen,gapextend, clipping):


	chromosome,start,end=s[0],s[1],s[2]
	cov_inbam=parser.Bamfile_Analyzer(bamfile,chromosome,start,end,coverage,out,processor)
	file=os.path.abspath(out + '/' + processor + '.unaligned.fa')
	empty=[]

	if not os.path.exists(file):

		open(os.path.abspath(out +'/' + processor + '.' + str(iteration +1) + '.srt.bam'), 'w').close()
		open(os.path.abspath(out +'/' + processor + '.' + str(iteration +1) + '.srt.bam.bai'), 'w').close()

		return empty, empty, empty, empty,empty

	else:

		if readtype == 'ONT':

			mmivar='map-ont'
			#consvar='' #CHECK IF USING DIFFERENT SCORES FOR THE 2 TECHNOLOGIES MAKES SENSE OR NOT

		else:

			mmivar='map-pb'
			#consvar='' #CHECK IF USING DIFFERENT SCORES FOR THE 2 TECHNOLOGIES MAKES SENSE OR NOT

		subprocess.call(['bash', SHCpath, out, Cpath, processor, os.path.basename(file), mmivar, chromind, str(match), str(mismatch), str(gapopen), str(gapextend)],stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'))
		c_bam=os.path.abspath(out + '/' + processor + '.cs.srt.bam')
		coords,seq,qual=finder.Get_Alignment_Positions(c_bam,clipping)

		if seq==[]:

			os.remove(c_bam)
			os.remove(c_bam + '.bai')
			open(os.path.abspath(out +'/' + processor + '.' + str(iteration +1) + '.srt.bam'), 'w').close()
			open(os.path.abspath(out +'/' + processor + '.' + str(iteration +1) + '.srt.bam.bai'), 'w').close()

			return empty, empty, empty, empty, empty

		else:

			repetitions=list(finder.RepeatsFinder(seq,regex,maxmotif))
			os.rename(c_bam,os.path.abspath(out +'/' + processor + '.' + str(iteration +1) + '.srt.bam'))
			os.rename(c_bam + '.bai',os.path.abspath(out +'/' + processor + '.' + str(iteration +1) + '.srt.bam.bai'))
			cor_coord_reps,consensus_string,consensus_coordinates=finder.corrector(refseq, seq,repetitions, coords, size, allowed)

			if set(consensus_coordinates).intersection(list(range(start,end))):

				return cor_coord_reps, consensus_string, consensus_coordinates, cov_inbam, qual

			else:

				return empty, empty, empty, empty, empty



def CleanResults(SHMpath,chromosome, out, bamfile1, bamfile2,cores):


	subprocess.call(['bash', SHMpath, os.path.abspath(out + '/haplotype1'), bamfile1,chromosome,str(cores-1)],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
	subprocess.call(['bash', SHMpath, os.path.abspath(out + '/haplotype2'), bamfile2,chromosome,str(cores-1)],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
