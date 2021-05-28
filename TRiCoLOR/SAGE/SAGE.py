#!/usr/bin/python3 env

#python 3 standard library

import sys
import os
import re
import subprocess
import math
import statistics
import itertools
import multiprocessing
from collections import defaultdict
from datetime import datetime,date

# additional modules

import pysam
import numpy as np
import editdistance
from cyvcf2 import VCF


class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	BAM = list()
	OUT = ''
	VCF= ''
	match=0
	mismatch=0
	gapopen=0
	gapextend=0
	coverage=0
	samplename=[]
	threads=0
	mendel=False

	#add other stuff that will be used later

	names=[]
	Cpath=os.path.abspath(os.path.dirname(os.path.dirname(__file__)) + '/REFER/consensus') #this is the default location. Should not be moved


def GenotypeCombos():

	'''
	Get possible combos
	'''

	alleles=['0', '1', '2', '.']
	parent1=parent2=[x+'|'+y for x,y in itertools.product(alleles, alleles)]
	combos=set([(x,y) for x,y in itertools.product(parent1,parent2)])

	return combos


def GetPossibleGenotypes(genchild,combos):

	'''
	Get possible genotypes
	'''

	correct=set()

	if genchild=='.|.':

		correct=combos

	else:

		if '.' in genchild:

			if genchild[0] == '.':

				for combo in combos:

					parent1=combo[0]
					parent2=combo[1]

					if genchild[-1] in parent1 or '.' in parent1 or genchild[-1] in parent2 or '.' in parent2:

						correct.add(combo)

			else:

				for combo in combos:

					parent1=combo[0]
					parent2=combo[1]

					if genchild[0] in parent1 or '.' in parent1 or genchild[0] in parent2 or '.' in parent2:

						correct.add(combo)

		else:

			for combo in combos:

				parent1=combo[0]
				parent2=combo[1]

				if genchild[0] in parent1 or '.' in parent1:

					if genchild[-1] in parent2 or '.' in parent2:

						correct.add(combo)

					else:

						if genchild[-1] in parent1 or '.' in parent1:

							if genchild[0] in parent2 or '.' in parent2:

								correct.add(combo)

				else:

					if genchild[-1] in parent1 or '.' in parent1:

						if genchild[0] in parent2 or '.' in parent2:

							correct.add(combo)

	return correct


def GenotypeDict(combos):

	'''
	Create a genotype dict for child
	'''

	childdict=dict()
	genchildtype=['1|0', '0|1', '1|1', '1|2', '1|.', '.|1', '0|.', '.|0', '.|.']

	for key in genchildtype:

		childdict[key] = GetPossibleGenotypes(key, combos)

	return childdict


def CheckMendelian(genchild,genparent1,genparent2,childdict):

	'''
	Check consistency
	'''


	if genchild not in childdict.keys(): #should not happen

		return '.'

	elif (genparent1,genparent2) in childdict[genchild]:

		return '0'

	else:

		return '1'


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


def VCF_HeaderModifier(rawheader,c):

	'''
	Modify header of the VCF file from REFER, adding infos from SAGE and new samplenames
	'''

	headlist=rawheader.split('\n')[:-1]
	newheader=''

	for el in headlist:

		if el.startswith('##bcftools') or el.startswith('##INFO=<ID=H') or el.startswith('##INFO=<ID=MAPQ'):

			continue

		elif el.startswith('##filedate'):

			newheader+= '##filedate=' + ''.join(str(date.today()).split('-')) +  '\n'

		elif el.startswith('##source'):

			newheader+='##source=TRiCoLOR\n'

		elif el.startswith('##FORMAT=<ID=DP2'):

			newheader+=el + '\n' + '##FORMAT=<ID=GS,Number=1,Type=Float,Description="Genotype Score (>=0 and <=2)">' + '\n'

		elif el.startswith('##SAMPLE'):

			newheader+=el + '\n'

			for name in c.names:

				newheader+= '##SAMPLE=<ID=' + name + '>' + '\n'

		elif el.startswith('#CHROM'):

			for name in c.names:

				el+='\t' + name.upper()

			newheader+=el + '\n'

		elif el.startswith('##INFO=<ID=AED'):

			newheader += el + '\n' + '##INFO=<ID=MISSR,Number=1,Type=Float,Description="Missing genotypes (ratio)">' +'\n' + '##INFO=<ID=MENDEL,Number=1,Type=Integer,Description="Mendelian consistency [.=Unknown;0=Consistent;1=Inconsistent]">' + '\n'

		else:

			newheader+=el + '\n'

		with open(os.path.abspath(c.OUT + '/TRiCoLOR.srt.vcf'), 'w') as vcfout:

			vcfout.write(newheader)


#the following 4 functions are inherited from REFER, slightly modified


def sub_none(coordinates):

	'''
	Substitute None (soft-clipped/inserted) coordinates with a negative value
	'''

	return [-9999999 if v is None else v for v in coordinates]


def find_nearest(array, value):

	'''
	Find closest match in array, given value
	'''

	return (np.abs(array - value)).argmin()


def Bamfile_Analyzer(bamfilein,chromosome,start,end,c,out,processor):

	'''
	Parse haplotype-resolved BAM and store sequences in FASTA for consensus calculation
	'''

	cov=0
	headers=[]
	sequences=[]
	lengths=[]

	bamfile=pysam.AlignmentFile(bamfilein,'rb')	

	for read in bamfile.fetch(chromosome,start,end):

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			read_start=read.reference_start
			read_end=read.reference_end

			if read_start <= start and read_end >= end:

				sequence=read.query_sequence #changed seq with query sequence
				coord=np.asarray(sub_none(read.get_reference_positions(full_length=True)))
				header=read.query_name

				#s_,e_=min(coord, key=lambda x:abs(x-start)), min(coord, key=lambda x:abs(x-end)) 
				#s_i,e_i = [i for i,e in enumerate(coord) if e == s_][0], [i for i,e in enumerate(coord) if e == e_][0]
				#finalseq=sequence[s_i:e_i+1]
				si,ei=find_nearest(coord,start),find_nearest(coord,end)
				finalseq=sequence[si:ei]

				headers.append(header)
				sequences.append(finalseq)
				lengths.append(len(finalseq))


	bamfile.close()

	if len(sequences) >= c.coverage:

		averagelen=statistics.mean(lengths)
		stdevlen=statistics.stdev(lengths)

		minbound=averagelen-3*stdevlen
		maxbound=averagelen+3*stdevlen

		fasta=''

		for head,seq,leng in zip(headers,sequences,lengths):

			if leng <= minbound or leng >= maxbound:

				continue #exclude outliers

			else:

				cov+=1
				
				if cov > 100: #arbitrary treshold to avoid problems with coverage too high
					
					break #do not add other sequences, as it slows down the subsequent consensus computation
					
				else:
				
					fasta+='>' + head + '\n' + seq + '\n'

		if cov >= c.coverage:

			with open(os.path.abspath(out+'/' + processor + '.unaligned.fa'),'w') as fastaout:

				fastaout.write(fasta)

	else:

		cov=len(sequences)

	return cov


def Bamfile_Analyzer_Single(bamfilein,chromosome,start,end,c,out,processor):

	'''
	Parse haplotype-tagged BAM and store sequences in FASTA for consensus calculation
	'''

	cov1=0
	cov2=0
	headers1=[]
	sequences1=[]
	lengths1=[]
	headers2=[]
	sequences2=[]
	lengths2=[]

	bamfile=pysam.AlignmentFile(bamfilein,'rb')	

	for read in bamfile.fetch(chromosome,start,end):

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			if not read.has_tag('HP'):

				continue

			else:


				read_start=read.reference_start
				read_end=read.reference_end

				if read_start <= start and read_end >= end:

					sequence=read.query_sequence #changed seq with query sequence
					coord=np.asarray(sub_none(read.get_reference_positions(full_length=True)))
					header=read.query_name

					#s_,e_=min(coord, key=lambda x:abs(x-start)), min(coord, key=lambda x:abs(x-end)) 
					#s_i,e_i = [i for i,e in enumerate(coord) if e == s_][0], [i for i,e in enumerate(coord) if e == e_][0]
					#finalseq=sequence[s_i:e_i+1]
					si,ei=find_nearest(coord,start),find_nearest(coord,end)
					finalseq=sequence[si:ei]

					if read.get_tag('HP') == 1: 

						headers1.append(header)
						sequences1.append(finalseq)
						lengths1.append(len(finalseq))

					else:

						headers2.append(header)
						sequences2.append(finalseq)
						lengths2.append(len(finalseq))

	bamfile.close()

	#hap1

	if len(sequences1) >= c.coverage:

		averagelen=statistics.mean(lengths1)
		stdevlen=statistics.stdev(lengths1)

		minbound=averagelen-3*stdevlen
		maxbound=averagelen+3*stdevlen

		fasta=''

		for head,seq,leng in zip(headers1,sequences1,lengths1):

			if leng <= minbound or leng >= maxbound:

				continue #exclude outliers

			else:

				cov1+=1
				
				if cov1 > 100: #arbitrary treshold to avoid problems with coverage too high
					
					break #do not add other sequences, as it slows down the subsequent consensus computation
					
				else:
				
					fasta+='>' + head + '\n' + seq + '\n'

		if cov1 >= c.coverage:

			with open(os.path.abspath(out+'/haplotype1/' + processor + '.unaligned.fa'),'w') as fastaout:

				fastaout.write(fasta)

	else:

		cov1=len(sequences1)


	#hap2

	if len(sequences2) >= c.coverage:

		averagelen=statistics.mean(lengths2)
		stdevlen=statistics.stdev(lengths2)

		minbound=averagelen-3*stdevlen
		maxbound=averagelen+3*stdevlen

		fasta=''

		for head,seq,leng in zip(headers2,sequences2,lengths2):

			if leng <= minbound or leng >= maxbound:

				continue #exclude outliers

			else:

				cov2+=1
				
				if cov2 > 100: #arbitrary treshold to avoid problems with coverage too high
					
					break #do not add other sequences, as it slows down the subsequent consensus computation
					
				else:
				
					fasta+='>' + head + '\n' + seq + '\n'

		if cov2 >= c.coverage:

			with open(os.path.abspath(out+'/haplotype2/' + processor + '.unaligned.fa'),'w') as fastaout:

				fastaout.write(fasta)

	else:

		cov2=len(sequences2)


	return cov1,cov2


def GetGTandGS(test,ref,alts):

	'''
	Most likely genotype and edit distance-based quality
	'''


	Redit=editdistance.eval(ref,test) #edit distance between test and reference 
	Rscore=Redit/len(max([test,ref], key=len)) #calculate an edit-distance score
	iRscore=1.0-Rscore #inverted score. max 1 (ref==test), min 0 if difference is max

	altscores=[]

	for alt in alts:

		Aedit=editdistance.eval(alt,test) #edit distance between test and alteration
		Ascore=Aedit/len(max([test,alt], key=len)) #calculate an edit-distance score
		iAscore=1.0-Ascore #inverted score. max 1 (alt==test), min 0 if difference is max
		altscores.append(iAscore)

	iAlt=altscores.index(max(altscores))

	if altscores[iAlt] > iRscore: #choose iAlt

		return str(iAlt +1),round(float(altscores[iAlt]),2)

	elif altscores[iAlt] < iRscore: #choose iRscore

		return '0', round(float(iRscore),2)

	else: #if the same, do not choose

		return '.', float(0)


def Runner(sli,name,bam_s,c,processor,Slist):

	'''
	Run genotyping in parallel
	'''
	out1=os.path.abspath(c.OUT + '/' + name +'/haplotype1')
	out2=os.path.abspath(c.OUT + '/' + name +'/haplotype2')
	out=os.path.dirname(out1)

	for s in sli:

		chromosome,start,end,ref,alt,genchild,raed,aed,dp1,dp2,length=s

		if type(bam_s) == str: #hp-tagged BAM

			cov_inbam1,cov_inbam2=Bamfile_Analyzer_Single(bam_s,chromosome,start,end,c,out,processor)

		else: #	

			cov_inbam1=Bamfile_Analyzer(bam_s[0],chromosome,start,end,c,out1,processor)
			cov_inbam2=Bamfile_Analyzer(bam_s[1],chromosome,start,end,c,out2,processor)

		
		file1=os.path.abspath(out1 + '/' + processor + '.unaligned.fa')
		file2=os.path.abspath(out2 + '/' + processor + '.unaligned.fa')

		if os.path.exists(file1):

			#generate consensus and skip re-mapping
			consfa1=os.path.abspath(out1 + '/' + processor + '.consensus.fa')
			
			with open(consfa1, 'w') as fout:

				subprocess.call([c.Cpath, str(c.match), str(c.mismatch), str(c.gapopen), str(c.gapextend), file1], stdout=fout, stderr=open(os.devnull, 'wb'))

			seq1=''

			with open(consfa1) as fin:

				for line in fin:

					if line[0] != '>':

						seq1+=line.rstrip()

			os.remove(file1)
			os.remove(consfa1)

		else:

			seq1=''

		if os.path.exists(file2):

			#generate consensus and skip re-mapping
			consfa2=os.path.abspath(out2 + '/' + processor + '.consensus.fa')
			
			with open(consfa2, 'w') as fout:

				subprocess.call([c.Cpath, str(c.match), str(c.mismatch), str(c.gapopen), str(c.gapextend), file2], stdout=fout, stderr=open(os.devnull, 'wb'))

			seq2=''

			with open(consfa2) as fin:

				for line in fin:

					if line[0] != '>':

						seq2+=line.rstrip()
			
			os.remove(file2)
			os.remove(consfa2)

		else:

			seq2=''
		
		if seq1 == '' and seq2 == '':

			genotype='.|.'
			quality= float(0)

		elif seq1 != '' and seq2 == '':

			gen1,qual1=GetGTandGS(seq1,ref,alt)
			genotype=gen1 + '|.'

			if qual1 == '.':

				quality=float(0)

			else:

				quality=qual1

		elif seq1 == '' and seq2 != '':

			gen2,qual2=GetGTandGS(seq2,ref,alt)
			genotype='.|' + gen2

			if qual2 == '.':

				quality=float(0)

			else:

				quality=qual2

		else:

			gen1,qual1=GetGTandGS(seq1,ref,alt)
			gen2,qual2=GetGTandGS(seq2,ref,alt)

			genotype= gen1 + '|' + gen2

			if qual1 == '.':

				if qual2 == '.':

					quality = float(0)

				else:

					quality=qual2

			else:

				if qual2 == '.':

					quality=qual1

				else:

					quality=round(float(qual1+qual2),2)

		Slist.append((chromosome,start,end,ref,alt,raed,aed,genchild,dp1,dp2,genotype,quality,seq1,seq2,cov_inbam1,cov_inbam2,length))


def run(parser, args):

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] TRiCoLOR SAGE v1.1')

	'''
	Check arguments, run functions
	'''

	#fill container

	c.BAM=args.bamfile
	c.OUT=os.path.abspath(args.output)
	c.VCF=os.path.abspath(args.vcffile)
	c.match=args.match
	c.mismatch=args.mismatch
	c.gapopen=args.gapopen
	c.gapextend=args.gapextend
	c.coverage=args.coverage
	c.samplename=args.samplename
	c.threads=args.threads
	c.mendel=args.mendel

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

	if len(c.BAM) != 2:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] SAGE strictly requires a couple of HP-tagged BAM or a comma-separated couple of haplotype-spliited BAM')
		sys.exit(1)


	for bam in c.BAM:

		if type(bam) == str: #HP-tagged couple

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

		else: #haplotype-resolved BAM

			for x in bam:

				try:

					pysam.quickcheck(x)

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Errror] BAM ' + x + ' does not exist, is not readable or is not a valid BAM')
					sys.exit(1)

				if not os.path.exists(x + '.bai'):

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Warning] Missing ' + x + ' index. Creating')

					try:

						pysam.index(x)

					except:

						now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
						print('[' + now + '][Errror] BAM ' + x + ' could not be indexed')
						sys.exit(1)	

	if c.samplename is None:

		c.names.append('SAMPLE1')
		c.names.append('SAMPLE2')

	else:

		if len(c.samplename[0]) != 2:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Errror] If provided, one samplename for each sample is required')
			sys.exit(1)

		else:

			c.names=c.samplename[0]

	if c.threads > multiprocessing.cpu_count():

		c.threads=multiprocessing.cpu_count()-1
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] Specified number of cores exceeds the number of cores available. Using all but one')

	try:

		it=VCF(c.VCF)

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Errror] VCF does not exist, is not readable or is not a valid VCF')
		sys.exit(1)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Processing input VCF')

	vheader=it.raw_header
	vcfall=[]

	for var in it:

		if var.ALT != []:

			vcfall.append((var.CHROM, var.start+1, var.INFO.get('TREND'), var.REF, var.ALT, var.genotypes[0], var.INFO.get('RAED'), var.INFO.get('AED'), int(var.format('DP1')), int(var.format('DP2')), var.INFO.get('TRLEN')))

	it.close()

	VCF_HeaderModifier(vheader,c)

	#preparing for multiprocessing

	entries=len(vcfall)
	chunk_size=entries/c.threads
	slices=Chunks(vcfall,math.ceil(chunk_size))
	manager = multiprocessing.Manager()
	Namesdict=dict() #results from parents

	for name,bam_s in zip(c.names,c.BAM):

		os.makedirs(os.path.abspath(c.OUT + '/' + name +'/haplotype1'))
		os.makedirs(os.path.abspath(c.OUT + '/' + name +'/haplotype2'))
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Genotyping ' + name)

		processes=[]
		Slist=manager.list()

		for i,sli in enumerate(slices):

			processor='p'+str(i+1)
			p=multiprocessing.Process(target=Runner, args=(sli,name,bam_s,c,processor,Slist))
			p.start()
			processes.append(p)

		for p in processes:

			p.join()

		sorted_Slist=sorted(Slist, key=lambda x: (natural_keys(x[0]),x[1])) #sort as required by bcftools
		Namesdict[name]=sorted_Slist

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Writing variants to VCF')
	
	FILTER='.'
	QUAL='.'
	ID='.'
	FORMAT='GT:DP1:DP2:GS'
	QUALCHILD=2.0

	if c.mendel:

		combos=GenotypeCombos()
		childdict=GenotypeDict(combos)

	with open(os.path.abspath(c.OUT + '/TRiCoLOR.srt.vcf'), 'a') as vcfout:

		for i in range(len(vcfall)):

			CHROM,POS,END,REF,ALT,RAED,AED,GENCHILD,DP1,DP2,LENGTH=str(Namesdict[c.names[0]][i][0]),str(Namesdict[c.names[0]][i][1]), str(Namesdict[c.names[0]][i][2]), Namesdict[c.names[0]][i][3], ','.join(x for x in Namesdict[c.names[0]][i][4]), str(Namesdict[c.names[0]][i][5]), str(Namesdict[c.names[0]][i][6]),'|'.join(str(x) for x in Namesdict[c.names[0]][i][7][:-1]).replace('-1', '.'), str(Namesdict[c.names[0]][i][8]), str(Namesdict[c.names[0]][i][9]),str(Namesdict[c.names[0]][i][16])

			toadd=[]
			seqs=[]
			missing=0

			for names in c.names:

				toadd.append(Namesdict[names][i][10] + ':' + str(Namesdict[names][i][14])  + ':' + str(Namesdict[names][i][15]) + ':' + repr(Namesdict[names][i][11]))
				seqs.append(names + ':' + Namesdict[names][i][12] + ',' + Namesdict[names][i][13])
				missing += Namesdict[names][i][10].count('.')

			MISSR=missing/((len(c.names)+1)*2)
			
			if c.mendel:

				MENDEL=CheckMendelian(GENCHILD, toadd[0].split(':')[0], toadd[1].split(':')[0],childdict)

			else:

				MENDEL='.'

			vcfout.write(CHROM + '\t' + POS + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER + '\t' + 'TREND='+END + ';TRLEN='+ LENGTH + ';RAED=' +RAED + ';AED=' + AED + ';MISSR=' + str(MISSR) + ';MENDEL=' + str(MENDEL) + '\t' + FORMAT + '\t' + GENCHILD + ':'+ DP1 + ':' + DP2 + ':' + str(QUALCHILD) + '\t' + '\t'.join(x for x in toadd) + '\n')

	pysam.tabix_index(os.path.abspath(c.OUT + '/TRiCoLOR.srt.vcf'), preset='vcf')

	#final clean-up
	for names in c.names:

		os.rmdir(os.path.abspath(c.OUT + '/' + names + '/haplotype1'))
		os.rmdir(os.path.abspath(c.OUT + '/' + names + '/haplotype2'))
		os.rmdir(os.path.abspath(c.OUT + '/' + names))

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Done')
	sys.exit(0)

