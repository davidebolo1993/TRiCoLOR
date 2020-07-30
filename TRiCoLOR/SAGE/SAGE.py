#!/usr/bin/python env

#python 3 standard library

import sys
import os
import re
import subprocess
import logging
import math
import statistics
import itertools
import multiprocessing
from shutil import which
from bisect import bisect_left,bisect_right
from collections import defaultdict

# additional modules

import pysam
import editdistance
from cyvcf2 import VCF


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

	command_dict= vars(args)
	notkey=['func']
	command_string= ' '.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey)
	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR.SAGE.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	print('Initialized .log file ' + os.path.abspath(args.output + '/TRiCoLOR.SAGE.log'))
	logging.info('main=TRiCoLOR ' + command_string)

	if which('samtools') is None:

		logging.error('samtools cannot be executed. Install samtools and re-run TRiCoLOR SAGE')
		exitonerror()

	for couples in args.bamfile:

		try:

			subprocess.call(['samtools','quickcheck', os.path.abspath(couples[0])],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

		except:

			logging.error(('BAM ' + couples[0] + ' does not exist, is not readable or is not a valid BAM'))
			exitonerror()

		if not os.path.exists(os.path.abspath(couples[0] + '.bai')):

			logging.warning('Missing index for BAM ' + couples[0] + '. Creating index ...')

			try:

				subprocess.call(['samtools', 'index', os.path.abspath(couples[0])], stdout=open(os.devnull,'wb'), stderr=open(os.devnull, 'wb'))

			except:

				logging.error('BAM ' + couples[0] + ' could not be indexed')
				exitonerror()

		try:

			subprocess.call(['samtools','quickcheck', os.path.abspath(couples[1])],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

		except:

			logging.error(('BAM ' + couples[1] + ' does not exist, is not readable or is not a valid BAM'))
			exitonerror()

		if not os.path.exists(os.path.abspath(couples[1] + '.bai')):

			logging.warning('Missing index for BAM ' + couples[1] + '. Creating index ...')

			try:

				subprocess.call(['samtools', 'index', os.path.abspath(couples[1])], stdout=open(os.devnull,'wb'), stderr=open(os.devnull, 'wb'))

			except:

				logging.error('BAM ' + couples[1] + ' could not be indexed')
				exitonerror()

	snames=[]

	if args.samplename is None:

		for i in range(len(args.bamfile)):

			snames.append('SAMPLE' + str(i+1))

	else:

		snames=args.samplename[0]

	if not os.path.exists(os.path.abspath(args.bcffile)):

		logging.error('BCF ' + args.bcffile + ' does not exist')
		exitonerror()

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>'))

	except:

		logging.error('Reference file does not exist, is not readable or is not a valid FASTA')
		exitonerror()

	Cpath=os.path.abspath(os.path.dirname(os.path.dirname(__file__)) + '/REFER/consensus')
	SHCpath=os.path.abspath(os.path.dirname(os.path.dirname(__file__)) + '/REFER/consensus.sh')

	if args.mmidir is None:

		gendir=os.path.abspath(os.path.dirname(args.genome))

	else:

		gendir=os.path.abspath(args.mmidir) #this should exist as this module is run after REFER

	logging.info('Related individuals to genotype: ' + str(len(snames)))

	if args.threads > multiprocessing.cpu_count():

		cores=multiprocessing.cpu_count()
		logging.warning(str(args.threads) + ' cores specified but only ' + str(cores) + ' available. Using all cores available')

	else:

		cores=args.threads

	logging.info('Match reward for consensus computation: ' +str(args.match))
	logging.info('Mismatch penalty for consensus computation: '+ str(args.mismatch))
	logging.info('Gap opening penalty for consensus computation: ' + str(args.gapopen))
	logging.info('Gap extending penalty for consensus computation: ' + str(args.gapextend))
	logging.info('Coverage treshold: ' + str(args.coverage))
	logging.info('Soft-clipping treshold: ' + str(args.softclipping))
	logging.info('Long reads type: ' + str(args.readstype))
	logging.info('Cores: ' + str(cores))

	if args.mendel:

		if len(snames) == 2:

			logging.info('Check mendelian consistency: True')

		else:

			logging.warning('Checking mendelian consistency is possible only when both parents are provided. Switching to False.')

	else:

		logging.info('Check mendelian consistency: False')

	logging.info('Parsing input BCF ...')

	infos,header=ParseBCF(os.path.abspath(args.bcffile))

	logging.info('Done')

	VCF_HeaderModifier(header,snames,os.path.abspath(args.output))

	entries=len(infos)
	chunk_size=entries/cores
	slices=Chunks(infos,math.ceil(chunk_size))

	Namesdict=dict()

	logging.info('Genotyping samples ...')

	for names,couples in zip(snames,args.bamfile):

		logging.info('Genotyping ' + names + ' ...')
		print('Genotyping ' + names + ' ...')

		os.makedirs(os.path.abspath(args.output) + '/' + names + '/haplotype1')
		os.makedirs(os.path.abspath(args.output) + '/' + names + '/haplotype2')

		manager = multiprocessing.Manager()
		PROC_ENTRIES=manager.dict()
		processes=[]

		for i,sli in enumerate(slices):

			processor='p'+str(i+1)
			p=multiprocessing.Process(target=Runner, args=(SHCpath,Cpath,gendir,processor,names,PROC_ENTRIES,sli,couples[0],couples[1],args.coverage,args.match, args.mismatch, args.gapopen, args.gapextend,os.path.abspath(args.output), args.readstype,args.softclipping))
			p.start()
			processes.append(p)

		for p in processes:

			p.join()

		Namesdict[names] = []

		for key in sorted(PROC_ENTRIES.keys(), key=natural_keys):

			Namesdict[names].extend(PROC_ENTRIES[key])

	logging.info('Done')	
	logging.info('Writing to BCF ...')

	FILTER='.'
	QUAL='.'
	ID='.'
	FORMAT='GT:DP1:DP2:GS'
	QUALCHILD=2.0

	if args.mendel:

		combos=GenotypeCombos()
		childdict=GenotypeDict(combos)

	with open(os.path.abspath(args.output + '/TRiCoLOR.vcf'), 'a') as vcfout:

		for i in range(len(infos)):

			CHROM,POS,END,REF,ALT,RAED,AED,GENCHILD,DP1,DP2,LENGTH=Namesdict[snames[0]][i][0], str(Namesdict[snames[0]][i][1]), str(Namesdict[snames[0]][i][2]), Namesdict[snames[0]][i][3], ','.join(x for x in Namesdict[snames[0]][i][4]), str(Namesdict[snames[0]][i][5]), str(Namesdict[snames[0]][i][6]),'|'.join(str(x) for x in Namesdict[snames[0]][i][7][:-1]).replace('-1', '.'), str(Namesdict[snames[0]][i][8]), str(Namesdict[snames[0]][i][9]),str(Namesdict[snames[0]][i][16])

			toadd=[]
			seqs=[]
			missing=0

			for names in snames:

				toadd.append(Namesdict[names][i][10] + ':' + str(Namesdict[names][i][14])  + ':' + str(Namesdict[names][i][15]) + ':' + repr(Namesdict[names][i][11]))
				seqs.append(names + ':' + Namesdict[names][i][12] + ',' + Namesdict[names][i][13])
				missing += Namesdict[names][i][10].count('.')

			MISSR=missing/((len(snames)+1)*2)
			
			if args.mendel:

				MENDEL=CheckMendelian(GENCHILD, toadd[0].split(':')[0], toadd[1].split(':')[0],childdict)

			else:

				MENDEL='.'

			vcfout.write(CHROM + '\t' + POS + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER + '\t' + 'TREND='+END + ';TRLEN='+ LENGTH + ';RAED=' +RAED + ';AED=' + AED + ';MISSR=' + str(MISSR) + ';MENDEL=' + str(MENDEL) + '\t' + FORMAT + '\t' + GENCHILD + ':'+ DP1 + ':' + DP2 + ':' + str(QUALCHILD) + '\t' + '\t'.join(x for x in toadd) + '\n')

			if args.store:

				with open(os.path.abspath(args.output + '/sequences.tsv'), 'a') as relatedout:

					relatedout.write(CHROM + '\t' + POS + '\t' + END + '\t' + '\t'.join(x for x in seqs) + '\n')

	subprocess.call(['bcftools', 'sort', '-o', os.path.abspath(args.output + '/TRiCoLOR.srt.bcf'), '-O', 'b', os.path.abspath(args.output + '/TRiCoLOR.vcf')],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
	subprocess.call(['bcftools', 'index', os.path.abspath(args.output + '/TRiCoLOR.srt.bcf')],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
	os.remove(os.path.abspath(args.output + '/TRiCoLOR.vcf'))
	
	for x in snames:

		os.rmdir(os.path.abspath(args.output + '/' + x + '/haplotype1'))
		os.rmdir(os.path.abspath(args.output + '/' + x + '/haplotype2'))
		os.rmdir(os.path.abspath(args.output + '/' + x ))

	logging.info('Done')
	print('Done')


#FUNCTIONS


def exitonerror():


	print('An error occured. Check the log file for more details')
	sys.exit(1)


def atoi(text):


	return int(text) if text.isdigit() else text


def natural_keys(text):


	return [ atoi(c) for c in re.split(r'(\d+)', text)]


def Chunks(l,n):


	return [l[i:i+n] for i in range(0, len(l), n)]


def sub_none(list_of_coord):


	return [-999999 if v is None else v for v in list_of_coord]


def Bamfile_Analyzer(bamfilein,chromosome,start,end, coverage, out, processor):

	cov=0
	
	start=start-350
	end=end+350

	headers=[]
	sequences=[]
	lengths=[]

	bamfile=pysam.AlignmentFile(bamfilein,'rb')	

	for read in bamfile.fetch(chromosome,start,end):

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			read_start=read.reference_start
			read_end=read.reference_end

			if read_start <= start and read_end >= end:

				sequence=read.seq
				coord=sub_none(read.get_reference_positions(full_length=True))
				header=read.query_name

				s_,e_=min(coord, key=lambda x:abs(x-start)), min(coord, key=lambda x:abs(x-end)) 
				s_i,e_i = [i for i,e in enumerate(coord) if e == s_][0], [i for i,e in enumerate(coord) if e == e_][0]
				finalseq=sequence[s_i:e_i+1]

				headers.append(header)
				sequences.append(finalseq)
				lengths.append(len(finalseq))


	bamfile.close()

	if len(sequences) >= coverage:

		averagelen=statistics.mean(lengths)
		stdevlen=statistics.stdev(lengths)

		minbound=averagelen-3*stdevlen
		maxbound=averagelen+3*stdevlen

		fasta=''

		for head,seq,leng in zip(headers,sequences,lengths):

			if leng <= minbound or leng >= maxbound:

				continue

			else:

				cov+=1
				
				if cov > 100: #same as REFER
					
					break
					
				else:
					
					fasta+='>' + head + '\n' + seq + '\n'

		if cov >= coverage:

			with open(os.path.abspath(out+'/' + processor + '.unaligned.fa'),'w') as fastaout:

				fastaout.write(fasta.rstrip())
	else:

		cov=len(sequences)

	return cov


def Get_Alignment_Positions(bamfilein,softclipped):

  
	coords=[]
	seq=[]

	bamfile=pysam.AlignmentFile(bamfilein,'rb')

	for read in bamfile.fetch():

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			coords = read.get_reference_positions(full_length=True)
			seq=read.seq
			cigar=read.cigartuples

			cigardict=defaultdict(int)

			for operation,length in cigar:

				cigardict[operation] += length

			if 4 in cigardict.keys():

				if (cigardict[4]/len(seq))*100 >= softclipped:

					coords=[]
					seq=[]

	bamfile.close()

	return coords,seq


def modifier(coordinates):

		
	coordinates=[el+1 if el is not None else el for el in coordinates]
	start=next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
			
		if ele is None:

			coordinates[ind] = start
			
		else:

			start = ele

	return coordinates


def ParseBCF(bcfile):


	infos=[]
	it=VCF(bcfile)
	header=it.raw_header

	for variant in it:

		if variant.ALT == []:

			continue

		else:

			infos.append((variant.CHROM, variant.start+1, variant.INFO.get('TREND'), variant.REF, variant.ALT, variant.genotypes[0], variant.INFO.get('RAED'), variant.INFO.get('AED'), int(variant.format('DP1')), int(variant.format('DP2')), variant.INFO.get('TRLEN')))

	it.close()

	return infos, header


def VCF_HeaderModifier(rawheader, samples, output):


	headlist=rawheader.split('\n')[:-1]
	newheader=''

	for el in headlist:

		if el.startswith('##filedate') or el.startswith('##bcftools') or el.startswith('##source') or el.startswith('##INFO=<ID=H') or el.startswith('##INFO=<ID=MAPQ'):

			continue

		elif el.startswith('##FORMAT=<ID=DP2'):

			newheader+=el + '\n' + '##FORMAT=<ID=GS,Number=1,Type=Float,Description="Genotype Score (>=0 and <=2)">' + '\n'

		elif el.startswith('##SAMPLE'):

			newheader+=el + '\n'

			for name in samples:

				newheader+= '##SAMPLE=<ID=' + name + '>' + '\n'

		elif el.startswith('#CHROM'):

			for name in samples:

				el+='\t' + name

			newheader+=el + '\n'

		elif el.startswith('##INFO=<ID=AED'):

			newheader += el + '\n' + '##INFO=<ID=MISSR,Number=1,Type=Float,Description="Missing genotypes (ratio)">' +'\n' + '##INFO=<ID=MENDEL,Number=1,Type=Integer,Description="Mendelian consistency [.=Unknown;0=Consistent;1=Inconsistent]">' + '\n'

		else:

			newheader+=el + '\n'

		with open(os.path.abspath(output + '/TRiCoLOR.vcf'), 'w') as vcfout:

			vcfout.write(newheader)


def GetGTandGS(test,ref,alts):


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



def Runner(SHCpath,Cpath,gendir,processor,name,PROC_ENTRIES,sli,bam1,bam2,coverage,match, mismatch, gapopen, gapextend, output, readstype,softclipping):


	Entries = PROC_ENTRIES[processor] = list()

	for s in sli:

		chromosome,start,end,ref,alt,genchild,raed,aed,dp1,dp2,length=s

		if readstype == 'ONT':

			mmivar='map-ont'
			chromind= os.path.abspath(gendir + '/' + chromosome + '.ont.mmi')

		else:

			mmivar='map-pb'
			chromind= os.path.abspath(gendir + '/' + chromosome + '.pb.mmi')

		out1=os.path.abspath(output + '/' + name +'/haplotype1')
		out2=os.path.abspath(output + '/' + name +'/haplotype2')

		cov1=Bamfile_Analyzer(bam1,chromosome,start,end, coverage, out1, processor)
		cov2=Bamfile_Analyzer(bam2,chromosome,start,end, coverage, out2, processor)

		file1=os.path.abspath(out1 + '/' + processor + '.unaligned.fa')
		file2=os.path.abspath(out2 + '/' + processor + '.unaligned.fa')

		if os.path.exists(file1):

			subprocess.call(['bash', SHCpath, out1, Cpath, processor, os.path.basename(file1), mmivar, chromind, str(match), str(mismatch), str(gapopen), str(gapextend)],stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'))
			c_bam1=os.path.abspath(out1 + '/' + processor + '.cs.srt.bam')
			coords,seq=Get_Alignment_Positions(c_bam1,softclipping)
			os.remove(c_bam1)
			os.remove(c_bam1+'.bai')

			if seq == []:

				seq1 = ''

			else:

				coords1=modifier(coords)
				si=bisect_left(coords1, start)
				ei=bisect_right(coords1,end)
				seq1=seq[si:ei]

		else:

			seq1=''

		if os.path.exists(file2):

			subprocess.call(['bash', SHCpath, out2, Cpath, processor, os.path.basename(file2), mmivar, chromind, str(match), str(mismatch), str(gapopen), str(gapextend)],stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'))
			c_bam2=os.path.abspath(out2 + '/' + processor + '.cs.srt.bam')
			coords,seq=Get_Alignment_Positions(c_bam2,softclipping)
			os.remove(c_bam2)
			os.remove(c_bam2+'.bai')

			if seq == []:

				seq2 = ''

			else:

				coords2=modifier(coords)
				si=bisect_left(coords2, start)
				ei=bisect_right(coords2, end)
				seq2=seq[si:ei]

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

		Entries.append((chromosome,start,end,ref,alt,raed,aed,genchild,dp1,dp2,genotype,quality,seq1,seq2,cov1,cov2,length))

	PROC_ENTRIES[processor]=Entries



def GenotypeCombos():

	alleles=['0', '1', '2', '.']
	parent1=parent2=[x+'|'+y for x,y in itertools.product(alleles, alleles)]
	combos=set([(x,y) for x,y in itertools.product(parent1,parent2)])

	return combos


def GetPossibleGenotypes(genchild,combos):

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

	childdict=dict()
	genchildtype=['1|0', '0|1', '1|1', '1|2', '1|.', '.|1', '0|.', '.|0', '.|.']

	for key in genchildtype:

		childdict[key] = GetPossibleGenotypes(key, combos)

	return childdict


def CheckMendelian(genchild,genparent1,genparent2,childdict):


	if genchild not in childdict.keys(): #should not happen

		return '.'

	elif (genparent1,genparent2) in childdict[genchild]:

		return '0'

	else:

		return '1'
