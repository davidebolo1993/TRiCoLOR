#python 3 standard library

import sys
import os
import math
import re
import subprocess
import logging
import multiprocessing
from shutil import which,rmtree
from bisect import bisect_left,bisect_right

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

			print('The output folder is not empty. Specify another output folder or clean the previsouly chosen')
			sys.exit(1)

	command_dict= vars(args)
	
	notkey=['func']
	command_string= ' '.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey)
	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR_SAGE.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	
	print('Initialized .log file ' + os.path.abspath(args.output + '/TRiCoLOR_SAGE.log'))

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

	SHCpath=os.path.abspath(os.path.dirname(os.path.dirname(__file__)) + '/REFER/consensus.sh')
	Cpath=os.path.abspath(os.path.dirname(os.path.dirname(__file__)) + '/REFER/alfred/bin/alfred') #? FASTER AND STAND-ALONE CONSENSUS. PRIMARY TASK
	gendir=os.path.abspath(os.path.dirname(args.genome))

	logging.info('Related individuals to genotype: ' + str(len(snames)))

	if args.threads > multiprocessing.cpu_count():

		cores=multiprocessing.cpu_count()
		logging.warning(str(args.threads) + ' cores specified but only ' + str(cores) + ' available. Using all cores available')

	else:

		cores=args.threads

	logging.info('Coverage treshold: ' + str(args.coverage))
	logging.info('Long reads type: ' + str(args.readstype))
	logging.info('Cores: ' + str(cores))

	logging.info('Parsing input BCF ...')

	infos,header=GetInfo(os.path.abspath(args.bcffile))

	logging.info('Done')

	VCF_HeaderModifier(header,snames,os.path.abspath(args.output))

	entries=len(infos)
	chunk_size=entries/cores
	slices=Chunks(infos,math.ceil(chunk_size))

	Namesdict=dict()

	logging.info('Genotype samples ...')

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
			p=multiprocessing.Process(target=Runner, args=(SHCpath,Cpath,gendir,processor,names,PROC_ENTRIES,sli,couples[0],couples[1],args.coverage,args.readstype,os.path.abspath(args.output)))
			p.start()
			processes.append(p)
		
		for p in processes:
		
			p.join()

		Namesdict[names] = []


		for key in sorted(PROC_ENTRIES.keys(), key=natural_keys):

			Namesdict[names].extend(PROC_ENTRIES[key])
			print(key)
			print(PROC_ENTRIES[key])

		rmtree(os.path.abspath(args.output) + '/' + names + '/haplotype1')
		rmtree(os.path.abspath(args.output) + '/' + names + '/haplotype2')
		os.rmdir(os.path.abspath(args.output) + '/' + names)

		logging.info('Done')

	logging.info('Writing to BCF ...')

	FILTER='.'
	QUAL='.'
	ID='.'
	FORMAT='GT:GS'
	QUALCHILD=str(1.0)

	with open(os.path.abspath(args.output + '/TRiCoLOR_SAGE.vcf'), 'a') as vcfout:

		for i in range(len(infos)):

			CHROM,POS,END,REF,ALT,GENCHILD=Namesdict[snames[0]][i][0], str(Namesdict[snames[0]][i][1]), str(Namesdict[snames[0]][i][2]), Namesdict[snames[0]][i][3], ','.join(x for x in Namesdict[snames[0]][i][4]), '|'.join(str(x) for x in Namesdict[snames[0]][i][5][:-1]).replace('-1', '.')

			toadd=[]

			for names in snames:

				toadd.append(Namesdict[names][i][6] + ':' + str(Namesdict[names][i][7]))

			vcfout.write(CHROM + '\t' + POS + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER + '\t' + 'END='+END + '\t' + FORMAT + '\t' + GENCHILD + ':' +  QUALCHILD+ '\t' + '\t'.join(x for x in toadd) + '\n')

	subprocess.call(['bcftools', 'sort', '-o', os.path.abspath(args.output + '/TRiCoLOR_SAGE.srt.bcf'), '-O', 'b', os.path.abspath(args.output + '/TRiCoLOR_SAGE.vcf')],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
	subprocess.call(['bcftools', 'index', os.path.abspath(args.output + '/TRiCoLOR_SAGE.srt.bcf')],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
	os.remove(os.path.abspath(args.output + '/TRiCoLOR_SAGE.vcf'))

	logging.info('Done')
	print('Done')

	sys.exit(0)


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


def GetInfo(bcfile):


	infos=[]
	it=VCF(bcfile)
	header=it.raw_header

	for variant in it:

		if variant.ALT == []:

			continue

		else:

			infos.append((variant.CHROM, variant.start+1, variant.end, variant.REF, variant.ALT, variant.genotypes[0]))

	it.close()

	return infos, header


def VCF_HeaderModifier(rawheader, samples, output):


	headlist=rawheader.split('\n')[:-1]
	newheader=''
	toadd='##FORMAT=<ID=GS,Number=1,Type=Float,Description="Genotype Similarity (ranges from 0.0 to 1.0)">' + '\n'

	for el in headlist:

		if el.startswith('##FORMAT'):

			newheader+=el + '\n'

			if (toadd not in x for x in headlist):

				newheader+= toadd

		elif el.startswith('##filedate') or el.startswith('##filedate') or el.startswith('##bcftools') or el.startswith('##source') or el.startswith('##INFO=<ID=H'):

			continue
			
		elif el.startswith('##SAMPLE'):

			newheader+=el + '\n'

			for name in samples:

				newheader+= '##SAMPLE=<ID=' + name + '>' + '\n'

		elif el.startswith('#CHROM'):

			for name in samples:

				el+='\t' + name

			newheader+=el + '\n'

		else:

			newheader+=el + '\n'

	with open(os.path.abspath(output + '/TRiCoLOR_SAGE.vcf'), 'w') as vcfout:

		vcfout.write(newheader)


def sub_none(list_of_coord):


	return [-999999 if v is None else v for v in list_of_coord]


def check_coverage(pysam_AlignmentFile, chromosome, start, end, coverage):


	counter = 0

	for read in pysam_AlignmentFile.fetch(chromosome, start, end):
		
		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			if read.reference_start <= start and read.reference_end >= end: 

				counter +=1

	if counter >= coverage:

		return True

	else:

		return False


def Bamfile_Analyzer(bamfilein,chromosome,start,end, coverage, out, processor):


	start=start-350
	end=end+350

	bamfile=pysam.AlignmentFile(bamfilein,'rb')	

	if check_coverage(bamfile, chromosome, start, end, coverage):

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

					with open(os.path.abspath(out+'/' + processor + '.unaligned.fa'),'a') as fastaout:

						fastaout.write('>' + header + '\n' + sequence[s_i:e_i+1] + '\n')
	
	bamfile.close()


def Get_Alignment_Positions(bamfilein):

  
    coords=[]
    seq=[]

    bamfile=pysam.AlignmentFile(bamfilein,'rb')

    for read in bamfile.fetch():

        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

            coords = read.get_reference_positions(full_length=True)
            seq=read.seq

    bamfile.close()

    return coords,seq


def GetGenotypeAndScore(test,ref,alt):

	refedit=editdistance.eval(ref,test)
	invertedRscore = (1.0-(refedit/len(test)))/2

	alts=[]

	for alteration in alt:

		altedit=editdistance.eval(alteration,test)
		invertedAscore = (1.0-(altedit/len(test)))/2
		alts.append(invertedAscore)

	ind = alts.index(max(alts))

	if alts[ind] >= invertedRscore:

		return str(ind +1),round(float(alts[ind]),2)

	elif alts[ind] < invertedRscore:

		return '0', round(float(invertedRscore),2)

	#else:

		#return '.', round(float(0))


def modifier(coordinates):

    
    coordinates=[el+1 if el is not None else el for el in coordinates]
    start=next(ele for ele in coordinates if ele is not None)

    for ind, ele in enumerate(coordinates):
        
        if ele is None:

            coordinates[ind] = start
        
        else:

            start = ele

    return coordinates


def Runner(SHCpath,Cpath,gendir,processor,name,PROC_ENTRIES,sli,bam1,bam2,coverage,readtype,output):


	Entries = PROC_ENTRIES[processor] = list()

	for s in sli:

		chromosome,start,end,ref,alt,gen=s

		if readtype == 'ONT':

			mmivar='map-ont'
			consvar='ont'
			chromind= os.path.abspath(gendir + '/' + chromosome + '.ont.mmi')

		else:

			mmivar='map-pb'
			consvar='pacbio'
			chromind= os.path.abspath(gendir + '/' + chromosome + '.pb.mmi')

		out1=os.path.abspath(output + '/' + name +'/haplotype1')
		out2=os.path.abspath(output + '/' + name +'/haplotype2')

		Bamfile_Analyzer(bam1,chromosome,start,end, coverage, out1, processor)
		Bamfile_Analyzer(bam2,chromosome,start,end, coverage, out2, processor)

		file1=os.path.abspath(out1 + '/' + processor + '.unaligned.fa')
		file2=os.path.abspath(out2 + '/' + processor + '.unaligned.fa')

		if os.path.exists(file1):

			subprocess.call(['bash', SHCpath, out1, Cpath, consvar, processor, os.path.basename(file1), mmivar, chromind],stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'))
			c_bam1=os.path.abspath(out1 + '/' + processor + '.cs.srt.bam')
			coords,seq=Get_Alignment_Positions(c_bam1)

			if seq == []:

				seq1 = ''

			else:

				coords1=modifier(coords)
				si=bisect_left(coords1, start)
				ei=bisect_right(coords1, end)
				seq1=seq[si:ei]

		else:

			seq1=''

		if os.path.exists(file2):

			subprocess.call(['bash', SHCpath, out2, Cpath, consvar, processor, os.path.basename(file2), mmivar, chromind],stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'))
			c_bam2=os.path.abspath(out2 + '/' + processor + '.cs.srt.bam')
			coords,seq=Get_Alignment_Positions(c_bam2)

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

			gen1,qual1=GetGenotypeAndScore(seq1,ref,alt)

			genotype=gen1 + '|.'

			if qual1 == '.':

				quality=float(0)

			else:

				quality=qual1

		elif seq1 == '' and seq2 != '':

			gen2,qual2=GetGenotypeAndScore(seq2,ref,alt)

			genotype='.|' + gen2

			if qual2 == '.':

				quality=float(0)

			else:

				quality=qual2

		else:

			gen1,qual1=GetGenotypeAndScore(seq1,ref,alt)
			gen2,qual2=GetGenotypeAndScore(seq2,ref,alt)

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

		Entries.append((chromosome,start,end,ref,alt,gen, genotype, quality))


	PROC_ENTRIES[processor]=Entries


if __name__ == '__main__':

	main()
