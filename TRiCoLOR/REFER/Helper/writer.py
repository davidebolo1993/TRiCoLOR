#!/usr/bin/python env


#python 3 standard library

import datetime
import os
from bisect import bisect_left, bisect_right
from operator import itemgetter


#additional libraries

import pysam
import editdistance
import pandas as pd

## FUNCTIONS


def BED_repswriter(chromosome,coordreps, out):


	seq=[el[0] for el in coordreps]
	start=[el[1] for el in coordreps]
	end=[el[2] for el in coordreps]
	rep=[el[3] for el in coordreps]
	chrom=[chromosome]*len(start)

	Table=pd.DataFrame({'Chromosome':chrom, 'Start':start,'End':end, 'Repeated Motif':seq,'Repetitions Number':rep},columns=['Chromosome', 'Start', 'End', 'Repeated Motif', 'Repetitions Number'])
	
	if os.path.exists(os.path.abspath(out + '/' + chromosome + '.repetitions.bed')):

		with open(os.path.abspath(out + '/' + chromosome + '.repetitions.bed'), 'a') as refout:

			Table.to_csv(refout ,sep='\t',index=False, header=False)

	else:

		with open(os.path.abspath(out + '/' + chromosome + '.repetitions.bed'), 'w') as refout:

			Table.to_csv(refout ,sep='\t',index=False)


def VCF_headerwriter(bamfile1, bamfile2, samplename, commandline, out, processor):

	bam1=pysam.AlignmentFile(bamfile1,'rb')
	header1=bam1.header
	chromosomes_info1=list(header1.items())[1][1]

	if bamfile2 is not None:

		bam2=pysam.AlignmentFile(bamfile2,'rb')
		header2=bam2.header
		chromosomes_info2=list(header2.items())[1][1]
		chromosomes_info = list({x['SN']:x for x in chromosomes_info1 + chromosomes_info2}.values())

	else:

		chromosomes_info = list({x['SN']:x for x in chromosomes_info1}.values())

	chromosomes=[]
	sizes=[]

	for infos in chromosomes_info:

		chromosomes.append(infos['SN'])
		sizes.append(infos['LN'])

	vcf_format='##fileformat=VCFv4.2'

	SVEND='##INFO=<ID=SVEND,Number=1,Type=Integer,Description="Repetition end">'
	RAED = '##INFO=<ID=RAED,Number=1,Type=Integer,Description="Edit distance between REF and most similar ALT allele">'
	AED = '##INFO=<ID=AED,Number=1,Type=Integer,Description="Edit distance between ALT alleles">'
	H1M='##INFO=<ID=H1M,Number=.,Type=String,Description="Haplotype1 Repeated Motif">'
	H1N='##INFO=<ID=H1N,Number=.,Type=Integer,Description="Haplotype1 Repetitions Number">'
	H2M='##INFO=<ID=H2M,Number=.,Type=String,Description="Haplotype2 Repeated Motif">'
	H2N='##INFO=<ID=H2N,Number=.,Type=Integer,Description="Haplotype2 Repetitions Number">'
	FORMAT1='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
	FORMAT2 = '##FORMAT=<ID=DP1,Number=1,Type=Integer,Description="Coverage depth for 1st haplotype">'
	FORMAT3 = '##FORMAT=<ID=DP2,Number=1,Type=Integer,Description="Coverage depth for 2nd haplotype">'

	classic_header='#CHROM' + '\t' + 'POS' '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' + '\t' + 'INFO' + '\t' + 'FORMAT' + '\t' + samplename.upper()

	with open(os.path.abspath(out + '/' + processor + '.TRiCoLOR.vcf'), 'w') as vcfout:

		vcfout.write(vcf_format + '\n' + '##filedate=' + str(datetime.date.today()) +  '\n' + '##source=' + commandline + '\n')

		for a,b in zip(chromosomes,sizes):

			vcfout.write('##contig=<ID='+str(a)+',length='+str(b)+'>'+'\n')

		vcfout.write(SVEND + '\n' + RAED + '\n' + AED + '\n' + H1M + '\n' + H1N + '\n' + H2M + '\n' + H2N + '\n')
		vcfout.write(FORMAT1 + '\n' + FORMAT2 + '\n' + FORMAT3 + '\n')
		vcfout.write('##SAMPLE=<ID=' + samplename +'>' + '\n' + classic_header + '\n')


def VCF_variantwriter(CHROM, POS, REF, ALT, INFO, FORMAT, out, processor):


	ID='.'
	FILTER='.'
	QUAL='.'
	GEN='GT:DP1:DP2'
	
	INFO_SVEND=str(INFO['SVEND'])
	INFO_RAED = str(INFO['RAED'])
	INFO_AED = str(INFO['AED'])
	INFO_H1M=INFO['H1M']

	if type(INFO_H1M) == list:

		INFO_H1M = ','.join(str(x) for x in INFO_H1M) 

	INFO_H1N=INFO['H1N']

	if type(INFO_H1N) == list:

		INFO_H1N = ','.join(str(x) for x in INFO_H1N) 

	INFO_H2M=INFO['H2M']

	if type(INFO_H2M) == list:

		INFO_H2M = ','.join(str(x) for x in INFO_H2M) 

	INFO_H2N=INFO['H2N']

	if type(INFO_H2N) == list:

		INFO_H2N = ','.join(str(x) for x in INFO_H2N) 

	FORMAT_GT=FORMAT['GT']
	FORMAT_DP1 = str(FORMAT['DP1'])
	FORMAT_DP2 = str(FORMAT['DP2'])

	with open(os.path.abspath(out + '/' + processor + '.TRiCoLOR.vcf'), 'a') as vcfout:

		vcfout.write(CHROM + '\t' + str(POS) + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER + '\t' + 'SVEND='+INFO_SVEND + ';'+ 'RAED='+ INFO_RAED + ';' + 'AED=' + INFO_AED  + ';' + 'H1M='+INFO_H1M + ';' + 'H1N='+INFO_H1N + ';' + 'H2M='+INFO_H2M + ';' + 'H2N='+INFO_H2N + '\t' + GEN + '\t' + FORMAT_GT + ':' + FORMAT_DP1 + ':' + FORMAT_DP2 + '\n')


def modifier2(seq,coords,POS,SVEND):

	NewSeq=''
	coords_purified=[]

	for i in range(len(coords)-1):

		if coords[i+1]-coords[i] > 1:

			coords_purified.append(coords[i])
			coords_purified.extend(list(range(coords[i]+1,coords[i+1])))
			NewSeq+=seq[i]
			NewSeq+="-"*(coords[i+1]-coords[i]-1)

		else:

			coords_purified.append(coords[i])
			NewSeq+=seq[i]

	coords_purified.append(coords[-1])
	NewSeq+=seq[-1]

	if not POS >= coords_purified[0]:

		number=POS
		how_many=coords_purified[0]-POS
		coords_purified = [number]*how_many + coords_purified
		NewSeq= '-'* how_many + NewSeq

	if not SVEND <= coords_purified[-1]:

		number=SVEND
		how_many=SVEND - coords_purified[-1]
		coords_purified = coords_purified + [number]*how_many 
		NewSeq= NewSeq + '-'* how_many

	return NewSeq,coords_purified


def GetIndex(start, end, coordinates):


	si=bisect_left(coordinates, start)
	ei=bisect_right(coordinates, end)

	return si,ei


def recursive_merge(sorted_int, list_, i):


	new_=(min(list_, key=itemgetter(1)), max(list_,key=itemgetter(2)))
	new_range=(new_[0][1], new_[-1][2])

	if i < len(sorted_int) -1:

		if sorted_int[i+1][1] <= new_range[1]:

			list_.append(sorted_int[i+1])
			recursive_merge(sorted_int, list_, i+1)


def Merger(sorted_int, refreps, h1reps, h2reps):


	sorted_ranges=[]

	ref_dict_number=dict()
	ref_dict_motif=dict()

	hap1_dict_number=dict()
	hap1_dict_motif=dict()

	hap2_dict_number=dict()
	hap2_dict_motif=dict()

	i=0

	while i < len(sorted_int):

		reps=sorted_int[i]
		to_int=sorted_int[i+1:]
		list_=[]
		l=0

		for elem in to_int:

			if elem[1] <= reps[2]:

				l+=1

				list_.append(elem)


		if len(list_) != 0:

			list_.append(reps)
			recursive_merge(sorted_int, list_, i+len(list_)-1) #? ANY BETTER IDEA THAN RECURSION? NOT PRIORITY.
			new_=(min(list_, key=itemgetter(1)), max(list_,key=itemgetter(2)))
			new_range=(new_[0][1], new_[-1][2])

			for el_ in list_:

				if el_ in refreps:

					if new_range not in ref_dict_motif:

						ref_dict_motif[new_range]= [el_[0]]
						ref_dict_number[new_range]= [el_[3]]

					else:

						ref_dict_motif[new_range].append(el_[0])
						ref_dict_number[new_range].append(el_[3])

				if el_ in h1reps:

					if new_range not in hap1_dict_motif:

						hap1_dict_motif[new_range]= [el_[0]]
						hap1_dict_number[new_range]= [el_[3]]

					else:

						hap1_dict_motif[new_range].append(el_[0])
						hap1_dict_number[new_range].append(el_[3])

				if el_ in h2reps:

					if new_range not in hap2_dict_motif:

						hap2_dict_motif[new_range]= [el_[0]]
						hap2_dict_number[new_range]= [el_[3]]

					else:

						hap2_dict_motif[new_range].append(el_[0])
						hap2_dict_number[new_range].append(el_[3])

			i+=len(list_)
			sorted_ranges.append(new_range)

		else:

			new_range=((reps[1], reps[2]))

			if reps in refreps:

				if new_range not in ref_dict_motif:

					ref_dict_motif[new_range]= [reps[0]]
					ref_dict_number[new_range]= [reps[3]]

				else:

					ref_dict_motif[new_range].append(reps[0])
					ref_dict_number[new_range].append(reps[3])

			if reps in h1reps:

				if new_range not in hap1_dict_motif:

					hap1_dict_motif[new_range]= [reps[0]]
					hap1_dict_number[new_range]= [reps[3]]

				else:

					hap1_dict_motif[new_range].append(reps[0])
					hap1_dict_number[new_range].append(reps[3])


			if reps in h2reps:

				if new_range not in hap2_dict_motif:

					hap2_dict_motif[new_range]= [reps[0]]
					hap2_dict_number[new_range]= [reps[3]]

				else:

					hap2_dict_motif[new_range].append(reps[0])
					hap2_dict_number[new_range].append(reps[3])

			if not new_range in sorted_ranges:

				sorted_ranges.append(new_range)

			i += 1

	return sorted_ranges,ref_dict_number,ref_dict_motif,hap1_dict_number,hap1_dict_motif,hap2_dict_number,hap2_dict_motif


def VCF_writer(chromosome, repref, reference_sequence, repsh1, seqh1, coordsh1, covh1, repsh2, seqh2, coordsh2, covh2, out, processor):

	intersection=list(set(repref+ repsh1 + repsh2))

	if len(intersection) != 0:

		sorted_intersection=sorted(intersection, key=itemgetter(1,2))
		sorted_ranges,ref_dict_number,ref_dict_motif,hap1_dict_number,hap1_dict_motif,hap2_dict_number,hap2_dict_motif=Merger(sorted_intersection, repref, repsh1, repsh2)

		for reps in sorted_ranges:

			CHROM=chromosome
			POS=reps[0]
			SVEND=reps[1]
			REF=reference_sequence[(POS-1):SVEND]

			if reps in hap1_dict_number.keys():

				seqh1_,coordsh1_=modifier2(seqh1,coordsh1,POS,SVEND)
				IS1,IE1=GetIndex(POS,SVEND,coordsh1_)
				ALT1=seqh1_[IS1:IE1].replace('-','')
				H1N=hap1_dict_number[reps]
				H1M=hap1_dict_motif[reps]
				DP1=covh1

			else:

				if seqh1 == []:

					ALT1='.'
					H1N='.'
					H1M='.'
					
					if covh1 == []:

						DP1='.'

					else:

						DP1=covh1

				else:

					seqh1_,coordsh1_=modifier2(seqh1,coordsh1,POS,SVEND)
					IS1,IE1=GetIndex(POS,SVEND,coordsh1_)
					ALT1=seqh1_[IS1:IE1].replace('-','')
					H1N='.'
					H1M='.'
					DP1=covh1

			if reps in hap2_dict_number.keys():

				seqh2_,coordsh2_=modifier2(seqh2,coordsh2,POS,SVEND)
				IS2,IE2=GetIndex(POS,SVEND,coordsh2_)
				ALT2=seqh2_[IS2:IE2].replace('-','')
				H2N=hap2_dict_number[reps]
				H2M=hap2_dict_motif[reps]
				DP2=covh2

			else:

				if seqh2 == []:

					ALT2='.'
					H2N='.'
					H2M='.'
					
					if covh2 == []:

						DP2='.'

					else:

						DP2=covh2

				else:

					seqh2_,coordsh2_=modifier2(seqh2,coordsh2,POS,SVEND)
					IS2,IE2=GetIndex(POS,SVEND,coordsh2_)
					ALT2=seqh2_[IS2:IE2].replace('-','')
					H2N='.'
					H2M='.'
					DP2=covh2

			if seqh1 == [] and seqh2 == []:

				GEN1='.'
				GEN2='.'
				ALT = '.'

			elif seqh1 != [] and seqh2 == []:

				GEN2='.'

				if ALT1 == REF:

					GEN1 = '0'
					ALT = '.'

				else:

					GEN1 = '1'
					ALT = ALT1

			elif seqh1 == [] and seqh2 != []:

				GEN1='.'

				if ALT2 == REF:

					GEN2 = '0'
					ALT = '.'

				else:

					GEN2 = '1'
					ALT = ALT2

			else:

				if ALT1 == REF and ALT2 == REF:

					#GEN1 = '0'
					#GEN2 = '0'
					#ALT= '.'
					continue

				elif ALT1 == REF and ALT2 != REF:

					GEN1= '0'
					GEN2 = '1'
					ALT = ALT2

				elif ALT1 != REF and ALT2 == REF:

					GEN1= '1'
					GEN2 = '0'
					ALT = ALT1

				else:

					if ALT1 == ALT2:

						GEN1 = '1'
						GEN2 = '1'
						ALT = ALT1

					else:

						GEN1 = '1'
						GEN2 = '2'
						ALT = ALT1 + ',' + ALT2

			GEN = GEN1 + '|' + GEN2
			INFO=dict()
			INFO['SVEND'] = SVEND

			if ALT1 == '.' and ALT2 == '.':

				RAED = '.'
				AED = '.'

			elif ALT1 != '.' and ALT2 == '.':

				RAED = editdistance.eval(REF, ALT1)
				AED= '.'

			elif ALT1 == '.' and ALT2 != '.':

				RAED = editdistance.eval(REF, ALT2)
				AED='.'

			else:

				RAED1=editdistance.eval(REF, ALT1)
				RAED2=editdistance.eval(REF, ALT2)

				if RAED1 < RAED2:

					RAED = RAED1

				else:

					RAED = RAED2

				AED = editdistance.eval(ALT1,ALT2)

			INFO['RAED'] = RAED
			INFO['AED'] = AED
			INFO['H1M'] = H1M
			INFO['H1N'] = H1N
			INFO['H2M'] = H2M 
			INFO['H2N'] = H2N


			FORMAT = dict()

			FORMAT['GT'] = GEN
			FORMAT['DP1'] = DP1
			FORMAT['DP2'] = DP2

			VCF_variantwriter(CHROM, POS, REF, ALT, INFO, FORMAT, out ,processor)
