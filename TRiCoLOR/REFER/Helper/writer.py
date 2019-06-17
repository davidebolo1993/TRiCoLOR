#!/usr/bin/python env


#python 3 standard library

import datetime
import os
from bisect import bisect_left, bisect_right
from collections import defaultdict
from operator import itemgetter
import subprocess


#additional libraries

import pysam


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

	END='##INFO=<ID=END,Number=1,Type=Integer,Description="Repetition end">'
	H1M='##INFO=<ID=H1M,Number=.,Type=String,Description="Haplotype1 Repeated Motif">'
	H1N='##INFO=<ID=H1N,Number=.,Type=String,Description="Haplotype1 Repetitions Number">'
	H2M='##INFO=<ID=H2M,Number=.,Type=String,Description="Haplotype2 Repeated Motif">'
	H2N='##INFO=<ID=H2N,Number=.,Type=String,Description="Haplotype2 Repetitions Number">'
	FORMAT='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
	classic_header='#CHROM' + '\t' + 'POS' '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' + '\t' + 'INFO' + '\t' + 'FORMAT' + '\t' + samplename.upper()

	with open(os.path.abspath(out + '/' + processor + '.TRiCoLOR.vcf'), 'w') as vcfout:

		vcfout.write(vcf_format + '\n' + '##filedate=' + str(datetime.date.today()) +  '\n' + '##source=' + commandline + '\n')

		for a,b in zip(chromosomes,sizes):

			vcfout.write('##contig=<ID='+str(a)+',length='+str(b)+'>'+'\n')

		vcfout.write(END + '\n' + H1M + '\n' + H1N + '\n' + H2M + '\n' + H2N + '\n')
		vcfout.write(FORMAT + '\n')
		vcfout.write('##SAMPLE=<ID=' + samplename +'>' + '\n' + classic_header + '\n')


def VCF_variantwriter(chrom, pos, ref, alt, info, form, out, processor):


	ID='.'
	FILTER='.'
	QUAL='.'
	GEN='GT'
	CHROM=chrom
	POS=str(pos)
	REF=ref
	ALT=alt
	INFO_END=str(info['END'])

	INFO_H1M=info['H1M']

	if type(INFO_H1M) == list:

		INFO_H1M = ','.join(str(x) for x in INFO_H1M) 

	INFO_H1N=info['H1N']

	if type(INFO_H1N) == list:

		INFO_H1N = ','.join(str(x) for x in INFO_H1N) 

	INFO_H2M=info['H2M']

	if type(INFO_H2M) == list:

		INFO_H2M = ','.join(str(x) for x in INFO_H2M) 

	INFO_H2N=info['H2N']

	if type(INFO_H2N) == list:

		INFO_H2N = ','.join(str(x) for x in INFO_H2N) 

	FORMAT=form


	with open(os.path.abspath(out + '/' + processor + '.TRiCoLOR.vcf'), 'a') as vcfout:

		vcfout.write(CHROM + '\t' + POS + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER + '\t' + 'END='+INFO_END + ';'+ 'H1M='+INFO_H1M + ';' + 'H1N='+INFO_H1N + ';' + 'H2M='+INFO_H2M + ';' + 'H2N='+INFO_H2N + '\t' + GEN + '\t' + FORMAT + '\n')


def modifier(coordinates):

	
	coordinates=[el+1 if el is not None else el for el in coordinates]
	start=next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:

			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates


def Modifier(list_of_coord,seq,reps):


	coords_without_insertions=modifier(list_of_coord)	
	NewSeq=''
	coords_purified=[]

	for i in range(len(coords_without_insertions)-1):

		if coords_without_insertions[i+1]-coords_without_insertions[i] > 1:

			coords_purified.append(coords_without_insertions[i])
			coords_purified.extend(list(range(coords_without_insertions[i]+1,coords_without_insertions[i+1])))
			NewSeq+=seq[i]
			NewSeq+="-"*(coords_without_insertions[i+1]-coords_without_insertions[i]-1)

		else:

			coords_purified.append(coords_without_insertions[i])
			NewSeq+=seq[i]

	coords_purified.append(coords_without_insertions[-1])
	NewSeq+=seq[-1]

	if not reps[0] >= coords_purified[0]:

		number=reps[0]
		how_many=coords_purified[0]-reps[0]
		coords_purified = [number]*how_many + coords_purified
		NewSeq= '-'* how_many + NewSeq


	if not reps[1] <= coords_purified[-1]:

		number=reps[1]
		how_many=reps[1] - coords_purified[-1]
		coords_purified = coords_purified + [number]*how_many 
		NewSeq= NewSeq + '-'* how_many

	return coords_purified,NewSeq


def Get_Seq_Pos(bamfilein,chromosome, start,end):
	

	seq=[]
	coords=[]
	
	if os.stat(os.path.abspath(bamfilein)).st_size == 0:

		return seq,coords

	else:

		if start==end:

			end=start+1

		bamfile=pysam.AlignmentFile(bamfilein,'rb')

		for read in bamfile.fetch(chromosome, start-1, end-1):

			if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

				coords = read.get_reference_positions(full_length=True)
				seq=read.seq

	return seq,coords


def GetIndex(start, end, coordinates):


	si=bisect_left(coordinates, start)
	ei=bisect_right(coordinates, end)-1

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


def VCF_writer(chromosome, repref, reference_sequence, repsh1, bamfile1, repsh2, bamfile2, out, processor):


	intersection=list(set(repref+ repsh1 + repsh2))

	if len(intersection) == 0:

		return

	else:

		sorted_intersection=sorted(intersection, key=itemgetter(1,2))
		sorted_ranges,ref_dict_number,ref_dict_motif,hap1_dict_number,hap1_dict_motif,hap2_dict_number,hap2_dict_motif=Merger(sorted_intersection, repref, repsh1, repsh2)

		for reps in sorted_ranges:

			if reps in ref_dict_number.keys() and reps in hap1_dict_number.keys() and reps in hap2_dict_number.keys():

				if ref_dict_number[reps] == hap1_dict_number[reps] and ref_dict_number[reps] == hap2_dict_number[reps]: 

					if  ref_dict_motif[reps] == hap1_dict_motif[reps] and ref_dict_motif[reps] == hap2_dict_motif[reps]: 

						continue 

					elif ref_dict_motif[reps] == hap1_dict_motif[reps] and ref_dict_motif[reps] != hap2_dict_motif[reps]: 

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]
						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps)
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

						if ref == alt2:

							continue

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]			
							form='0|1' 

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

					elif ref_dict_motif[reps] != hap1_dict_motif[reps] and ref_dict_motif[reps] == hap2_dict_motif[reps]: 

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]
						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1])
						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps)
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

						if ref == alt1:

							continue

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]			
							form='1|0' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

					elif ref_dict_motif[reps] != hap1_dict_motif[reps] and ref_dict_motif[reps] != hap2_dict_motif[reps]: 

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]
						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 
						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps)
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps)
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

						if ref == alt1 and ref == alt2:

							continue

						elif ref == alt1 and ref != alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = ref_dict_motif[reps] 
							info['H1N'] = ref_dict_number[reps] 
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]			
							form='0|1' 

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

						elif ref != alt1 and ref == alt2:

							info=dict()			

							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = ref_dict_motif[reps] 
							info['H2N'] = ref_dict_number[reps] 
							form='1|0' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

						elif ref != alt1 and ref != alt2:

							if alt1 != alt2 :

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]			
								form='1|2' 

								VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]			
								form='1|1' 

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

				elif ref_dict_number[reps] == hap1_dict_number[reps] and ref_dict_number[reps] != hap2_dict_number[reps]: 

					if ref_dict_motif[reps] == hap1_dict_motif[reps]:

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]
						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps)
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

						if ref == alt2:

							continue

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]			
							form='0|1'

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

					else:

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]
						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 
						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps)
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps)
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

						if ref == alt1 and ref == alt2:

							continue

						elif ref == alt1 and ref != alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = ref_dict_motif[reps] 
							info['H1N'] = ref_dict_number[reps] 
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]			
							form='0|1' 

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

						elif ref != alt1 and ref == alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = ref_dict_motif[reps]  
							info['H2N'] = ref_dict_number[reps] 
							form='1|0' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)
						
						elif ref != alt1 and ref != alt2:

							if alt1 != alt2 :

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|2' 

								VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|1' 

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)


				elif ref_dict_number[reps] != hap1_dict_number[reps] and ref_dict_number[reps] == hap2_dict_number[reps]: 

					if ref_dict_motif[reps] == hap2_dict_motif[reps]: 

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1])
						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')


						if ref == alt1: 

							continue

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='1|0' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

					else: 


						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 


						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)

						alt1=seq_h1[si_1:(ei_1+1)].replace('-','') 
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','') 


						if ref == alt1 and ref == alt2: 

							continue

						elif ref == alt1 and ref != alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = ref_dict_motif[reps] 
							info['H1N'] = ref_dict_number[reps] 
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='0|1' 

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

						elif ref != alt1 and ref == alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = ref_dict_motif[reps] 
							info['H2N'] = ref_dict_number[reps] 
			
							form='1|0' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

						elif ref != alt1 and ref != alt2:

							if alt1 != alt2 :

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|2' 

								VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|1' 

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)


				else: 

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1])
					seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])


					coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
					coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 

					si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
					si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)

					alt1=seq_h1[si_1:(ei_1+1)].replace('-','') 
					alt2=seq_h2[si_2:(ei_2+1)].replace('-','') 


					if ref == alt1 and ref == alt2: 

						continue

					elif ref == alt1 and ref != alt2:

						info=dict()
			
						info['END'] = reps[1]
						info['H1M'] = ref_dict_motif[reps] 
						info['H1N'] = ref_dict_number[reps] 
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
			
						form='0|1' 

						VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

					elif ref != alt1 and ref == alt2:

						info=dict()
			
						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = ref_dict_motif[reps] 
						info['H2N'] = ref_dict_number[reps] 
			
						form='1|0' 

						VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

					elif ref != alt1 and ref != alt2:

						if alt1 != alt2 :

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='1|2' 

							VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='1|1' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)


			elif reps in ref_dict_number.keys() and reps in hap1_dict_number.keys() and reps not in hap2_dict_number.keys(): 

				if ref_dict_number[reps] == hap1_dict_number[reps] and ref_dict_motif[reps] == hap1_dict_motif[reps]: 

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 

					if len(seq_h2) == 0: 

						info=dict()
			
						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = '.'
						info['H2N'] = '.'
			
						form='0|.'  

						VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

					else: 

						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)						
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','') 

						if alt2 == '':

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = '.'
							info['H2N'] = '.'
			
							form='0|.' 

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

						else:

							if ref==alt2:

								continue 

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = '.'
								info['H2N'] = '.'
			
								form='0|1' 

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

				
				else: 

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
					coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
					si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
					alt1=seq_h1[si_1:(ei_1+1)].replace('-','')
					
					seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 

					if ref==alt1:

						if len(seq_h2) == 0: 

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = ref_dict_motif[reps] 
							info['H1N'] = ref_dict_number[reps] 
							info['H2M'] = '.'
							info['H2N'] = '.'
			
							form='0|.' 

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

						else: 

							coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
							si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
							alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

							if alt2 =='': 

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = ref_dict_motif[reps] 
								info['H1N'] = ref_dict_number[reps] 
								info['H2M'] = '.'
								info['H2N'] = '.'
			
								form='0|.' 

								VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

							else:

								if ref==alt2:

									continue

								else:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = ref_dict_motif[reps] 
									info['H1N'] = ref_dict_number[reps] 
									info['H2M'] = '.'
									info['H2N'] = '.'
			
									form='0|1' 

									VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

					else:

						if len(seq_h2) == 0: 

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = '.'
							info['H2N'] = '.'
			
							form='1|.' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

						else: 

							coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
							si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
							alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

							if alt2=='': 

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = '.'
								info['H2N'] = '.'
			
								form='1|.' 


								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

							else:

								if ref==alt2:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = hap1_dict_motif[reps]
									info['H1N'] = hap1_dict_number[reps]
									info['H2M'] = '.'
									info['H2N'] = '.'
			
									form='1|0' 

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

								else:

									if alt1==alt2:

										info=dict()
			
										info['END'] = reps[1]
										info['H1M'] = hap1_dict_motif[reps]
										info['H1N'] = hap1_dict_number[reps]
										info['H2M'] = '.'
										info['H2N'] = '.'
			
										form='1|1' 


										VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

									else:

										info=dict()
			
										info['END'] = reps[1]
										info['H1M'] = hap1_dict_motif[reps]
										info['H1N'] = hap1_dict_number[reps]
										info['H2M'] = '.'
										info['H2N'] = '.'
			
										form='1|2' 


										VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)


			elif reps in ref_dict_number.keys() and reps not in hap1_dict_number.keys() and reps in hap2_dict_number.keys(): 

				if ref_dict_number[reps] == hap2_dict_number[reps] and ref_dict_motif[reps] == hap2_dict_motif[reps]: 

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1])

					if len(seq_h1) == 0: 

						info=dict()
			
						info['END'] = reps[1]
						info['H1M'] = '.'
						info['H1N'] = '.'
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
			
						form='.|0' 

						VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

					else: 

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')


						if alt1=='':

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='.|0' 

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

						else:

							if ref==alt1:

								continue

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|0' 

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)
				
				else: 

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])
					coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
					si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
					alt2=seq_h2[si_2:(ei_2+1)].replace('-','')
					
					seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 

					if ref==alt2:

						if len(seq_h1) == 0: 

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = ref_dict_motif[reps] 
							info['H2N'] = ref_dict_number[reps] 
			
							form='.|0' 

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

						else: 

							coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
							si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
							alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

							if alt1=='': 

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = ref_dict_motif[reps] 
								info['H2N'] = ref_dict_number[reps] 
			
								form='.|0' 

								VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

							else:

								if ref==alt1:

									continue

								else:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = ref_dict_motif[reps] 
									info['H2N'] = ref_dict_number[reps] 
			
									form='1|0' 

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

					else:
						
						if len(seq_h1) == 0: 

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='.|1' 

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

						else: 

							coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
							si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
							alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

							if alt1=='':

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='.|1' 

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

							else:

								if ref==alt1:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = hap2_dict_motif[reps]
									info['H2N'] = hap2_dict_number[reps]
			
									form='0|1' 

									VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

								else:

									if alt1==alt2:

										info=dict()
			
										info['END'] = reps[1]
										info['H1M'] = '.'
										info['H1N'] = '.'
										info['H2M'] = hap2_dict_motif[reps]
										info['H2N'] = hap2_dict_number[reps]
			
										form='1|1' 


										VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

									else:

										info=dict()
			
										info['END'] = reps[1]
										info['H1M'] = '.'
										info['H1N'] = '.'
										info['H2M'] = hap2_dict_motif[reps]
										info['H2N'] = hap2_dict_number[reps]


										form='1|2' 

										VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)
			

			elif reps not in ref_dict_number.keys() and reps in hap1_dict_number.keys() and reps in hap2_dict_number.keys(): 

				pos=reps[0]
				ref=reference_sequence[(reps[0]-1):reps[1]]

				seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
				coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
				si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
				alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

				seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 
				coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
				si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
				alt2=seq_h2[si_2:(ei_2+1)].replace('-','')


				if ref==alt1 and ref==alt2:

					continue

				elif ref==alt1 and ref != alt2:

					info=dict()

					info['END'] = reps[1]
					info['H1M'] = hap1_dict_motif[reps] 
					info['H1N'] = hap1_dict_number[reps] 
					info['H2M'] = hap2_dict_motif[reps]
					info['H2N'] = hap2_dict_number[reps]
									
					form='0|1' 

					VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)


				elif ref != alt1 and ref == alt2:

					info=dict()

					info['END'] = reps[1]
					info['H1M'] = hap1_dict_motif[reps]
					info['H1N'] = hap1_dict_number[reps]
					info['H2M'] = hap2_dict_motif[reps] 
					info['H2N'] = hap2_dict_number[reps] 
									
					form='1|0' 

					VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)


				elif ref != alt1 and ref != alt2:

					if alt1==alt2:

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
									
						form='1|1' 

						VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

					else:

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
									
						form='1|2' 

						VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)


			elif reps in ref_dict_number.keys() and reps not in hap1_dict_number.keys() and reps not in hap2_dict_number.keys(): 


				pos=reps[0]
				ref=reference_sequence[(reps[0]-1):reps[1]]

				seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
				seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 

				if len(seq_h1) == 0 and len(seq_h2) == 0: 

					continue


				elif len(seq_h1) != 0 and len(seq_h2) == 0: 

					coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
					si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
					alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

					if alt1=='':

						continue 

					else:

						if ref == alt1:

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.' 
							info['H1N'] = '.' 
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='0|.' 

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

						else:

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='1|.' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

				
				elif len(seq_h1) == 0 and len(seq_h2) != 0: 

					coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
					si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
					alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

					if alt2 == '':

						continue

					else:

						if ref == alt2:

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.' 
							info['H1N'] = '.'
							info['H2M'] = '.' 
							info['H2N'] = '.' 
									
							form='.|0' 

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

						else:

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='.|1' 

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

					
				elif len(seq_h1) != 0 and len(seq_h2) != 0: 

					coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
					si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
					alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

					coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
					si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
					alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

					if alt1 == '':

						if alt2 == '':

							continue 

						else:

							if ref==alt2:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='.|0' 

								VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

							else:

								info=dict()


								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='.|1' 

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

					else:

						if alt2 == '':

							if ref==alt1:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='0|.' 

								VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

							else:

								info=dict()


								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='1|.'  

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)


						else: 


							if ref == alt1 and ref == alt2: 

								continue

							elif ref == alt1 and ref != alt2:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.' 
								info['H1N'] = '.' 
								info['H2M'] = '.'
								info['H2N'] = '.'
			
								form='0|1' 

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

							elif ref != alt1 and ref == alt2:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='1|0' 

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

							elif ref != alt1 and ref != alt2:

								if alt1 != alt2 :

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = '.'
									info['H2N'] = '.'
			
									form='1|2' 

									VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)

								else:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = '.'
									info['H2N'] = '.'
			
									form='1|1' 

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)


			elif reps not in ref_dict_number.keys() and reps in hap1_dict_number.keys() and reps not in hap2_dict_number.keys(): 

				pos=reps[0]
				ref=reference_sequence[(reps[0]-1):reps[1]]

				seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 

				coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
				si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
				alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

				seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 

				if ref==alt1:

					if len(seq_h2) == 0: 

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = '.'
						info['H2N'] = '.'
									
						form='0|.' 

						VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

					else:

						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','')


						if alt2=='':

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='0|.' 

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

						else:

							if ref==alt2:

								continue

							else:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = '.'
								info['H2N'] = '.'
									
								form='0|1' 

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

				else:

					if len(seq_h2) == 0: 

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = '.'
						info['H2N'] = '.'
									
						form='1|.' 

						VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)


					else:


						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

						if alt2=='':

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='1|.' 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

						else:

							if ref==alt2:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = '.'
								info['H2N'] = '.'
									
								form='1|0' 

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)
							
							else:
								
								if alt1==alt2:

									info=dict()

									info['END'] = reps[1]
									info['H1M'] = hap1_dict_motif[reps]
									info['H1N'] = hap1_dict_number[reps]
									info['H2M'] = '.'
									info['H2N'] = '.'
									
									form='1|1' 

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)


								else:

									info=dict()

									info['END'] = reps[1]
									info['H1M'] = hap1_dict_motif[reps]
									info['H1N'] = hap1_dict_number[reps]
									info['H2M'] = '.'
									info['H2N'] = '.'
									
									form='1|2' 
							
									VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)

			else: 

				pos=reps[0]
				ref=reference_sequence[(reps[0]-1):reps[1]]

				seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])
				coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) 
				si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
				alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

				seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 

				if ref==alt2:

					if len(seq_h1) == 0: 

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = '.'
						info['H1N'] = '.'
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
									
						form='.|0' 

						VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

					else:

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

						if alt1=='':

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
									
							form='.|0' 

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out,processor)

						else:

							if ref==alt1:

								continue

							else:

								info=dict()
							
								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
									
								form='1|0' 

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

				else:


					if len(seq_h1) == 0: 

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = '.'
						info['H1N'] = '.'
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
									
						form='.|1' 

						VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

					else:

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) 
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

						if alt1=='':

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
									
							form='.|1' 


							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)

						else:

							if ref==alt1:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
								
								form='0|1' 


								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out,processor)
							
							else:

								if alt1==alt2:

									info=dict()

									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = hap2_dict_motif[reps]
									info['H2N'] = hap2_dict_number[reps]
									
									form='1|1' 

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out,processor)

								else:

									info=dict()

									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = hap2_dict_motif[reps]
									info['H2N'] = hap2_dict_number[reps]
									
									form='1|2' 
							
									VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out,processor)
