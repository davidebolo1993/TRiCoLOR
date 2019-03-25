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



def VCF_headerwriter(bamfile1, bamfile2, samplename, commandline, out):

	#Header built following instructions at 'https://samtools.github.io/hts-specs/VCFv4.2.pdf'

	bam1=pysam.AlignmentFile(bamfile1,'rb')
	header1=bam1.header
	bam2=pysam.AlignmentFile(bamfile2,'rb')
	header2=bam2.header

	chromosomes_info1=list(header1.items())[1][1]
	chromosomes_info2=list(header2.items())[1][1]

	chromosomes_info = list({x['SN']:x for x in chromosomes_info1 + chromosomes_info2}.values()) #just in case one header misses a contig information

	classic_chrs = ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y']]
	chromosomes=[]
	sizes=[]

	for infos in chromosomes_info:

		if infos['SN'] in classic_chrs:

			chromosomes.append(infos['SN'])
			sizes.append(infos['LN'])

	vcf_format='##fileformat=VCFv4.2'

	#INFO field appear after CONTIG field

	END='##INFO=<ID=END,Number=1,Type=Integer,Description="Repetition end">'
	H1M='##INFO=<ID=H1M,Number=.,Type=String,Description="Haplotype1 Repeated Motif">' #it is possible to have multiple repetition motifs for the same region, if they were splitted
	H1N='##INFO=<ID=H1N,Number=.,Type=String,Description="Haplotype1 Repetitions Number">' #it is possible to have multiple repetition numbers for the same region, if they were splitted
	H2M='##INFO=<ID=H2M,Number=.,Type=String,Description="Haplotype2 Repeated Motif">' #it is possible to have multiple repetition motifs for the same region, if they were splitted
	H2N='##INFO=<ID=H2N,Number=.,Type=String,Description="Haplotype2 Repetitions Number">' #it is possible to have multiple repetition numbers for the same region, if they were splitted

	#FORMAT field appear after INFO field

	FORMAT='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'

	classic_header='#CHROM' + '\t' + 'POS' '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' + '\t' + 'INFO' + '\t' + 'FORMAT' + '\t' + samplename.upper()

	with open(os.path.abspath(out + '/TRiCoLOR.vcf'), 'a') as vcfout:

		vcfout.write(vcf_format + '\n' + '##filedate=' + str(datetime.date.today()) +  '\n' + '##source=' + commandline + '\n') #filedate and source are not strictly required, but looks nice ! Adding command line to source can be also helpful

		for a,b in zip(chromosomes,sizes):

			vcfout.write('##contig=<ID='+str(a)+',length='+str(b)+'>'+'\n')


		vcfout.write(END + '\n' + H1M + '\n' + H1N + '\n' + H2M + '\n' + H2N + '\n')
		vcfout.write(FORMAT + '\n')
		vcfout.write('##SAMPLE=<ID=' + samplename +'>' + '\n' + classic_header + '\n')



def VCF_variantwriter(chrom, pos, ref, alt, info, form, out):


	ID='.' #this filed is by default missing
	FILTER='.' #this filed is by default missing
	QUAL='.' #this filed is by default missing
	GEN='GT'
	CHROM=chrom
	POS=str(pos)
	REF=ref
	ALT=alt
	INFO_END=str(info['END'])


	#deal with possible multiple numbers/motifs, as we can have overlapping ones

	
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


	with open(os.path.abspath(out + '/TRiCoLOR.vcf'), 'a') as vcfout:

		vcfout.write(CHROM + '\t' + POS + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER + '\t' + 'END='+INFO_END + ';'+ 'H1M='+INFO_H1M + ';' + 'H1N='+INFO_H1N + ';' + 'H2M='+INFO_H2M + ';' + 'H2N='+INFO_H2N + '\t' + GEN + '\t' + FORMAT + '\n')


def modifier(coordinates): #fast way to remove None and substitute with closest number in list

	
	coordinates=[el+1 if el is not None else el for el in coordinates] #get true coordinates
	start=next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:

			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates


def Modifier(list_of_coord,seq,reps):


	coords_without_insertions=modifier(list_of_coord)

	#Modify deletions
	
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


	if not reps[0] >= coords_purified[0]: #add fake positions that will be removed 


		number=reps[0]
		how_many=coords_purified[0]-reps[0]
		coords_purified = [number]*how_many + coords_purified
		NewSeq= '-'* how_many + NewSeq


	if not reps[1] <= coords_purified[-1]: #add fake positions that will be removed

		number=reps[1]
		how_many=reps[1] - coords_purified[-1]
		coords_purified = coords_purified + [number]*how_many 
		NewSeq= NewSeq + '-'* how_many


	return coords_purified,NewSeq



def Get_Seq_Pos(bamfilein,chromosome, start,end): #as the consensus sequence is supposed to generate just one sequence aligned to the reference, all other alignments are removes
	

	seq=[]
	coords=[]

	
	if os.stat(os.path.abspath(bamfilein)).st_size == 0: #if is empty, return empty seqs

		return seq,coords

	else:

		if start==end:

			end=start+1

		bamfile=pysam.AlignmentFile(bamfilein,'rb')

		for read in bamfile.fetch(chromosome, start-1, end-1): #fetch on zero-based, as clipped regions many not be in coordinates

			if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

				coords = read.get_reference_positions(full_length=True)
				seq=read.seq


	return seq,coords



def GetIndex(start, end, coordinates):

	si=bisect_left(coordinates, start)
	ei=bisect_right(coordinates, end) -1

	return si,ei


def recursive_merge(sorted_int, list_, i):

	new_=(min(list_, key=itemgetter(1)), max(list_,key=itemgetter(2))) #get extended range
	new_range=(new_[0][1], new_[-1][2])

	if i < len(sorted_int) -1:

		if sorted_int[i+1][1] <= new_range[1]:

			list_.append(sorted_int[i+1])
			recursive_merge(sorted_int, list_, i+1)




def Merger(sorted_int, refreps, h1reps, h2reps): #return non overlapping-ranges and dictionaries that keeps the original infos for the ranges.

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

			recursive_merge(sorted_int, list_, i+len(list_)-1) #check if new interval overlaps other intervals

			new_=(min(list_, key=itemgetter(1)), max(list_,key=itemgetter(2))) #get extended range
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




def VCF_writer(chromosome, reference_repetitions, reference_sequence, haplotype1_repetitions, bamfile1, haplotype2_repetitions, bamfile2, out):


	if not os.path.exists(os.path.abspath(bamfile1 + '.bai')):

		subprocess.call(['samtools', 'index', os.path.abspath(bamfile1)],stderr=open(os.devnull, 'wb'))
	

	if not os.path.exists(os.path.abspath(bamfile2 + '.bai')):

		subprocess.call(['samtools', 'index', os.path.abspath(bamfile2)],stderr=open(os.devnull, 'wb'))


	repref=reference_repetitions
	repsh1=list(haplotype1_repetitions)
	repsh2=list(haplotype2_repetitions)

	intersection=list(set(repref+ repsh1 + repsh2)) # get unique repetitions between the three lists

	if len(intersection) == 0:

		return

	else:

		sorted_intersection=sorted(intersection, key=itemgetter(1,2)) #sort repetitions by start and then by end

		sorted_ranges,ref_dict_number,ref_dict_motif,hap1_dict_number,hap1_dict_motif,hap2_dict_number,hap2_dict_motif=Merger(sorted_intersection, repref, repsh1, repsh2)

		for reps in sorted_ranges:

			if reps in ref_dict_number.keys() and reps in hap1_dict_number.keys() and reps in hap2_dict_number.keys():

				if ref_dict_number[reps] == hap1_dict_number[reps] and ref_dict_number[reps] == hap2_dict_number[reps]: #all same number

					if  ref_dict_motif[reps] == hap1_dict_motif[reps] and ref_dict_motif[reps] == hap2_dict_motif[reps]: #all same motif

						continue #write nothing, as there is no variant

					elif ref_dict_motif[reps] == hap1_dict_motif[reps] and ref_dict_motif[reps] != hap2_dict_motif[reps]: #hap2 different motif 

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','') #sequece where is supposed to be alteration

						if ref == alt2: #re-check if, for any reasons, the two sequences are the same. If so, skip to next iteration

							continue

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='0|1' #second allele variant. Unique variant. The sequence is known. 

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)


					elif ref_dict_motif[reps] != hap1_dict_motif[reps] and ref_dict_motif[reps] == hap2_dict_motif[reps]: #hap1 different motif 


						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1])
						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','') #sequece where is supposed to be alteration

						if ref == alt1: #re-check if, for any reasons, the two sequences are the same. If so, skip to next.

							continue

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='1|0' #first allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


					elif ref_dict_motif[reps] != hap1_dict_motif[reps] and ref_dict_motif[reps] != hap2_dict_motif[reps]: #haplos have both different motifs


						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables

						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)

						alt1=seq_h1[si_1:(ei_1+1)].replace('-','') #sequece where is supposed to be alteration
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','') #sequece where is supposed to be alteration


						if ref == alt1 and ref == alt2: #re-check if, for any reasons, the three sequences are the same. If so, skip to next iteration

							continue

						elif ref == alt1 and ref != alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = ref_dict_motif[reps] ########## hap1 and ref are the same
							info['H1N'] = ref_dict_number[reps] ########## hap1 and ref are the same
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='0|1' #second allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

						elif ref != alt1 and ref == alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = ref_dict_motif[reps] ########## hap2 and ref are the same
							info['H2N'] = ref_dict_number[reps] ########## hap2 and ref are the same
			
							form='1|0' #first allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


						elif ref != alt1 and ref != alt2:

							if alt1 != alt2 :

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|2' #both alleles variant. Double variant. The sequence is known

								VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|1' #both alleles variant. Unique variant. The sequence is known

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


				elif ref_dict_number[reps] == hap1_dict_number[reps] and ref_dict_number[reps] != hap2_dict_number[reps]: # h2 different number, already alteration

					if ref_dict_motif[reps] == hap1_dict_motif[reps]: #only h2 is a variant

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','') #sequece where is supposed to be alteration

						if ref == alt2: #re-check if, for any reasons, the two sequences are the same. If so, skip to next iteration

							continue

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='0|1' #second allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

					else: #motif is different, even h1 is a variant


						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) 

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables

						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)

						alt1=seq_h1[si_1:(ei_1+1)].replace('-','') #sequece where is supposed to be alteration
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','') #sequece where is supposed to be alteration


						if ref == alt1 and ref == alt2: #re-check if, for any reasons, the three sequences are the same. If so, skip to next.

							continue

						elif ref == alt1 and ref != alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = ref_dict_motif[reps] ### hap1 and ref are the same
							info['H1N'] = ref_dict_number[reps] ### hap1 and ref are the same
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='0|1' #second allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

						elif ref != alt1 and ref == alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = ref_dict_motif[reps] ### hap2 and ref are the same
							info['H2N'] = ref_dict_number[reps] ### hap2 and ref are the same
			
							form='1|0' #first allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

						
						elif ref != alt1 and ref != alt2:

							if alt1 != alt2 :

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|2' #both alleles variant. Double variant. The sequence is known

								VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|1' #both alleles variant. Unique variant. The sequence is known

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


				elif ref_dict_number[reps] != hap1_dict_number[reps] and ref_dict_number[reps] == hap2_dict_number[reps]: # h1 different number, already alteration

					if ref_dict_motif[reps] == hap2_dict_motif[reps]: #only h1 is a variant

						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1])
						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')


						if ref == alt1: #re-check if, for any reasons, the two sequences are the same. If so, skip to next.

							continue

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='1|0' #first allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

					else: #motif is different, even h2 is a variant


						pos=reps[0]
						ref=reference_sequence[(reps[0]-1):reps[1]]

						seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) 
						seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables


						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)

						alt1=seq_h1[si_1:(ei_1+1)].replace('-','') #sequece where is supposed to be alteration
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','') #sequece where is supposed to be alteration


						if ref == alt1 and ref == alt2: #re-check if, for any reasons, the three sequences are the same. If so, skip to next.

							continue

						elif ref == alt1 and ref != alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = ref_dict_motif[reps] ### hap1 and ref are the same
							info['H1N'] = ref_dict_number[reps] ### hap1 and ref are the same
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='0|1' #second allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

						elif ref != alt1 and ref == alt2:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = ref_dict_motif[reps] ### hap2 and ref are the same
							info['H2N'] = ref_dict_number[reps] ### hap2 and ref are the same
			
							form='1|0' #first allele variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

						elif ref != alt1 and ref != alt2:

							if alt1 != alt2 :

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|2' #both alleles variant. Double variant. The sequence is known

								VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='1|1' #both alleles variant. Single variant. The sequence is known

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


				else: # h1 and h2 different number, alterations

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1])
					seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])


					coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
					coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables

					si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1) 
					si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)

					alt1=seq_h1[si_1:(ei_1+1)].replace('-','') #sequece where is supposed to be alteration
					alt2=seq_h2[si_2:(ei_2+1)].replace('-','') #sequece where is supposed to be alteration


					if ref == alt1 and ref == alt2: #re-check if, for any reasons, the three sequences are the same. If so, skip to next.

						continue

					elif ref == alt1 and ref != alt2:

						info=dict()
			
						info['END'] = reps[1]
						info['H1M'] = ref_dict_motif[reps] ### hap1 and ref are the same
						info['H1N'] = ref_dict_number[reps] ### hap1 and ref are the same
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
			
						form='0|1' #second allele variant. Unique variant. The sequence is known

						VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

					elif ref != alt1 and ref == alt2:

						info=dict()
			
						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = ref_dict_motif[reps] ### hap2 and ref are the same
						info['H2N'] = ref_dict_number[reps] ### hap2 and ref are the same
			
						form='1|0' #first allele variant. Unique variant. The sequence is known

						VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

					elif ref != alt1 and ref != alt2:

						if alt1 != alt2 :

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='1|2' #both alleles variant. Double variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)

						else:

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='1|1' #both alleles variant. Unique variant. The sequence is known

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


			elif reps in ref_dict_number.keys() and reps in hap1_dict_number.keys() and reps not in hap2_dict_number.keys(): #range not present in hap2, can be a difference

				if ref_dict_number[reps] == hap1_dict_number[reps] and ref_dict_motif[reps] == hap1_dict_motif[reps]: #only hap2 can be variant

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) #this coordinates may not exist

					if len(seq_h2) == 0: #we don't have coverage on that region

						info=dict()
			
						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = '.'
						info['H2N'] = '.'
			
						form='0|.'  #first allele not variant. Second allele not known

						VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

					else: #we have coverage

						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)						
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','') #sequece where is supposed to be alteration

						if alt2 == '':

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = '.'
							info['H2N'] = '.'
			
							form='0|.' #first allele not variant. Second allele not known

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

						else:

							if ref==alt2:

								continue #skip to next iteration

							else:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = '.'
								info['H2N'] = '.'
			
								form='0|1' #first allele not variant. Second allele known

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

				
				else: #different number or motif for haplotype 1

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) #must exist
					coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
					si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
					alt1=seq_h1[si_1:(ei_1+1)].replace('-','')
					
					seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) #this may not exist

					if ref==alt1:

						if len(seq_h2) == 0: #we don't have coverage on that region

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = ref_dict_motif[reps] ############# ref and hap1 are the same
							info['H1N'] = ref_dict_number[reps] ############ ref and hap1 are the same
							info['H2M'] = '.'
							info['H2N'] = '.'
			
							form='0|.' #first allele not variant. Second allele not known

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

						else: #we have coverage

							coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
							si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
							alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

							if alt2 =='': 

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = ref_dict_motif[reps] ############# ref and hap1 are the same
								info['H1N'] = ref_dict_number[reps] ############ ref and hap1 are the same
								info['H2M'] = '.'
								info['H2N'] = '.'
			
								form='0|.' #first allele not variant. Second allele not known

								VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

							else:

								if ref==alt2:

									continue

								else:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = ref_dict_motif[reps] ############# ref and hap1 are the same
									info['H1N'] = ref_dict_number[reps] ############ ref and hap1 are the same
									info['H2M'] = '.'
									info['H2N'] = '.'
			
									form='0|1' #second allele variant. Second allele known

									VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

					else:

						if len(seq_h2) == 0: #we don't have coverage on that region

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = '.'
							info['H2N'] = '.'
			
							form='1|.' #first allele variant. Second allele not known

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

						else: #we have coverage

							coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
							si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
							alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

							if alt2=='': 

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = '.'
								info['H2N'] = '.'
			
								form='1|.' #first allele variant. Second allele not known


								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

							else:

								if ref==alt2:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = hap1_dict_motif[reps]
									info['H1N'] = hap1_dict_number[reps]
									info['H2M'] = '.'
									info['H2N'] = '.'
			
									form='1|0' #first allele variant. Second allele known

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

								else:

									if alt1==alt2:

										info=dict()
			
										info['END'] = reps[1]
										info['H1M'] = hap1_dict_motif[reps]
										info['H1N'] = hap1_dict_number[reps]
										info['H2M'] = '.'
										info['H2N'] = '.'
			
										form='1|1' #both alleles variant. Unique variant


										VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

									else:

										info=dict()
			
										info['END'] = reps[1]
										info['H1M'] = hap1_dict_motif[reps]
										info['H1N'] = hap1_dict_number[reps]
										info['H2M'] = '.'
										info['H2N'] = '.'
			
										form='1|2' #both alleles are variant. Double variant


										VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)


			elif reps in ref_dict_number.keys() and reps not in hap1_dict_number.keys() and reps in hap2_dict_number.keys(): #range not present in hap1, can be a difference

				if ref_dict_number[reps] == hap2_dict_number[reps] and ref_dict_motif[reps] == hap2_dict_motif[reps]: #only hap1 can be variant

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1])

					if len(seq_h1) == 0: #we don't have coverage on that region

						info=dict()
			
						info['END'] = reps[1]
						info['H1M'] = '.'
						info['H1N'] = '.'
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
			
						form='.|0' #second allele not variant. first allele not known

						VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

					else: #we have coverage

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')


						if alt1=='':

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='.|0' #second allele not variant. first allele not known

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

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
			
								form='1|0' #first allele variant. Sequence known

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)
				
				else: #different number or motif for haplotype 2

					pos=reps[0]
					ref=reference_sequence[(reps[0]-1):reps[1]]

					seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])
					coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
					si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
					alt2=seq_h2[si_2:(ei_2+1)].replace('-','')
					
					seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) #this may not exist

					if ref==alt2:

						if len(seq_h1) == 0: #we don't have coverage on that region

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = ref_dict_motif[reps] ############# ref and hap2 are the same
							info['H2N'] = ref_dict_number[reps] ############# ref and hap2 are the same
			
							form='.|0' #second allele not variant. First allele not known

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

						else: #we have coverage

							coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
							si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
							alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

							if alt1=='': #start coordinate or end coordinate not present, same as we don't have coverage

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = ref_dict_motif[reps] ############# ref and hap2 are the same
								info['H2N'] = ref_dict_number[reps] ############# ref and hap2 are the same
			
								form='.|0' #second allele is not variant. First allele not known

								VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

							else:

								if ref==alt1:

									continue

								else:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = ref_dict_motif[reps] ############# ref and hap2 are the same
									info['H2N'] = ref_dict_number[reps] ############# ref and hap2 are the same
			
									form='1|0' #first allele variant.

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

					else:
						
						if len(seq_h1) == 0: #we don't have coverage on that region

							info=dict()
			
							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
			
							form='.|1' #second allele is variant. First allele not known

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

						else: #we have coverage

							coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
							si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
							alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

							if alt1=='':

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
			
								form='.|1' #second allele is variant. First allele not known

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

							else:

								if ref==alt1:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = hap2_dict_motif[reps]
									info['H2N'] = hap2_dict_number[reps]
			
									form='0|1' #second allele is variant. First allele known

									VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

								else:

									if alt1==alt2:

										info=dict()
			
										info['END'] = reps[1]
										info['H1M'] = '.'
										info['H1N'] = '.'
										info['H2M'] = hap2_dict_motif[reps]
										info['H2N'] = hap2_dict_number[reps]
			
										form='1|1' #both alleles variant. Unique variant


										VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

									else:

										info=dict()
			
										info['END'] = reps[1]
										info['H1M'] = '.'
										info['H1N'] = '.'
										info['H2M'] = hap2_dict_motif[reps]
										info['H2N'] = hap2_dict_number[reps]


										form='1|2' #both alleles variant. Double variant

										VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)
			

			elif reps not in ref_dict_number.keys() and reps in hap1_dict_number.keys() and reps in hap2_dict_number.keys(): #range not present in ref, can be a difference

				pos=reps[0]
				ref=reference_sequence[(reps[0]-1):reps[1]]

				seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) #this coordinates may not exist
				coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
				si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
				alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

				seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) #this coordinates may not exist
				coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
				si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
				alt2=seq_h2[si_2:(ei_2+1)].replace('-','')


				if ref==alt1 and ref==alt2:

					continue

				elif ref==alt1 and ref != alt2:

					info=dict()

					info['END'] = reps[1]
					info['H1M'] = hap1_dict_motif[reps] #cannot change to reference dict, as we have not such information
					info['H1N'] = hap1_dict_number[reps] #cannot change to reference dict, as we have not such information
					info['H2M'] = hap2_dict_motif[reps]
					info['H2N'] = hap2_dict_number[reps]
									
					form='0|1' #second allele variant.

					VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)


				elif ref != alt1 and ref == alt2:

					info=dict()

					info['END'] = reps[1]
					info['H1M'] = hap1_dict_motif[reps]
					info['H1N'] = hap1_dict_number[reps]
					info['H2M'] = hap2_dict_motif[reps] #cannot change to reference dict, as we have not such information
					info['H2N'] = hap2_dict_number[reps] #cannot change to reference dict, as we have not such information
									
					form='1|0' #first allele variant.

					VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


				elif ref != alt1 and ref != alt2:

					if alt1==alt2:

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
									
						form='1|1' #both alleles variants. Unique variant

						VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

					else:

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
									
						form='1|2' #both alleles variants. Double variant

						VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)


			elif reps in ref_dict_number.keys() and reps not in hap1_dict_number.keys() and reps not in hap2_dict_number.keys(): #range present just in reference


				pos=reps[0]
				ref=reference_sequence[(reps[0]-1):reps[1]]

				seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) #this coordinates may not exist
				seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) #this coordinates may not exist

				if len(seq_h1) == 0 and len(seq_h2) == 0: #we don't have coverage for the two regions

					continue


				elif len(seq_h1) != 0 and len(seq_h2) == 0: #we don't have coverage for hap2

					coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
					si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
					alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

					if alt1=='':

						continue 

					else:

						if ref == alt1:

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.' ###
							info['H1N'] = '.' ###
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='0|.' #first allele not variant. Second allele not known

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

						else:

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='1|.' #first allele variant. Second allele not known

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

				
				elif len(seq_h1) == 0 and len(seq_h2) != 0: #we don't have coverage for hap1

					coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
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
							info['H2M'] = '.' ###
							info['H2N'] = '.' ###
									
							form='.|0' #second allele not variant. First allele not known

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

						else:

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='.|1' #second allele variant. First allele not known

							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

					
				elif len(seq_h1) != 0 and len(seq_h2) != 0: #coverage for both

					coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
					si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
					alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

					coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
					si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
					alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

					if alt1 == '':

						if alt2 == '':

							continue #range outside both haplotypes

						else:

							if ref==alt2:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='.|0' #second allele not variant. first one not known

								VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

							else:

								info=dict()


								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='.|1' #second allele variant. first one not known

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

					else:

						if alt2 == '':

							if ref==alt1:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='0|.' #first allele not variant. Second not known

								VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

							else:

								info=dict()


								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='1|.'  #first allele variant. Second not known

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


						else: #no haplotype out of range


							if ref == alt1 and ref == alt2: 

								continue

							elif ref == alt1 and ref != alt2:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.' ### 
								info['H1N'] = '.' ### 
								info['H2M'] = '.'
								info['H2N'] = '.'
			
								form='0|1' #second allele variant.

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

							elif ref != alt1 and ref == alt2:

								info=dict()
			
								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = '.' 
								info['H2N'] = '.' 
			
								form='1|0' #first allele variant.

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

							elif ref != alt1 and ref != alt2:

								if alt1 != alt2 :

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = '.'
									info['H2N'] = '.'
			
									form='1|2' #both alleles variant. Double variant

									VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)

								else:

									info=dict()
			
									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = '.'
									info['H2N'] = '.'
			
									form='1|1' #both alleles variant. Unique variant

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


			elif reps not in ref_dict_number.keys() and reps in hap1_dict_number.keys() and reps not in hap2_dict_number.keys(): #range present just in haplotype 1

				pos=reps[0]
				ref=reference_sequence[(reps[0]-1):reps[1]]

				seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) #this coordinates must exist

				coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
				si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
				alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

				seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1]) #this coordinates may not exist

				if ref==alt1:

					if len(seq_h2) == 0: #no coverage for the second allele

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = '.'
						info['H2N'] = '.'
									
						form='0|.' #first allele not variant. Second allele not known

						VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

					else:

						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','')


						if alt2=='':

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='0|.' #first allele not variant. Second allele not known

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

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
									
								form='0|1' #second allele variant

								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

				else:

					if len(seq_h2) == 0: #no coverage for the second allele

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = hap1_dict_motif[reps]
						info['H1N'] = hap1_dict_number[reps]
						info['H2M'] = '.'
						info['H2N'] = '.'
									
						form='1|.' #first allele variant. Second one not known 

						VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


					else:


						coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
						si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
						alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

						if alt2=='':

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = hap1_dict_motif[reps]
							info['H1N'] = hap1_dict_number[reps]
							info['H2M'] = '.'
							info['H2N'] = '.'
									
							form='1|.' #first allele variant. Second one not known 

							VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

						else:

							if ref==alt2:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = hap1_dict_motif[reps]
								info['H1N'] = hap1_dict_number[reps]
								info['H2M'] = '.'
								info['H2N'] = '.'
									
								form='1|0' #first allele variant.

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)
							
							else:
								
								if alt1==alt2:

									info=dict()

									info['END'] = reps[1]
									info['H1M'] = hap1_dict_motif[reps]
									info['H1N'] = hap1_dict_number[reps]
									info['H2M'] = '.'
									info['H2N'] = '.'
									
									form='1|1' #both alleles variant. Unique variant

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)


								else:

									info=dict()

									info['END'] = reps[1]
									info['H1M'] = hap1_dict_motif[reps]
									info['H1N'] = hap1_dict_number[reps]
									info['H2M'] = '.'
									info['H2N'] = '.'
									
									form='1|2' #both alleles variant. Double variant
							
									VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)

			else: #range present just in haplotype 2

				pos=reps[0]
				ref=reference_sequence[(reps[0]-1):reps[1]]

				seq_h2,coord_h2=Get_Seq_Pos(bamfile2,chromosome, reps[0], reps[1])#this coordinates must exist
				coord_h2,seq_h2=Modifier(coord_h2,seq_h2,reps) #overwrite previous variables
				si_2,ei_2=GetIndex(reps[0],reps[1],coord_h2)
				alt2=seq_h2[si_2:(ei_2+1)].replace('-','')

				seq_h1,coord_h1=Get_Seq_Pos(bamfile1,chromosome, reps[0], reps[1]) #this coordinates may not exist

				if ref==alt2:

					if len(seq_h1) == 0: #no coverage for the second allele

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = '.'
						info['H1N'] = '.'
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
									
						form='.|0' #second allele not a variant. First one not known

						VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

					else:

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
						si_1,ei_1=GetIndex(reps[0],reps[1],coord_h1)
						alt1=seq_h1[si_1:(ei_1+1)].replace('-','')

						if alt1=='':

							info=dict()

							info['END'] = reps[1]
							info['H1M'] = '.'
							info['H1N'] = '.'
							info['H2M'] = hap2_dict_motif[reps]
							info['H2N'] = hap2_dict_number[reps]
									
							form='.|0' #second allele not a variant. First one not known

							VCF_variantwriter(chromosome, pos, ref, '.', info, form, out)

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
									
								form='1|0' #first allele variant

								VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

				else:


					if len(seq_h1) == 0: #no coverage for the second allele

						info=dict()

						info['END'] = reps[1]
						info['H1M'] = '.'
						info['H1N'] = '.'
						info['H2M'] = hap2_dict_motif[reps]
						info['H2N'] = hap2_dict_number[reps]
									
						form='.|1' #second allele variant. First one not known

						VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

					else:

						coord_h1,seq_h1=Modifier(coord_h1,seq_h1,reps) #overwrite previous variables
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


							VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)

						else:

							if ref==alt1:

								info=dict()

								info['END'] = reps[1]
								info['H1M'] = '.'
								info['H1N'] = '.'
								info['H2M'] = hap2_dict_motif[reps]
								info['H2N'] = hap2_dict_number[reps]
								
								form='0|1' #second allele variant.


								VCF_variantwriter(chromosome, pos, ref, alt2, info, form, out)
							
							else:

								if alt1==alt2:

									info=dict()

									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = hap2_dict_motif[reps]
									info['H2N'] = hap2_dict_number[reps]
									
									form='1|1' #both alleles variant. Unique variant

									VCF_variantwriter(chromosome, pos, ref, alt1, info, form, out)

								else:

									info=dict()

									info['END'] = reps[1]
									info['H1M'] = '.'
									info['H1N'] = '.'
									info['H2M'] = hap2_dict_motif[reps]
									info['H2N'] = hap2_dict_number[reps]
									
									form='1|2' #both alleles variant. Double variant
							
									VCF_variantwriter(chromosome, pos, ref, alt1 + ',' + alt2, info, form, out)