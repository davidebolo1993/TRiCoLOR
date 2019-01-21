import pysam
import datetime
import os
from collections import defaultdict
from operator import itemgetter

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
	RRM='##INFO=<ID=RRM,Number=1,Type=String,Description="Reference Repeated Motif">'
	RRN='##INFO=<ID=RRN,Number=1,Type=Integer,Description="Reference Repetitions Number">'
	ARM='##INFO=<ID=ARM,Number=1,Type=String,Description="Alternative-allele Repeated Motif">'
	ARN='##INFO=<ID=ARN,Number=1,Type=Integer,Description="Alternative-allele Repetitions Number">'

	#FORMAT field appear after INFO field

	FORMAT='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'

	classic_header='#CHROM' + '\t' + 'POS' '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' + '\t' + 'INFO' + '\t' + 'FORMAT' + '\t' + samplename.upper()

	with open(os.path.abspath(out + '/TRiCoLOR.vcf'), 'a') as vcfout:

		vcfout.write(vcf_format + '\n' + '##filedate=' + str(datetime.date.today()) +  '\n' + '##source=' + commandline + '\n') #filedate and source are not strictly required, but looks nice ! Adding command line to source can be also helpful

		for a,b in zip(chromosomes,sizes):

			vcfout.write('##contig=<ID='+str(a)+',length='+str(b)+'>'+'\n')


		vcfout.write(END + '\n' + RRM + '\n' + RRN + '\n' + ARM + '\n' + ARN + '\n')
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
	INFO_RRM=info['RRM']
	INFO_RRN=str(info['RRN'])
	INFO_ARM=info['ARM']
	INFO_ARN=str(info['ARN'])
	FORMAT=form


	with open(os.path.abspath(out + '/TRiCoLOR.vcf'), 'a') as vcfout:

		vcfout.write(CHROM + '\t' + POS + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER + '\t' + 'END='+INFO_END + ';'+ 'RRM='+INFO_RRM + ';' + 'RRN='+INFO_RRN + ';' + 'ARM='+INFO_ARM + ';' + 'ARN='+INFO_ARN + '\t' + GEN + '\t' + FORMAT + '\n')



def modifier(coordinates): #fast way to remove None and substitute with closest number in list

	
	coordinates=[el+1 if el is not None else el for el in coordinates] #get true coordinates
	start = next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:

			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates


def list_duplicates(list_of_seq):

	a_=defaultdict(list)

	for i,item in enumerate(list_of_seq):

		a_[item].append(i)

	return ((key,locs) for key,locs in a_.items() if len(locs) > 1)


def Modifier(list_of_coord,seq):


	coords_without_insertions=modifier(list_of_coord)

	where_dup=[]

	for dup in list_duplicates(coords_without_insertions):

		where_dup.append(dup)

	mod_dup=[]

	for dups in where_dup:

		dup_num=dups[0]
		to_add=1/(len(dups[1]))
		new_=[dup_num+(i*to_add) for i in range(len(dups[1]))]
		mod_dup.append((dup_num,new_))

	for i in range(len(mod_dup)):

		coords_without_insertions[min(where_dup[i][1]):max(where_dup[i][1])+1]=mod_dup[i][1]

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

	return coords_purified,NewSeq


def GetIndex(start, end, coordinates):

	si=[i for i,e in enumerate(coordinates) if e==start]
	ei=[i for i,e in enumerate(coordinates) if e==end]

	return si[0],ei[-1]



def VCF_writer(chromosome, reference_repetitions, reference_sequence, haplotype1_repetitions, sequence_coordinates_haplotype1, haplotype2_repetitions, sequence_coordinates_haplotype2, out):


	repref=reference_repetitions
	repsh1=list(haplotype1_repetitions)
	repsh2=list(haplotype2_repetitions)

	coord_h1,seq_h1=Modifier(sequence_coordinates_haplotype1[1],sequence_coordinates_haplotype1[0])
	coord_h2,seq_h2=Modifier(sequence_coordinates_haplotype2[1],sequence_coordinates_haplotype2[0])

	intersection=list(set(repref+ repsh1 + repsh2)) # get unique repetitions between the three lists
	sorted_intersection=sorted(intersection, key=itemgetter(1)) #sort repetitions


	for reps in sorted_intersection:

		if reps in repref and reps in repsh1 and reps in repsh2: #shared repetition, don't write to .vcf file as it is not a variant

			continue

		elif reps in repref and reps in repsh1 and reps not in repsh2: #repetition shared between reference and allele 1. Write variant for allele 2.

			pos=reps[1]
			ref=reference_sequence[(reps[1]-1):reps[2]]
			si,ei=GetIndex(reps[1],reps[2],coord_h2)
			alt=seq_h2[si:(ei+1)].replace('-','')

			if ref==alt: #if that repetition exists but was not called for some reason, exclude, as it is not a variant

				continue

			else:


				info=dict()
			
				info['END'] = reps[2]
				info['RRM'] = reps[0]
				info['RRN'] = reps[3]
				info['ARM'] = '.'
				info['ARN'] = '.'
			
				form='0|1' #second allele variant


				VCF_variantwriter(chromosome, pos, ref, alt, info, form, out)
		
		elif reps in repref and reps not in repsh1 and reps in repsh2: #repetition shared between reference and allele 2. Write variant for allele 1.


			pos=reps[1]
			ref=reference_sequence[(reps[1]-1):reps[2]]
			si,ei=GetIndex(reps[1],reps[2],coord_h1)
			alt=seq_h1[si:(ei+1)].replace('-','')
			
			if ref==alt: #if that repetition exists but was not called for some reason, exclude, as it is not a variant

				continue

			else:

				info=dict()
			
				info['END'] = reps[2]
				info['RRM'] = reps[0]
				info['RRN'] = reps[3]
				info['ARM'] = '.'
				info['ARN'] = '.'
			
				form='1|0' #first allele variant

				VCF_variantwriter(chromosome, pos, ref, alt, info, form, out)


		elif reps in repref and reps not in repsh1 and reps not in repsh2:

			pos=reps[1]
			ref=reference_sequence[(reps[1]-1):reps[2]]
			si_1,ei_1=GetIndex(reps[1],reps[2],coord_h1)
			si_2,ei_2=GetIndex(reps[1],reps[2],coord_h2)
			alt=seq_h1[si_1:(ei_1+1)].replace('-','') + ',' + seq_h2[si_2:(ei_2+1)].replace('-','')
			

			if ref==alt.split(',')[0] and ref != alt.split(',')[1]: #check if ref sequence is the same of the first allele

				info=dict()
			
				info['END'] = reps[2]
				info['RRM'] = reps[0]
				info['RRN'] = reps[3]
				info['ARM'] = '.'
				info['ARN'] = '.'
			
				form='0|1' #second allele variant

				VCF_variantwriter(chromosome, pos, ref, alt.split(',')[1], info, form, out)


			elif ref==alt.split(',')[1] and ref != alt.split(',')[0]: #check if ref sequence is the same of the second allele

				info=dict()
			
				info['END'] = reps[2]
				info['RRM'] = reps[0]
				info['RRN'] = reps[3]
				info['ARM'] = '.'
				info['ARN'] = '.'
			
				form='1|0' #first allele variant

				VCF_variantwriter(chromosome, pos, ref, alt.split(',')[0], info, form, out)


			elif ref==alt.split(',')[0] and ref==alt.split(',')[1]: #check if ref sequence is the same of both alleles

			
				continue


			else:


				if alt.split(',')[0] == alt.split(',')[1]:

					alt=alt.split(',')[0]


				info=dict()
			
				info['END'] = reps[2]
				info['RRM'] = reps[0]
				info['RRN'] = reps[3]
				info['ARM'] = '.'
				info['ARN'] = '.'
			
				form='1|1' #both variants

				VCF_variantwriter(chromosome, pos, ref, alt, info, form, out)


		elif reps not in repref and reps in repsh1 and reps not in repsh2: #repetition found in haplotype 1 only. Write variant for allele 1.


			pos=reps[1]
			ref=reference_sequence[(reps[1]-1):reps[2]]
			si,ei=GetIndex(reps[1],reps[2],coord_h1)
			alt=seq_h1[si:(ei+1)].replace('-','')

			if ref==alt: #check if ref sequence is the same of the allele


				continue

			else:
				
				info=dict()
			
				info['END'] = reps[2]
				info['RRM'] = '.'
				info['RRN'] = '.'
				info['ARM'] = reps[0]
				info['ARN'] = reps[3]
			
				form='1|0' #allele 1 is a variant
				
				VCF_variantwriter(chromosome, pos, ref, alt, info, form, out)


		elif reps not in repref and reps not in repsh1 and reps in repsh2: #repetition found in haplotype 2 only. Write variant for allele 2.


			pos=reps[1]
			ref=reference_sequence[(reps[1]-1):reps[2]]
			si,ei=GetIndex(reps[1],reps[2],coord_h2)
			alt=seq_h2[si:(ei+1)].replace('-','')

			if ref==alt:

				continue

			else:

				info=dict()
				info['END'] = reps[2]
				info['RRM'] = '.'
				info['RRN'] = '.'
				info['ARM'] = reps[0]
				info['ARN'] = reps[3]
			
				form='0|1' #allele 2 is a variant


				VCF_variantwriter(chromosome, pos, ref, alt, info, form, out)

		else: #repetition found in haplotype 1 and haplotype 2 but not in reference


			pos=reps[1]
			ref=reference_sequence[(reps[1]-1):reps[2]]
			si_1,ei_1=GetIndex(reps[1],reps[2],coord_h1)
			si_2,ei_2=GetIndex(reps[1],reps[2],coord_h2)
			alt=seq_h1[si_1:(ei_1+1)].replace('-','') + ',' + seq_h2[si_2:(ei_2+1)].replace('-','')
			info=dict()

			if ref==alt.split(',')[0] and ref != alt.split(',')[1]: #check if ref sequence is the same of the first allele

			
				info['END'] = reps[2]
				info['RRM'] = '.'
				info['RRN'] = '.'
				info['ARM'] = reps[0]
				info['ARN'] = reps[3]
			
				form='0|1' #allele 2 is variant

				VCF_variantwriter(chromosome, pos, ref, alt, info, form, out)

			elif ref==alt.split(',')[1] and ref != alt.split(',')[0]: #check if ref sequence is the same of the second allele


				info['END'] = reps[2]
				info['RRM'] = '.'
				info['RRN'] = '.'
				info['ARM'] = reps[0]
				info['ARN'] = reps[3]
			
				form='1|0' #allele 1 is variant

				VCF_variantwriter(chromosome, pos, ref, alt, info, form, out)

			elif ref==alt.split(',')[0] and ref==alt.split(',')[1]: #check if ref sequence is the same of both the alleles


				continue


			else:


				if alt.split(',')[0] == alt.split(',')[1]:

					alt=alt.split(',')[0]

				info['END'] = reps[2]
				info['RRM'] = '.'
				info['RRN'] = '.'
				info['ARM'] = reps[0]
				info['ARN'] = reps[3]
			
				form='1|1' #both alleles variant

				VCF_variantwriter(chromosome, pos, ref, alt, info, form, out)




