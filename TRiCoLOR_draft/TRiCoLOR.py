import sys
import os
import re
import pyfaidx
import pysam
import glob
import subprocess
import itertools
#import bisect
import csv
import math
from collections import defaultdict
from operator import itemgetter
from multiprocessing import Process
from shutil import which
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot
from RepFinder import *
from BamParser import *
import argparse



def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Tandem Repeats Caller fOr LOng Reads''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''')
	parser.add_argument('-g','--genome', help='reference genome', metavar='',required=True)
	parser.add_argument('-b', '--bed', help='.bed file generated during the sensing step or proprietary .bed file in the same format to use for STR searching', metavar='',required=True)
	parser.add_argument('-b1', '--bam1', help='haplotype-resolved bam1 file',metavar='',required=True)
	parser.add_argument('-b2', '--bam2', help='haplotype-resolved bam2 file',metavar='',required=True)
	parser.add_argument('-m','--motif', type=int, help='size of the motif of the repetition: use 0 for any. This parameter modify the search algorithm.',metavar='',default=0)
	parser.add_argument('-t', '--times', type=int, help='consencutive times a repetition must occur at least to be detected: use 0 for any. This parameter modify the search algorithm.',metavar='',default=3)
	parser.add_argument('-s', '--size', type=int, help='what size repetitions - even fuzzy ones - must have at least to be called',metavar='',required=True)	
	parser.add_argument('-mmi', '--mmiref', default=None, help='path to minimap2 .mmi chromosome-resolved references folder',metavar='')
	parser.add_argument('-o', '--output', help='path to where results will be saved',metavar='',required=True)
	args = parser.parse_args()


	#check the presence of needed external tools


	import timeit

	start_t=timeit.default_timer()

	external_tools=['samtools', 'minimap2', 'htsbox']

	for tools in external_tools:

		try:

			assert(which(tools) is not None)

		except:

			sys.exit(tools + ' was not found as an executable command. Install ' + tools + ' and restart TRiCoLOR')


	#check inputs validity

	#check if the genome file exists, is readable and is in .fasta format

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().split('>')[0] == '')

	except:

		sys.exit('The reference file does not exist, is not readable or is not in .fasta format')


	#check if the two .bam files exist, are readable and are valid .bam files


	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bam1)],stderr=open(os.devnull, 'wb'))

	except:

		sys.exit('The .bam1 file does not exist, is not readable or is not a valid .bam file')


	try:

		subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bam2)],stderr=open(os.devnull, 'wb'))
		
	except:

		sys.exit('The .bam2 file does not exist, is not readable or is not a valid .bam file')


	#look for the index of the .bam files


	if not os.path.exists(os.path.abspath(args.bam1 + '.bai')):

		print('Creating index for bam1, as no index was found in folder...')

		subprocess.call(['samtools', 'index', os.path.abspath(args.bam1)])

	if not os.path.exists(os.path.abspath(args.bam2 + '.bai')):

		print('Creating index for bam2, as no index was found in folder...')

		subprocess.call(['samtools', 'index', os.path.abspath(args.bam2)])

	#check write permissions on the path where the output will be saved

	try:

		assert(os.access(os.path.dirname(os.path.abspath(args.output)),os.W_OK))

	except:

		sys.exit('You do not have write permissions on the directory in which the directory result will be created: please use for which you have write permissions')


	if args.mmiref is not None:

		mmi_abspath=os.path.abspath(args.mmiref)

	else:

		mmi_abspath=os.path.dirname(os.path.abspath(args.genome))


	ref=pyfaidx.Fasta(args.genome)
	b_in=Bed_Reader(args.bed)
	chromosomes_seen=set()
	it_ = iter(b_in)

	for i in range(b_in.length()):

		chromosome, start, end=next(it_)

		if chromosome not in chromosomes_seen:

			chrom=ref[chromosome]
			ref_seq=chrom[:len(chrom)].seq
			chromosomes_seen.add(chromosome)

			#check if the .mmi reference exists, otherwise create it for the wanted chromosome

			if not os.path.exists(os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi')):

			#check if we have write permissions on reference directory

				try:

					assert(os.access(mmi_abspath, os.W_OK))

				except:

					sys.exit('You do not have write permissions on the reference folder: please create a new folder with a copy of the reference and use that folder')

				print('Creating a .mmi index for ' + chromosome + '...')

				try:

					subprocess.check_call(['samtools','faidx',os.path.abspath(args.genome),chromosome],stderr=open(os.devnull, 'wb'),stdout=os.path.abspath(mmi_abspath + '/' + chromosome + '.fa'))

				except:

					sys.exit('Something went wrong with the creation of the .mmi chromosome index, most likely your genome is in .fasta format but does not contain the information for chromosome ' + chromosome)


				subprocess.call(['minimap2','-d', os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi'),os.path.abspath(mmi_abspath + '/' + chromosome + '.fa')]) #must work if the previous step was done correctly
				
			mmi_ref=os.path.abspath(mmi_abspath + '/' + chromosome + '.mmi')

		print(start, end)

		#try:

		further = Ref_Repeats(ref_seq, chromosome, start, end, args.motif, args.times, args.size, args.output,i)

		if further:

			p1=Process(target=Haplo1_Repeats, args=(os.path.abspath(args.bam1), chromosome, start, end, args.motif, args.times, args.size, ref_seq, mmi_ref, args.output, i))
			p2=Process(target=Haplo2_Repeats, args=(os.path.abspath(args.bam2), chromosome, start, end, args.motif, args.times, args.size, ref_seq, mmi_ref, args.output, i))
			p1.start()
			p2.start()
			p1.join()
			p2.join()
			CompareTables(args.output,i)

		else:

			continue

		#except:

		#with open (os.path.abspath(args.output + '/Log.txt'),'a') as logout:

			#logout.write('Something wrong in region ' + str(start) + '-' + str(end) + '\n')

	CleanResults(args.output, os.path.abspath(args.bam1), os.path.abspath(args.bam2))


	end_t=timeit.default_timer()
	elapsed=end_t-start_t

	print('Done in', elapsed, 'seconds')



class Bed_Reader():

	def __init__(self,bedfile):

		self.bedfile=bedfile

	def __iter__(self):

		with open (self.bedfile, 'r') as bedin:

			for line in csv.reader(bedin, delimiter='\t'):

				if not line[0].startswith('#') and line !=[]: 

					if len(line) != 3:

						raise TypeError('Input must be a .bed file with 3 fields: chromosome, start and end')

					else:

						yield (line[0], int(line[1]), int(line[2]))

	def length(self):

		with open(self.bedfile, 'r') as bedin:

			size=sum(1 for _ in bedin if not _.startswith('#') and not _.strip()=='')

			return size


class EmptyTable():

	def __init__ (self, tablepath):

		self.tablepath=tablepath


	def write(self):

		if os.path.exists(os.path.abspath(self.tablepath)):

			return

		else:

			with open(self.tablepath, 'w') as refout:

				Empty=pd.DataFrame(columns=['Chromosome','Repetition Motif', 'Start', 'End', 'Repetitions Number'])
				Empty.to_csv(refout, index=False, sep='\t')



def isNaN(value):

	return value != value


def isEmpty(list_obj): #recursive function to check for empty list

    if list_obj == []:

        return True

    else:

        return all((isinstance(sli, list) and isEmpty(sli)) for sli in list_obj)


def Concat_Tables(list_of_paths):

	List_of_tables=[]

	for tab in list_of_paths:

		Tab=pd.read_csv(tab, sep='\t')
		List_of_tables.append(Tab)

	By_Row=pd.concat(List_of_tables, axis=0,ignore_index=True)

	return By_Row


def ResultsWriter_FromRef(chromosome,repetitions,realstart,out, iteration):

	seq=[el[0] for el in repetitions]
	start=[realstart+el[1] for el in repetitions]
	end=[realstart+el[2] for el in repetitions]
	rep=[el[3] for el in repetitions]
	chrom=[chromosome]*len(start)
	Table=pd.DataFrame({'Chromosome':chrom,'Repetition Motif':seq, 'Start':start,'End':end, 'Repetitions Number':rep},columns=['Chromosome','Repetition Motif', 'Start', 'End', 'Repetitions Number'])

	with open(os.path.abspath(out + '/' + str(iteration +1) + '_RepetitionsTable.tsv'), 'w') as refout:

		Table.to_csv(refout ,sep='\t',index=False)
	


def ResultsWriter_FromCons(chromosome,repetitions_with_coord, out, iteration): 

	seq=[el[0] for el in repetitions_with_coord]
	start=[el[1] for el in repetitions_with_coord]
	end=[el[2] for el in repetitions_with_coord]
	rep=[el[3] for el in repetitions_with_coord]
	chrom=[chromosome]*len(start)
	Table=pd.DataFrame({'Chromosome':chrom,'Repetition Motif':seq, 'Start':start,'End':end, 'Repetitions Number':rep},columns=['Chromosome','Repetition Motif', 'Start', 'End', 'Repetitions Number'])

	
	if os.path.exists(os.path.abspath(out + '/' + str(iteration+1) + '_RepetitionsTable.tsv')):

		with open(os.path.abspath(out + '/' + str(iteration+1) + '_RepetitionsTable.tsv'), 'a') as refout:

			Table.to_csv(refout ,sep='\t',index=False, header=False)

	else:

		with open(os.path.abspath(out + '/' + str(iteration +1) + '_RepetitionsTable.tsv'), 'a') as refout:

			Table.to_csv(refout ,sep='\t',index=False)



def Ref_Repeats(reference_seq, chromosome, start, end, kmer, times, size, out, iteration):

	out_=os.path.abspath(out+'/reference')

	if not os.path.exists(out_):

		os.makedirs(out_)

	wanted=reference_seq[start-1:end] #pyfaidx way to get start-end


	if 'N' in wanted: #region with ambiguous bases

		Table=EmptyTable(os.path.abspath(out_ + '/' + str(iteration +1) + '_RepetitionsTable.tsv'))
		Table.write()

		return False

	else:

		repetitions=list(RepeatsFinder(wanted,kmer,times))

		filtered=[rep for rep in repetitions if len(rep[0])*rep[3] >= size]

		if isEmpty(repetitions):

			Table=EmptyTable(os.path.abspath(out_ + '/' + str(iteration +1) + '_RepetitionsTable.tsv'))
			Table.write()

			return True

		else:

			ResultsWriter_FromRef(chromosome,filtered,start,out_, iteration)

			return True



def Haplo1_Repeats(bamfile1, chromosome, start, end, kmer, times, size ,ref_seq, mmi_ref, out, iteration):

	out_=os.path.abspath(out+'/haplotype1')

	if not os.path.exists(out_):

		os.makedirs(out_)

	seq,coord = Bamfile_Analyzer(bamfile1,chromosome,start,end)
	filseq,filcoord=InCommon(seq,coord)

	if isEmpty(filseq): #no informations for that region in this haplotype

		Table=EmptyTable(os.path.abspath(out_ +'/' + str(iteration +1) + '_RepetitionsTable.tsv'))
		Table.write()

		return

	else:

		Fasta_Generator(filseq,filcoord,out_)
		MA(out_,mmi_ref)
		consensus_bams=glob.glob(os.path.abspath(out_+'/*consensus.srt.bam'))

		for bam in consensus_bams:

			coords,seq=Get_Alignment_Positions(bam)

			if isEmpty(seq):

				Table=EmptyTable(os.path.abspath(out_ +'/' + str(iteration +1) + '_RepetitionsTable.tsv'))
				Table.write()

				continue

			else:

				repetitions=list(RepeatsFinder(seq,kmer,times))

				if isEmpty(repetitions): #no repetitions found in the region

					if bam != consensus_bams[-1]: #not the last one

						continue

					else:

						Table=EmptyTable(os.path.abspath(out_ +'/' + str(iteration +1) + '_RepetitionsTable.tsv'))
						Table.write()

				else:

					cor_coord_reps=corrector(ref_seq, seq, repetitions, coords, size, allowed=1) #probably an exception here is needed
					ResultsWriter_FromCons(chromosome, cor_coord_reps,out_,iteration)

		if len(consensus_bams) == 1:

			os.rename(consensus_bams[0], consensus_bams[0].replace('.'.join(consensus_bams[0].split('.',2)[:2]),os.path.abspath(out_+"/" +str(iteration+1))))
			os.remove(consensus_bams[0].replace('.bam', '.bam.bai'))

		else:

			with open(os.path.abspath(out_ + '/FileToMerge.txt'), 'a') as fin:

				for file in consensus_bams:

					fin.write(file + '\n')

			subprocess.call(['samtools', 'merge', '-b', os.path.abspath(out_ + '/FileToMerge.txt'), os.path.abspath(out_+'/'+str(iteration+1) + '.srt.bam')])
			os.remove(os.path.abspath(out_ + '/FileToMerge.txt'))

			for bams in consensus_bams:

				os.remove(bams)
				os.remove(bams.replace('.bam', '.bam.bai'))



def Haplo2_Repeats(bamfile2, chromosome, start, end, kmer, times, size, ref_seq, mmi_ref, out, iteration):

	out_=os.path.abspath(out+'/haplotype2')

	if not os.path.exists(out_):

		os.makedirs(out_)

	seq,coord = Bamfile_Analyzer(bamfile2,chromosome,start,end)
	filseq,filcoord=InCommon(seq,coord)

	if isEmpty(filseq): #no informations for that region in this haplotype

		Table=EmptyTable(os.path.abspath(out_ +'/' + str(iteration +1) + '_RepetitionsTable.tsv'))
		Table.write()

		return

	else:

		Fasta_Generator(filseq,filcoord,out_)
		MA(out_,mmi_ref)
		consensus_bams=glob.glob(os.path.abspath(out_+'/*consensus.srt.bam'))

		for bam in consensus_bams:

			coords,seq=Get_Alignment_Positions(bam)

			if isEmpty(seq):

				Table=EmptyTable(os.path.abspath(out_ +'/' + str(iteration +1) + '_RepetitionsTable.tsv'))
				Table.write()

				continue

			else:

				repetitions=list(RepeatsFinder(seq,kmer,times))

				if isEmpty(repetitions): #no repetitions found in the region

					if bam != consensus_bams[-1]: #not the last one

						continue

					else:

						Table=EmptyTable(os.path.abspath(out_ +'/' + str(iteration +1) + '_RepetitionsTable.tsv'))
						Table.write()

				else:

					cor_coord_reps=corrector(ref_seq, seq, repetitions, coords, size, allowed=1) #probably an exception here is needed
					ResultsWriter_FromCons(chromosome, cor_coord_reps,out_,iteration)

		if len(consensus_bams) == 1:

			os.rename(consensus_bams[0], consensus_bams[0].replace('.'.join(consensus_bams[0].split('.',2)[:2]),os.path.abspath(out_+"/" +str(iteration+1))))
			os.remove(consensus_bams[0].replace('.bam', '.bam.bai'))

		else:

			with open(os.path.abspath(out_ + '/FileToMerge.txt'), 'a') as fin:

				for file in consensus_bams:

					fin.write(file + '\n')

			subprocess.call(['samtools', 'merge', '-b', os.path.abspath(out_ + '/FileToMerge.txt'), os.path.abspath(out_+'/'+str(iteration+1) + '.srt.bam')])
			os.remove(os.path.abspath(out_ + '/FileToMerge.txt'))

			for bams in consensus_bams:

				os.remove(bams)
				os.remove(bams.replace('.bam', '.bam.bai'))


def CompareTables(out, iteration):

	out=os.path.abspath(out)
	out_=[os.path.abspath(out+j) for j in ['/reference', '/haplotype1', '/haplotype2']]

	Reference_Tsv=os.path.abspath(out_[0] + '/' + str(iteration +1) + '_RepetitionsTable.tsv')
	Haplo_One_Tsv=os.path.abspath(out_[1] + '/' + str(iteration +1) + '_RepetitionsTable.tsv')
	Haplo_Two_Tsv=os.path.abspath(out_[2] + '/' + str(iteration +1) + '_RepetitionsTable.tsv')

	Reference_Tab=pd.read_csv(Reference_Tsv, sep='\t')
	Haplo_One_Tab=pd.read_csv(Haplo_One_Tsv, sep='\t')
	Haplo_Two_Tab=pd.read_csv(Haplo_Two_Tsv, sep='\t')

	#Compare repetitions found in the reference and in the 2 haplotype-resolved .bam files for the region

	MergeHaploTab=Haplo_One_Tab.merge(Haplo_Two_Tab, how='outer', indicator=True)
	MergeHaploTab.sort_values(by=['Chromosome','Start'], inplace=True)
	MergeHaploTab.replace(to_replace={'_merge':{'left_only':'haplotype_1', 'right_only':'haplotype_2'}}, inplace=True)
	MergeHaploTab.rename({'_merge': 'haplo_differences'}, axis='columns', inplace=True) #rename to avoid problem with the next indicator
	MergedRefHaplo=MergeHaploTab=Reference_Tab.merge(MergeHaploTab, how='outer', indicator=True)
	MergedRefHaplo.sort_values(by=['Chromosome','Start'], inplace=True) #don't know if it's needed
	MergedRefHaplo.reset_index(inplace=True)


	FinalCol=[]


	if not '_merge' in MergedRefHaplo.columns:

		MergedRefHaplo.to_csv(os.path.abspath(out+'/Comparison_RepetitionsTable.tsv'),sep='\t',index=False)

	else:

		for i in range(len(MergedRefHaplo['_merge'])):

			if not isNaN(MergedRefHaplo['haplo_differences'][i]):

				if MergedRefHaplo['haplo_differences'][i] == 'both' and  MergedRefHaplo['_merge'][i] == 'both':

					FinalCol.append('all')

				elif MergedRefHaplo['haplo_differences'][i] != 'both' and  MergedRefHaplo['_merge'][i] == 'both':

					FinalCol.append('reference and '+str(MergedRefHaplo['haplo_differences'][i]))

				elif MergedRefHaplo['haplo_differences'][i] != 'both' and  MergedRefHaplo['_merge'][i] != 'both':

					if  MergedRefHaplo['_merge'][i] == 'left_only':

						FinalCol.append('reference')

					else:

						FinalCol.append(str(MergedRefHaplo['haplo_differences'][i]))

				elif MergedRefHaplo['haplo_differences'][i] == 'both' and  MergedRefHaplo['_merge'][i] != 'both':

					FinalCol.append('haplotype_1 and haplotype_2')

			else:

				FinalCol.append('reference')

	MergedRefHaplo.drop(columns=['index','_merge','haplo_differences'], inplace=True)
	MergedRefHaplo['Where']=FinalCol

	if os.path.exists(os.path.abspath(out+'/Comparison_RepetitionsTable.tsv')):

		with open(os.path.abspath(out+'/Comparison_RepetitionsTable.tsv'), 'a') as compareout:

			MergedRefHaplo.to_csv(compareout, sep='\t',index=False, header=False)

	else:

		with open(os.path.abspath(out+'/Comparison_RepetitionsTable.tsv'), 'a') as compareout:

			MergedRefHaplo.to_csv(compareout,sep='\t',index=False)



def CleanResults(out, bam1, bam2):

	#merge all the .srt.bam files for the two haplotypes

	out_=[os.path.abspath(out+j) for j in ['/reference', '/haplotype1', '/haplotype2']]

	Hap1_Bams=glob.glob(os.path.abspath(out_[1])+'/*.srt.bam')

	if isEmpty(Hap1_Bams):

		pass

	else:


		try:

			subprocess.check_call(['sh', '/home/bolognin/TRiCoLOR_py/Merging.sh', bam1, os.path.abspath(out_[1]),os.path.abspath(out_[1]+'/Haplotype1.merged')],stderr=open(os.devnull, 'wb'))

		except:

			print('Something wrong in final merging for haplotype 1. Aborted last step')

			return


	Hap2_Bams=glob.glob(os.path.abspath(out_[2])+'/*.srt.bam')


	if isEmpty(Hap2_Bams):

		pass

	else:

		try:

			subprocess.check_call(['sh', '/home/bolognin/TRiCoLOR_py/Merging.sh', bam2, os.path.abspath(out_[2]),os.path.abspath(out_[2]+'/Haplotype2.merged')],stderr=open(os.devnull, 'wb'))

		except:

			print('Something wrong in final merging for haplotype 2. Aborted last step')

			return


	RefTables=glob.glob(os.path.abspath(out_[0])+'/*.tsv')
	Hap1_Tables=glob.glob(os.path.abspath(out_[1])+'/*.tsv')
	Hap2_Tables=glob.glob(os.path.abspath(out_[2])+'/*.tsv')

	Table0=Concat_Tables(RefTables)
	Table1=Concat_Tables(Hap1_Tables)
	Table2=Concat_Tables(Hap2_Tables)


	with open(os.path.abspath(out_[0] + '/RepetitionsTable.tsv'), 'w') as refout:

		Table0.to_csv(refout, sep='\t', index=False)

	with open(os.path.abspath(out_[1] + '/RepetitionsTable.tsv'), 'w') as hap1out:

		Table1.to_csv(hap1out, sep='\t',index=False)

	with open(os.path.abspath(out_[2] + '/RepetitionsTable.tsv'), 'w') as hap2out:

		Table2.to_csv(hap2out, sep='\t',index=False)


	for tabs in RefTables:

		os.remove(tabs)

	for tabs in Hap1_Tables:

		os.remove(tabs)

	for tabs in Hap2_Tables:

		os.remove(tabs)



if __name__ == main():

	main()
