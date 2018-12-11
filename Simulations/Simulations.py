#!/usr/bin/python

import subprocess
import random
import glob
from itertools import chain
import os
import pandas as pd
import scipy
import numpy as np
import matplotlib.pyplot as plt


def Simulate(bed, reference, number_of_simulations, sim_type, model_qc, accuracy, out):


	for i in range(1,number_of_simulations+1):

		acceptable=5 #when looking for the wanted repetitions, accepts a repetitions that will have start or end close to the supposed one (in case of deletions, start can be after the true one by the aligner)

		#chromosome chrX-chrY are not taken into account

		number_of_chr = 22
		chromosomes=['chr%s'%i for i in chain(range(1, number_of_chr+1))]
		chromosome=random.sample(chromosomes,k=1)[0]
		point=random.randint(16000000,20000000) #4Mb region that should be common to all the chromosomes and avoid telomeres/centromeres

		if not os.path.exists(os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa')): #checks for the presence of a .mmi ref for the chromosome in the reference folder. If not present, creates it. Will save time during alignments.

			with open(os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'),'w') as fout: 

				subprocess.call(['samtools', 'faidx', reference, chromosome], stdout=fout)

			subprocess.call(['minimap2', '-d', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.mmi'), os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa')])

		new_mmi=os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.mmi')

		if sim_type=='contraction': 

			label1='ChrSim_h1' #label for haplotype 1
			label2='ChrSim_h2' #label for haplotype 2

			folder=os.path.abspath(out + '/Contractions/Sim'+str(i)) #output for contractions

			#simulate haplotype 1, with alterations

			subprocess.call(['python', '/home/bolognin/TRiCoLOR_py/ModifyTR.py', '--reference',reference, '--bed', bed, '--chromosome', chromosome, '--mode', 'contraction','--position', str(point), '--size', '5', '--label',label1, '--output',folder])
			
			#simulate haplotype 2, without alterations
			
			subprocess.call(['python', '/home/bolognin/TRiCoLOR_py/ModifyTR.py', '--reference',reference, '--bed', bed, '--chromosome', chromosome, '--mode', 'contraction','--position', str(point), '--size', '0', '--label',label2, '--output',folder])

			#resize haplotypes to the region in which contraction is: save time when simulate reads in the region

			with open(os.path.abspath(folder+'/'+label1+'.fa')) as fin1, open(os.path.abspath(folder+'/'+label2+'.fa')) as fin2, open(os.path.abspath(folder+'/'+label1.replace('ChrSim_h1','RegSim_h1')+'.fa'),'w') as fout1,open(os.path.abspath(folder+'/'+label2.replace('ChrSim_h2','RegSim_h2')+'.fa'),'w') as fout2:

				#get informations from fasta header

				header1=fin1.readline().strip().split('>')[1]
				header2=fin2.readline().strip().split('>')[1]
				kmer=len(header1.split('_',3)[2].split(':')[1])
				times1= int(header1.split('_',5)[4].split(':')[1])
				times2=int(header2.split('_',5)[4].split(':')[1])
				true_point=int(header1.split(':',1)[1].split('_',1)[0])

				subprocess.call(['samtools', 'faidx', os.path.abspath(folder+'/'+label1+'.fa'), header1+':'+str(true_point-51000)+'-'+str(true_point+51000)],stdout=fout1) #rescale region, simulation will be faster
				subprocess.call(['samtools', 'faidx', os.path.abspath(folder+'/'+label2+'.fa'), header2+':'+str(true_point-51000)+'-'+str(true_point+51000)],stdout=fout2) #rescale region, simulation will be faster


			os.remove(os.path.abspath(folder+'/'+label1+'.fa'))
			os.remove(os.path.abspath(folder+'/'+label2+'.fa'))
			os.remove(os.path.abspath(folder+'/'+label1+'.fa.fai'))
			os.remove(os.path.abspath(folder+'/'+label2+'.fa.fai'))

			#simualte ONT reads in region for hap1 and hap2

			os.chdir(folder)

			#reads for haplotype 1

			subprocess.call(['pbsim', '--prefix', 'Simh1', '--depth', '20', '--difference-ratio', '30:30:40','--data-type', 'CLR', '--length-mean', '8000', '--accuracy-mean', str(accuracy), '--model_qc', model_qc, os.path.abspath(folder+'/'+'RegSim_h1.fa')])
			
			#reads for haplotype 2

			subprocess.call(['pbsim', '--prefix', 'Simh2', '--depth', '20', '--difference-ratio', '30:30:40','--data-type', 'CLR', '--length-mean', '8000', '--accuracy-mean', str(accuracy), '--model_qc', model_qc, os.path.abspath(folder+'/'+'RegSim_h2.fa')])

			os.remove('Simh1_0001.ref')
			os.remove('Simh1_0001.maf')
			os.remove('Simh2_0001.ref')
			os.remove('Simh2_0001.maf')

			with open('Simh1.sam', 'w') as fh1:

				subprocess.call(['minimap2', '-ax', 'map-ont', new_mmi, 'Simh1_0001.fastq'], stdout=fh1)

			os.remove('Simh1_0001.fastq')

			with open(os.path.abspath('Simh2.sam'), 'w') as fh2:

				subprocess.call(['minimap2', '-ax', 'map-ont', new_mmi,'Simh2_0001.fastq'], stdout=fh2)

			os.remove(os.path.abspath('Simh2_0001.fastq'))

			with open(os.path.abspath('Simh1.bam'), 'w') as fh1:

				subprocess.call(['htsbox', 'samview', '-bS','Simh1.sam'], stdout=fh1)

			with open(os.path.abspath('Simh2.bam'), 'w') as fh2:

				subprocess.call(['htsbox', 'samview', '-bS','Simh2.sam'], stdout=fh2)

			with open(os.path.abspath('Simh1.srt.bam'), 'w') as fh1:

				subprocess.call(['samtools', 'sort', 'Simh1.bam'], stdout=fh1)

			with open(os.path.abspath('Simh2.srt.bam'), 'w') as fh2:

				subprocess.call(['samtools', 'sort', 'Simh2.bam'], stdout=fh2)

			os.remove('Simh1.sam')
			os.remove('Simh2.sam')
			os.remove('Simh1.bam')
			os.remove('Simh2.bam')
			subprocess.call(['samtools', 'index', 'Simh1.srt.bam']) #indexed haplotype1
			subprocess.call(['samtools', 'index', 'Simh2.srt.bam']) #indexed haplotype2

			#if any error occur when calling TRiCoLOR (.bams that are not readable for some reason) skip to next iteration (does not seem to be necessary)


			with open ('Sim.bed','w') as bedout:

				bedout.write(chromosome + '\t' + str(true_point-1000) + '\t' +  str(true_point+1000) + '\n')


			try:

				subprocess.check_call(['python','/home/bolognin/TRiCoLOR_py/TRiCoLOR.py','-g', reference, '-b','Sim.bed', '-b1', 'Simh1.srt.bam', '-b2', 'Simh2.srt.bam', '-m', str(kmer), '-t','4', '-o',os.path.abspath(out + '/Contractions/Sim'+str(i)+'/Res'+str(i)), '-mmi', new_mmi])
			
			except:

				continue


			#can we predict the exact number of tandem repeat contraction ?

			reference_repetitions=int(header1.split('_',4)[3].split(':')[1])
			repetition_type=header1.split('_',3)[2].split(':')[1]

			Hap1_Tab=pd.read_csv(os.path.abspath(out +'/Contractions/Sim'+str(i)+'/Res'+str(i)+'/haplotype1/RepetitionsTable.tsv'), sep='\t')
			Hap2_Tab=pd.read_csv(os.path.abspath(out +'/Contractions/Sim'+str(i)+'/Res'+str(i)+'/haplotype2/RepetitionsTable.tsv'), sep='\t')

			Hap1=[]
			Hap2=[]

			#look for the wanted repetition in range

			for l in range(len(Hap1_Tab['Start'])):

				if abs(Hap1_Tab['Start'][l] - true_point) <= acceptable or abs (Hap1_Tab['End'][l] -  (true_point+len(repetition_type)*reference_repetitions)) <= acceptable:

					Hap1.append((Hap1_Tab['Repetition Motif'][l],Hap1_Tab['Repetitions Number'][l]))


			for l in range(len(Hap2_Tab['Start'])):

				if abs(Hap2_Tab['Start'][l] - true_point) <= acceptable or abs (Hap2_Tab['End'][l] -  (true_point+len(repetition_type)*reference_repetitions)) <= acceptable:

					Hap2.append((Hap2_Tab['Repetition Motif'][l],Hap2_Tab['Repetitions Number'][l]))


			#sum the occurences, even if is only one: something left from the first time we simulate looking for general expansions and contractions

			Rep_Hap1=sum(n for _, n in Hap1)
			Rep_Hap2=sum(n for _, n in Hap2)


			Hap1_Res=[]
			Hap2_Res=[]
			Hap1_Res.append('Reference:'+str(reference_repetitions)+'-Haplo:'+str(Rep_Hap1)) 
			Hap2_Res.append('Reference:'+str(reference_repetitions)+'-Haplo:'+str(Rep_Hap2)) 
			Hap1_Type=['Reference:'+repetition_type+'-Haplo:'+Hap1[0][0]]
			Hap2_Type=['Reference:'+repetition_type+'-Haplo:'+Hap2[0][0]]
			TableRes=pd.DataFrame({'Hap1_Num':Hap1_Res,'Hap2_Num':Hap2_Res, 'Hap1_Core': Hap1_Type, 'Hap2_Core':Hap2_Type})
			TableRes.to_csv('Sim_'+str(i)+'_Results.tsv',sep='\t',index=False)


			os.chdir(os.path.abspath(out))


		if sim_type=='expansion': 

			label1='ChrSim_h1' #label for haplotype 1
			label2='ChrSim_h2' #label for haplotype 2

			folder=os.path.abspath(out + '/Expansions/Sim'+str(i)) #output for contractions

			#simulate haplotype 1, with alterations

			subprocess.call(['python', '/home/bolognin/TRiCoLOR_py/ModifyTR.py', '--reference',reference, '--bed', bed, '--chromosome', chromosome, '--mode', 'expansion','--position', str(point), '--size', '5', '--label',label1, '--output',folder])
			
			#simulate haplotype 2, without alterations
			
			subprocess.call(['python', '/home/bolognin/TRiCoLOR_py/ModifyTR.py', '--reference',reference, '--bed', bed, '--chromosome', chromosome, '--mode', 'expansion','--position', str(point), '--size', '0', '--label',label2, '--output',folder])

			#resize haplotypes to the region in which contraction is: save time when simulate reads in the region

			with open(os.path.abspath(folder+'/'+label1+'.fa')) as fin1, open(os.path.abspath(folder+'/'+label2+'.fa')) as fin2, open(os.path.abspath(folder+'/'+label1.replace('ChrSim_h1','RegSim_h1')+'.fa'),'w') as fout1,open(os.path.abspath(folder+'/'+label2.replace('ChrSim_h2','RegSim_h2')+'.fa'),'w') as fout2:

				#get informations from fasta header

				header1=fin1.readline().strip().split('>')[1]
				header2=fin2.readline().strip().split('>')[1]
				kmer=len(header1.split('_',3)[2].split(':')[1])
				times1= int(header1.split('_',5)[4].split(':')[1])
				times2=int(header2.split('_',5)[4].split(':')[1])
				true_point=int(header1.split(':',1)[1].split('_',1)[0])

				subprocess.call(['samtools', 'faidx', os.path.abspath(folder+'/'+label1+'.fa'), header1+':'+str(true_point-51000)+'-'+str(true_point+51000)],stdout=fout1) #rescale region, simulation will be faster
				subprocess.call(['samtools', 'faidx', os.path.abspath(folder+'/'+label2+'.fa'), header2+':'+str(true_point-51000)+'-'+str(true_point+51000)],stdout=fout2) #rescale region, simulation will be faster


			os.remove(os.path.abspath(folder+'/'+label1+'.fa'))
			os.remove(os.path.abspath(folder+'/'+label2+'.fa'))
			os.remove(os.path.abspath(folder+'/'+label1+'.fa.fai'))
			os.remove(os.path.abspath(folder+'/'+label2+'.fa.fai'))

			#simualte ONT reads in region for hap1 and hap2

			os.chdir(folder)

			#reads for haplotype 1

			subprocess.call(['pbsim', '--prefix', 'Simh1', '--depth', '20', '--difference-ratio', '30:30:40','--data-type', 'CLR', '--length-mean', '8000', '--accuracy-mean', str(accuracy), '--model_qc', model_qc, os.path.abspath(folder+'/'+'RegSim_h1.fa')])
			
			#reads for haplotype 2

			subprocess.call(['pbsim', '--prefix', 'Simh2', '--depth', '20', '--difference-ratio', '30:30:40','--data-type', 'CLR', '--length-mean', '8000', '--accuracy-mean', str(accuracy), '--model_qc', model_qc, os.path.abspath(folder+'/'+'RegSim_h2.fa')])

			os.remove('Simh1_0001.ref')
			os.remove('Simh1_0001.maf')
			os.remove('Simh2_0001.ref')
			os.remove('Simh2_0001.maf')

			with open('Simh1.sam', 'w') as fh1:

				subprocess.call(['minimap2', '-ax', 'map-ont', new_mmi, 'Simh1_0001.fastq'], stdout=fh1)

			os.remove('Simh1_0001.fastq')

			with open(os.path.abspath('Simh2.sam'), 'w') as fh2:

				subprocess.call(['minimap2', '-ax', 'map-ont', new_mmi,'Simh2_0001.fastq'], stdout=fh2)

			os.remove(os.path.abspath('Simh2_0001.fastq'))

			with open(os.path.abspath('Simh1.bam'), 'w') as fh1:

				subprocess.call(['htsbox', 'samview', '-bS','Simh1.sam'], stdout=fh1)

			with open(os.path.abspath('Simh2.bam'), 'w') as fh2:

				subprocess.call(['htsbox', 'samview', '-bS','Simh2.sam'], stdout=fh2)

			with open(os.path.abspath('Simh1.srt.bam'), 'w') as fh1:

				subprocess.call(['samtools', 'sort', 'Simh1.bam'], stdout=fh1)

			with open(os.path.abspath('Simh2.srt.bam'), 'w') as fh2:

				subprocess.call(['samtools', 'sort', 'Simh2.bam'], stdout=fh2)

			os.remove('Simh1.sam')
			os.remove('Simh2.sam')
			os.remove('Simh1.bam')
			os.remove('Simh2.bam')
			subprocess.call(['samtools', 'index', 'Simh1.srt.bam']) #indexed haplotype1
			subprocess.call(['samtools', 'index', 'Simh2.srt.bam']) #indexed haplotype2

			
			#if any error occur when calling TRiCoLOR (.bams that are not readable for some reason) skip to next iteration (does not seem to be necessary)

			with open ('Sim.bed','w') as bedout:

				bedout.write(chromosome + '\t' + str(true_point-1000) + '\t' +  str(true_point+1000) + '\n')


			try:

				subprocess.check_call(['python','/home/bolognin/TRiCoLOR_py/TRiCoLOR.py','-g', reference, '-b','Sim.bed', '-b1', 'Simh1.srt.bam', '-b2', 'Simh2.srt.bam', '-m', str(kmer), '-t','4', '-o',os.path.abspath(out + '/Expansions/Sim'+str(i)+'/Res'+str(i)), '-mmi', new_mmi])
			
			except:

				continue


			#can we predict the exact number of tandem repeat contraction ?

			reference_repetitions=int(header1.split('_',4)[3].split(':')[1])
			repetition_type=header1.split('_',3)[2].split(':')[1]

			Hap1_Tab=pd.read_csv(os.path.abspath(out +'/Expansions/Sim'+str(i)+'/Res'+str(i)+'/haplotype1/RepetitionsTable.tsv'), sep='\t')
			Hap2_Tab=pd.read_csv(os.path.abspath(out +'/Expansions/Sim'+str(i)+'/Res'+str(i)+'/haplotype2/RepetitionsTable.tsv'), sep='\t')

			Hap1=[]
			Hap2=[]

			#look for the wanted repetition in range

			for l in range(len(Hap1_Tab['Start'])):

				if abs(Hap1_Tab['Start'][l] - true_point) <= acceptable or abs (Hap1_Tab['End'][l] -  (true_point+len(repetition_type)*reference_repetitions)) <= acceptable:

					Hap1.append((Hap1_Tab['Repetition Motif'][l],Hap1_Tab['Repetitions Number'][l]))


			for l in range(len(Hap2_Tab['Start'])):

				if abs(Hap2_Tab['Start'][l] - true_point) <= acceptable or abs (Hap2_Tab['End'][l] -  (true_point+len(repetition_type)*reference_repetitions)) <= acceptable:

					Hap2.append((Hap2_Tab['Repetition Motif'][l],Hap2_Tab['Repetitions Number'][l]))


			#sum the occurences, even if is only one: something left from the first time we simulate looking for general expansions and contractions

			Rep_Hap1=sum(n for _, n in Hap1)
			Rep_Hap2=sum(n for _, n in Hap2)


			Hap1_Res=[]
			Hap2_Res=[]
			Hap1_Res.append('Reference:'+str(reference_repetitions)+'-Haplo:'+str(Rep_Hap1)) 
			Hap2_Res.append('Reference:'+str(reference_repetitions)+'-Haplo:'+str(Rep_Hap2)) 
			Hap1_Type=['Reference:'+repetition_type+'-Haplo:'+Hap1[0][0]]
			Hap2_Type=['Reference:'+repetition_type+'-Haplo:'+Hap2[0][0]]
			TableRes=pd.DataFrame({'Hap1_Num':Hap1_Res,'Hap2_Num':Hap2_Res, 'Hap1_Core': Hap1_Type, 'Hap2_Core':Hap2_Type})
			TableRes.to_csv('Sim_'+str(i)+'_Results.tsv',sep='\t',index=False)


			os.chdir(os.path.abspath(out))


def Concat_Tables(list_of_paths):

	List_of_tables=[]

	for tab in list_of_paths:

		Tab=pd.read_csv(tab, sep='\t')
		List_of_tables.append(Tab)

	By_Row=pd.concat(List_of_tables, axis=0,ignore_index=True)

	return By_Row


def Shifted(word): #check if word is shifted by one

	if len(word) == 2:

		return word[1:] + word[0]

	else:

		possibilities=[]

		possibilities.append(word[1:] + word[0])
		possibilities.append(word[2] + word[:-1])

		return possibilities


def Precision_Recall_and_Exact_Core(number_of_simulations, out):


	size=5 #size of the simulated contraction/expansion 

	#Contractions

	TablesContractions=[]
	Simulations=number_of_simulations
	SimFolder=[os.path.abspath(out+'/Contractions/Sim' + str(i)) for i in range(1,Simulations+1)]


	for i in range(len(SimFolder)):

		Table_For_Sim=glob.glob(os.path.abspath(SimFolder[i]+'/*.tsv'))
		TablesContractions.extend(Table_For_Sim)

	TableDel=Concat_Tables(TablesContractions)



	#Exact Number of Repetitions +-1

	True_Positives1=0 #derived from hap1
	False_Negatives1=0

	True_Negatives1=0 #derived from hap2
	False_Positives1=0


	Correct_Core1=0
	Incorrect_Core1=0
	Transpose_Core1=0



	for i in range(len(TableDel['Hap1_Num'])):

		if abs(int(TableDel['Hap1_Num'][i].split(':')[1].split('-')[0])-int(TableDel['Hap1_Num'][i].split(':')[2])) <= size +1: #simulate a contraction of 5: see a contraction of  5 +- 1
			
			True_Positives1+=1

		else:

			False_Negatives1+=1


		if TableDel['Hap1_Core'][i].split(':')[1].split('-')[0] == TableDel['Hap1_Core'][i].split(':')[2]:


			Correct_Core1 +=1


		elif TableDel['Hap1_Core'][i].split(':')[2] in Shifted(TableDel['Hap1_Core'][i].split(':')[1].split('-')[0]):

			Transpose_Core1 += 1

		else:

			Incorrect_Core1 +=1

			
		if abs(int(TableDel['Hap2_Num'][i].split(':')[1].split('-')[0])-int(TableDel['Hap2_Num'][i].split(':')[2])) <= 1: #simulate nothing: see nothning +-1

			True_Negatives1+=1

		else:

			False_Positives1+=1

		if TableDel['Hap2_Core'][i].split(':')[1].split('-')[0] == TableDel['Hap2_Core'][i].split(':')[2]:

			Correct_Core1 +=1


		elif TableDel['Hap2_Core'][i].split(':')[2] in Shifted(TableDel['Hap2_Core'][i].split(':')[1].split('-')[0]):

			Transpose_Core1 += 1

		else:

			Incorrect_Core1 +=1



	DeletionRecall1=True_Positives1/(True_Positives1+False_Negatives1)
	DeletionPrecision1=True_Negatives1/(True_Negatives1+False_Positives1)
	Deletions_Precision_Recall1=(DeletionPrecision1,DeletionRecall1)


	#Exact_Core_Del=[[Correct_Core1,Transpose_Core1],[Incorrect_Core1,.0]] #fraction of correct cores. Others are shifted cores. For other deletions is the same, of course


	#Exact Number of Repetitions +-2


	True_Positives2=0 #derived from hap1
	False_Negatives2=0

	True_Negatives2=0 #derived from hap2
	False_Positives2=0



	for i in range(len(TableDel['Hap1_Num'])):

		if abs(int(TableDel['Hap1_Num'][i].split(':')[1].split('-')[0])-int(TableDel['Hap1_Num'][i].split(':')[2])) <= size +2: #simulate a contraction of 5: see a contraction of  5 +- 2
			
			True_Positives2+=1

		else:

			False_Negatives2+=1


		if abs(int(TableDel['Hap2_Num'][i].split(':')[1].split('-')[0])-int(TableDel['Hap2_Num'][i].split(':')[2])) <= 2: #simulate nothing: see nothning +-2

			True_Negatives2+=1

		else:

			False_Positives2+=1


	DeletionRecall2=True_Positives2/(True_Positives2+False_Negatives2)
	DeletionPrecision2=True_Negatives2/(True_Negatives2+False_Positives2)
	Deletions_Precision_Recall2=(DeletionPrecision2,DeletionRecall2)


	#Exact Number of Repetitions +-3


	True_Positives3=0 #derived from hap1
	False_Negatives3=0

	True_Negatives3=0 #derived from hap2
	False_Positives3=0



	for i in range(len(TableDel['Hap1_Num'])):

		if abs(int(TableDel['Hap1_Num'][i].split(':')[1].split('-')[0])-int(TableDel['Hap1_Num'][i].split(':')[2])) <= size +3: #simulate a contraction of 5: see a contraction of  5 +- 3
			
			True_Positives3+=1

		else:

			False_Negatives3+=1


		if abs(int(TableDel['Hap2_Num'][i].split(':')[1].split('-')[0])-int(TableDel['Hap2_Num'][i].split(':')[2])) <= 3: #simulate nothing: see nothning +-3

			True_Negatives3+=1

		else:

			False_Positives3+=1


	DeletionRecall3=True_Positives3/(True_Positives3+False_Negatives3)
	DeletionPrecision3=True_Negatives3/(True_Negatives3+False_Positives3)
	Deletions_Precision_Recall3=(DeletionPrecision3,DeletionRecall3)


	#Exact Number of Repetitions +-4


	True_Positives4=0 #derived from hap1
	False_Negatives4=0

	True_Negatives4=0 #derived from hap2
	False_Positives4=0



	for i in range(len(TableDel['Hap1_Num'])):

		if abs(int(TableDel['Hap1_Num'][i].split(':')[1].split('-')[0])-int(TableDel['Hap1_Num'][i].split(':')[2])) <= size +4: #simulate a contraction of 5: see a contraction of  5 +- 4
			
			True_Positives4+=1

		else:

			False_Negatives4+=1


		if abs(int(TableDel['Hap2_Num'][i].split(':')[1].split('-')[0])-int(TableDel['Hap2_Num'][i].split(':')[2])) <= 4: #simulate nothing: see nothning +-4

			True_Negatives4+=1

		else:

			False_Positives4+=1


	DeletionRecall4=True_Positives4/(True_Positives4+False_Negatives4)
	DeletionPrecision4=True_Negatives3/(True_Negatives4+False_Positives4)
	Deletions_Precision_Recall4=(DeletionPrecision4,DeletionRecall4)


	AllDeletions=[Deletions_Precision_Recall1,Deletions_Precision_Recall2,Deletions_Precision_Recall3,Deletions_Precision_Recall4]


	#Expansions


	TablesExpansions=[]


	Simulations=number_of_simulations
	SimFolder=[os.path.abspath(out+'/Expansions/Sim' + str(i)) for i in range(1,Simulations+1)]


	for i in range(len(SimFolder)):

		Table_For_Sim=glob.glob(os.path.abspath(SimFolder[i]+'/*.tsv'))
		TablesExpansions.extend(Table_For_Sim)

	TableExp=Concat_Tables(TablesExpansions)


	True_Positives1=0
	False_Negatives1=0

	True_Negatives1=0
	False_Positives1=0

	Correct_Core2=0
	Incorrect_Core2=0
	Transpose_Core2=0



	for i in range(len(TableExp['Hap1_Num'])):

		if abs(int(TableExp['Hap1_Num'][i].split(':')[1].split('-')[0])-int(TableExp['Hap1_Num'][i].split(':')[2])) <= size +1: #simulate a contraction of 5: see a contraction of  5 +- 1

			True_Positives1+=1

		else:

			False_Negatives1+=1

		if TableExp['Hap1_Core'][i].split(':')[1].split('-')[0] == TableExp['Hap1_Core'][i].split(':')[2]:


			Correct_Core2 +=1


		elif TableExp['Hap1_Core'][i].split(':')[2] in Shifted(TableExp['Hap1_Core'][i].split(':')[1].split('-')[0]):

			Transpose_Core2 +=1


		else:

			Incorrect_Core2 +=1


		if abs(int(TableExp['Hap2_Num'][i].split(':')[1].split('-')[0])-int(TableExp['Hap2_Num'][i].split(':')[2])) <= 1: #simulate nothing: see nothning +-1

			True_Negatives1+=1

		else:

			False_Positives1+=1

		if TableExp['Hap2_Core'][i].split(':')[1].split('-')[0] == TableExp['Hap2_Core'][i].split(':')[2]:


			Correct_Core2 +=1


		elif TableExp['Hap2_Core'][i].split(':')[2] in Shifted(TableExp['Hap2_Core'][i].split(':')[1].split('-')[0]):

			Transpose_Core2 +=1


		else:

			Incorrect_Core2+=1



	ExpansionRecall1=True_Positives1/(True_Positives1+False_Negatives1)
	ExpansionPrecision1=True_Negatives1/(True_Negatives1+False_Positives1)
	Expansions_Precision_Recall1=(ExpansionPrecision1,ExpansionRecall1)



	Exact_Core=[[(Correct_Core1+Correct_Core2)/(Correct_Core1+Correct_Core2+Transpose_Core1+Transpose_Core2+Incorrect_Core1+Incorrect_Core2),(Transpose_Core1+Transpose_Core2)/(Correct_Core1+Correct_Core2+Transpose_Core1+Transpose_Core2+Incorrect_Core1+Incorrect_Core2)],[(Incorrect_Core1+Incorrect_Core2)/(Correct_Core1+Correct_Core2+Transpose_Core1+Transpose_Core2+Incorrect_Core1+Incorrect_Core2),0]] #fraction of correct cores. Others are shifted cores. For other deletions is the same, of course



	True_Positives2=0
	False_Negatives2=0

	True_Negatives2=0
	False_Positives2=0


	for i in range(len(TableExp['Hap1_Num'])):

		if abs(int(TableExp['Hap1_Num'][i].split(':')[1].split('-')[0])-int(TableExp['Hap1_Num'][i].split(':')[2])) <= size +2: #simulate a contraction of 5: see a contraction of  5 +- 2

			True_Positives2+=1

		else:

			False_Negatives2+=1


		if abs(int(TableExp['Hap2_Num'][i].split(':')[1].split('-')[0])-int(TableExp['Hap2_Num'][i].split(':')[2])) <= 2: #simulate a contraction of 5: see a contraction of  5 +- 2
			
			True_Negatives2+=1

		else:

			False_Positives2+=1


	ExpansionRecall2=True_Positives2/(True_Positives2+False_Negatives2)
	ExpansionPrecision2=True_Negatives2/(True_Negatives2+False_Positives2)
	Expansions_Precision_Recall2=(ExpansionPrecision2,ExpansionRecall2)


	True_Positives3=0
	False_Negatives3=0

	True_Negatives3=0
	False_Positives3=0



	for i in range(len(TableExp['Hap1_Num'])):

		if abs(int(TableExp['Hap1_Num'][i].split(':')[1].split('-')[0])-int(TableExp['Hap1_Num'][i].split(':')[2])) <= size +3: #simulate a contraction of 5: see a contraction of  5 +- 3

			True_Positives3+=1

		else:

			False_Negatives3+=1


		if abs(int(TableExp['Hap2_Num'][i].split(':')[1].split('-')[0])-int(TableExp['Hap2_Num'][i].split(':')[2])) <= 3: #simulate a contraction of 5: see a contraction of  5 +- 3
			
			True_Negatives3+=1

		else:

			False_Positives3+=1


	ExpansionRecall3=True_Positives3/(True_Positives3+False_Negatives3)
	ExpansionPrecision3=True_Negatives3/(True_Negatives3+False_Positives3)
	Expansions_Precision_Recall3=(ExpansionPrecision3,ExpansionRecall3)


	True_Positives4=0
	False_Negatives4=0

	True_Negatives4=0
	False_Positives4=0



	for i in range(len(TableExp['Hap1_Num'])):

		if abs(int(TableExp['Hap1_Num'][i].split(':')[1].split('-')[0])-int(TableExp['Hap1_Num'][i].split(':')[2])) <= size +4: #simulate a contraction of 5: see a contraction of  5 +- 4

			True_Positives4+=1

		else:

			False_Negatives4+=1


		if abs(int(TableExp['Hap2_Num'][i].split(':')[1].split('-')[0])-int(TableExp['Hap2_Num'][i].split(':')[2])) <= 4: #simulate a contraction of 5: see a contraction of  5 +- 4
			
			True_Negatives4+=1

		else:

			False_Positives4+=1


	ExpansionRecall4=True_Positives4/(True_Positives4+False_Negatives4)
	ExpansionPrecision4=True_Negatives4/(True_Negatives4+False_Positives4)
	Expansions_Precision_Recall4=(ExpansionPrecision4,ExpansionRecall4)



	AllExpansions=[Expansions_Precision_Recall1,Expansions_Precision_Recall2,Expansions_Precision_Recall3,Expansions_Precision_Recall4]



	return [AllDeletions,AllExpansions],Exact_Core


def fmeasure(p, r):

	"""Calculates the f1 score for precision p and recall r."""

	return float(2*p*r / (p+r))


def fmeasureCurve(f, p):

	"""For a given f1 value and precision get the recall value."""

	return float(f * p / (2 * p - f))


def plotFMeasures(fstepsize=.1, stepsize=0.001):

	"""Plots 10 f1 score curves."""

	p = scipy.arange(0., 1., stepsize)[1:]

	for f in scipy.arange(0., 1., fstepsize)[1:]:

		points = [(x, fmeasureCurve(f, x)) for x in p if 0 < fmeasureCurve(f, x) <= 1.5]
		xs, ys = zip(*points)
		curve, = plt.plot(xs, ys, "--", color="lightgray", linewidth=0.5)
		plt.annotate(r"$f=%.1f$" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.05, ys[-10] - 0.035), size="small", color="lightgray")



def plotPrecisionRecallDiagram(title=None, points=None, labels=None,colors=None,lines=None, linecolors=None,loc=0):

	ax = plt.gca()   
	plt.title(title)
	plt.xlabel("Precision")
	plt.ylabel("Recall")
	plotFMeasures()
	scps = [] 

	for h in range(len(points)):

		for i, (x,y) in enumerate(points[h]):
			
			label = labels[h][i]
			color=colors[h][i]
			scp = ax.scatter(x, y, label=label, c=color,s=50, linewidths=0.75, alpha=0.75)
			scps.append(scp)

		if lines is not None:
			
			if linecolors is not None:

				lines = ax.plot([el[0] for el in points[h]],[el[1] for el in points[h]],":", color=linecolors[h], linewidth=0.5)

			else:

				lines= ax.plot([el[0] for el in points[h]],[el[1] for el in points[h]],":", color="black", linewidth=0.5)
				
		plt.legend(loc=loc, scatterpoints=1, numpoints=1, fancybox=False) #loc=0 guess the best location for legend
	
	plt.axis([-0.02, 1.02, -0.02, 1.02])


def ConstructPie(title=None, points=None):


	fig, ax = plt.subplots()
	size = 0.25
	vals = np.array(points)
	cmap = plt.get_cmap("tab20c")
	outer_colors = ['#0C5B29', '#A6A6A6']
	inner_colors = ['#479B42', '#18EA1C', '#A6A6A6','#A6A6A6']
	ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,wedgeprops=dict(width=size,edgecolor='w'))
	ax.pie(vals.flatten(), radius=1-size, colors=inner_colors,wedgeprops=dict(width=size,edgecolor='w'))
	ax.set(aspect="equal", title=title)
	colors=['#0C5B29', '#479B42', '#18EA1C', '#A6A6A6']
	pre=['Correct Motif: ', 'Incorrect Motif: ', 'Perfect Motif: ', 'Slipped Motif: ', ]

	labels=[]

	for i in vals.tolist():

		labels.append(str(sum(i)))

	for i in vals.tolist():

		for el in i:

			labels.append(str(el))

	labels=labels[:-2]

	finals=[pre[i] + labels[i][:5] for i in range(len(labels))]

	plt.legend(labels=finals, loc='center')




#Executing codes...


for simtype in ['contraction', 'expansion']:

	Simulate(bed, reference, number_of_simulations, sim_type=simtype, model_qc, accuracy, out)

	#or run in parallel in 2 different ipython idles
	#Simulate('/home/bolognin/DataForSim/KnownRepetitions.bed', '/home/bolognin/TRiCoLOR/GenomeDir/GRCh38Decoy.fa', 100, 'contraction', '/home/bolognin/DataForSim/model_qc_clr', 0.9, '/home/bolognin/Simulations_Accuracy_09')
	#Simulate('/home/bolognin/DataForSim/KnownRepetitions.bed', '/home/bolognin/TRiCoLOR/GenomeDir/GRCh38Decoy.fa', 100, 'expansion', '/home/bolognin/DataForSim/model_qc_clr', 0.9, '/home/bolognin/Simulations_Accuracy_09')

pr,core=Precision_Recall_and_Exact_Core(number_of_simulations, out)
colors=[['#8FC3C9','#309FE7','#2B35AB','#26263B'],['#D38D8D','#EB7F7F','#DC4444', '#E21515']]
labels=[['EC+-1','EC+-2','EC+-3','EC+-4'],['EE+-1','EE+-2','EE+-3','EE+-4']]


plotPrecisionRecallDiagram('Precision and Recall, Accuracy = .9',pr, labels=labels, colors=colors)
plt.show()

ConstructPie('Motif Matching, Accuracy = .9', core)
plt.show()
