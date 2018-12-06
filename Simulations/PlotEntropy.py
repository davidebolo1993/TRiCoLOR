#!/usr/bin/python


#Plot entropy results for simulations. Useful to identify a treshold (mean - 3 stddevs, the one that we used) that allows to discriminate between sequences with a repetition (lower entropy) and sequences without repetitions (higher entropy). The scanning size is 20 bp. The resulting treshold was 1.1 for consensus .bams and 1.3 for raw .bams


import os
import glob
import math
import random
import plotly.graph_objs as go
from plotly.offline import plot
import plotly.figure_factory as ff


def mean(data):

	n = len(data)

	return sum(data)/n 

def _ss(data):

	c = mean(data)
	ss = sum((x-c)**2 for x in data)
	
	return ss

def stddev(data, ddof=0):

	n = len(data)
	ss = _ss(data)
	pvar = ss/(n-ddof)

	return pvar**0.5


def entropy(string):

	prob = [float(string.count(c)) / len(string) for c in dict.fromkeys(list(string))]
	entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])

	return entropy


def BamScanner(bamfile,scansize=20):

	BamFile=pysam.AlignmentFile(bamfile,"rb")

	ent__=[]

	for read in BamFile.fetch():

		if not read.is_unmapped and not read.is_secondary:

			sequence=read.seq

			start=0
			end=scansize

			while len(sequence)-1 >= start:


				if len(sequence)-start >= scansize:

					ent__.append(entropy(sequence[start:end]))
					start+=scansize
					end+=scansize

				else:

					ent__.append(entropy(sequence[start:len(sequence)]))
					break

	return ent__



def Get_Consensus_Entropy_From_Simulations(path, number_of_simulations):


	Insertions_bam=[glob.glob(os.path.abspath(path + '/Insertions/Sim'+str(i)+'/Res'+str(i)+'/haplotype1')+'/*.srt.bam') for i in range(1,number_of_simulations+1)]
	Insertions_bam=[x for x in Insertions_bam if x != []]


	Normals_bam=[glob.glob(os.path.abspath(path+ '/Insertions/Sim'+str(i)+'/Res'+str(i)+'/haplotype2')+'/*.srt.bam') for i in range(1,number_of_simulations+1)]
	Normals_bam=[x for x in Normals_bam if x != []]


	Deletions_bam=[glob.glob(os.path.abspath(path+'/Deletions/Sim'+str(i)+'/Res'+str(i)+'/haplotype1')+'/*.srt.bam') for i in range(1,number_of_simulations+1)]
	Deletions_bam=[x for x in Deletions_bam if x != []]

	
	#Insertions

	In_en=[]

	for bam in Insertions_bam:

		en_=BamScanner(bam[0])		
		
    	if en_ != []:
			
		In_en.append(en_)

	In_distributions=[]

	for el in In_en:

		In_distributions.extend(el[0:-1])

	In_MS=mean(In_distributions) - 3 *stddev(In_distributions)


	#Normals

	No_en=[]

	for bam in Normals_bam:

		en_=BamScanner(bam[0])
		
		if en_ != []:
			
			No_en.append(en_)

	No_distributions=[]

	for el in No_en:

		No_distributions.extend(el[0:-1])

	No_MS=mean(No_distributions) - 3 *stddev(No_distributions)


	#Deletions

	Del_en=[]

	for bam in Deletions_bam:

		en_=BamScanner(bam[0])
    
		if en_ != []:
			
			Del_en.append(en_)

	Del_distributions=[]

	for el in Del_en:

		Del_distributions.extend(el[0:-1])

	Del_MS=mean(Del_distributions) - 3 *stddev(Del_distributions)



	colors = ['#333F44', '#37AA9C', '#94F3E4']
	fig = ff.create_distplot([In_distributions,No_distributions,Del_distributions], ["I: "+str(In_MS), "N: "+str(No_MS),"D: "+str(Del_MS)], show_hist=False, colors=colors)


	plot(fig, filename=os.path.abspath(path + '/ConsensusEntropy.html'), auto_open=False)


#When we are scanning the genome for the first time, we use raw reads, so the treshold must be calculated using raw reads. The previous was just an example. 
	
def Get_Raw_Entropy_From_Simulations(path, number_of_simulations):


	Insertions_bam=[os.path.abspath(path + '/Insertions/Sim'+str(i)+'/InsSimh1_'+str(i)+'.srt.bam') for i in range(1,number_of_simulations+1)][:10]
	Insertions_bam=[x for x in Insertions_bam if x != []]


	Normals_bam=[os.path.abspath(path + '/Insertions/Sim'+str(i)+'/InsSimh2_'+str(i)+'.srt.bam') for i in range(1,number_of_simulations+1)][:10]
	Normals_bam=[x for x in Normals_bam if x != []]


	Deletions_bam=[os.path.abspath(path + '/Deletions/Sim'+str(i)+'/DelSimh1_'+str(i)+'.srt.bam') for i in range(1,number_of_simulations+1)][:10]
	Deletions_bam=[x for x in Deletions_bam if x != []]

	
	#Insertions

	In_en=[]

	for bam in Insertions_bam:

		en_=BamScanner(bam)	
		In_en.append(en_)

	In_distributions=[]

	for el in In_en:

		In_distributions.extend(el[0:-1])

	In_MS=mean(In_distributions) - 3 *stddev(In_distributions)


	#Normal

	No_en=[]

	for bam in Normals_bam:

		en_=BamScanner(bam)
		No_en.append(en_)

	No_distributions=[]

	for el in No_en:

		No_distributions.extend(el[0:-1])


	No_MS=mean(No_distributions) - 3 *stddev(No_distributions)


	#Deletions

	Del_en=[]

	for bam in Deletions_bam:


		en_=BamScanner(bam)	
		Del_en.append(en_)

	Del_distributions=[]

	for el in Del_en:

		Del_distributions.extend(el[0:-1])

	Del_MS=mean(Del_distributions) - 3 *stddev(Del_distributions)


	#As there are too many points, they must be rescaled before plotting. 10000 for each is fine.

	In_distributions=random.sample(In_distributions,10000)
	No_distributions=random.sample(No_distributions,10000)
	Del_distributions=random.sample(Del_distributions,10000)



	colors = ['#333F44', '#37AA9C', '#94F3E4']
	fig = ff.create_distplot([In_distributions,No_distributions,Del_distributions], ["I: "+str(In_MS), "N: "+str(No_MS),"D: "+str(Del_MS)], show_hist=False, colors=colors)


	plot(fig, filename=os.path.abspath(path + '/RawEntropy.html'), auto_open=False)

