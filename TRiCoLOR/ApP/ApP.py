#!/usr/bin/python env

#Python 3 standard library

import os
import glob
import math
import logging
from collections import defaultdict
from shutil import which
import sys
import subprocess
import csv

# additional modules

import pyfaidx
import pysam
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot


def run(parser, args):


	if not os.path.exists(os.path.abspath(args.output)): #folder does not exist

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('Cannot create the output folder') #no write permission, probably			
			sys.exit(1)

	else: #path already exists

		if not os.access(os.path.abspath(args.output),os.W_OK): #folder exists but no write permissions

			print('Missing write permissions on the output folder')			
			sys.exit(1)
			
		elif os.listdir(os.path.abspath(args.output)): #folder exists but isn't empty. Do not overwrite results.

			print('The output folder is not empty. Specify another output folder or clean the previsouly chosen')
			sys.exit(1)

	command_dict= vars(args)
	
	notkey=['func']
	command_string= ' '.join("{}={}".format(key,val) for key,val in command_dict.items() if key not in notkey)

	logging.basicConfig(filename=os.path.abspath(args.output + '/TRiCoLOR_ApP.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

	logging.info('main=TRiCoLOR ' + command_string)

	if which('samtools') is None:

		logging.error('samtools cannot be executed. Install samtools and re-run TRiCoLOR ApP')
		sys.exit(1)

	#check if the genome file exists, is readable and is in .fasta format

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>')) #genome .file starts with '>'

	except:

		logging.error('Reference file does not exist, is not readable or is not a valid FASTA')
		sys.exit(1)

	bams=args.bamfile[0]

	if len(bams) > 2:

		logging.error('TRiCoLOR supports haploid and diploid genomes only')

	for bam in bams:

		try:

			subprocess.check_call(['samtools','quickcheck', os.path.abspath(bam)],stderr=open(os.devnull, 'wb'))

		except:

			logging.error('BAM ' + bam + ' does not exist, is not readable or is not a valid BAM')
			sys.exit(1)


		if not os.path.exists(os.path.abspath(bam + '.bai')):

			logging.error('Missing index for BAM ' + bam + '. Not a BAM generated with REFER')
			sys.exit(1)

	logging.info('Ploidy: ' + str(len(bams)))

	b_in=Bed_Reader(args.bedfile)
	it_ = iter(b_in)
	labels_set=set()

	logging.info('Analyzing ...')
	logging.info('Regions in BED: ' + str(b_in.length()))

	skippedlabel=0
	skippederror=0

	for i in range(b_in.length()):

		chromosome, start, end,label = next(it_)

		if label in labels_set:

			logging.warning('Skipped region with duplicated label')
			skippedlabel+=1
			
		else:

			try:

				Generate_Alignment_ToPlot(args.genome,bams,chromosome,start,end,args.genomebed,args.haplotypebed,label,args.output)

			except:

				logging.exception('Unexpected error while analyzing region ' + chromosome + ':' + str(start) + '-' +str(end) + ". Log is below")
				skippederror+=1

	logging.info('Regions with duplicated labels skipped: ' + str(skippedlabel))
	logging.info('Regions with unexpected errors: ' + str(skippederror))
	logging.info('Done')


class Bed_Reader():


	def __init__(self,bedfile):

		self.bedfile=bedfile

	def __iter__(self):

		with open (self.bedfile, 'r') as bedin:

			for line in csv.reader(bedin, delimiter='\t'):

				if not line[0].startswith('#') and line !=[]: 

					if len(line) < 4:

						logging.error('BED to TRiCoLOR ApP -bed/--bedfile must contain at least: chromosome, start, end, label (other fields are ignored)')
						sys.exit(1)

					else:

						yield (line[0], int(line[1]), int(line[2]), line[3]) # exclude other fields if they are present

	def length(self):

		with open(self.bedfile, 'r') as bedin:

			size=sum(1 for _ in bedin if not _.startswith('#') and not _.strip()=='')

			return size


def Get_Alignment_Positions(bamfile,chromosome,start,end):


	coords=[]
	seq=[]

	BamFile=pysam.AlignmentFile(bamfile,'rb')

	if start==end:

		end=start+1


	for read in BamFile.fetch(chromosome, start-1, end-1): #fetch on zero based positions

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			coords = read.get_reference_positions(full_length=True)
			seq=read.seq

	return coords,seq


def list_duplicates(list_of_seq):


	a_=defaultdict(list)

	for i,item in enumerate(list_of_seq):

		a_[item].append(i)

	return ((key,locs) for key,locs in a_.items() if len(locs) > 1)


def modifier(coordinates): #fast way to remove None and substitute with closest number in list

    
    coordinates=[el+1 if el is not None else el for el in coordinates] #get true coordinates
    start=next(ele for ele in coordinates if ele is not None)


    for ind, ele in enumerate(coordinates):
        
        if ele is None:
                
            coordinates[ind] = start
        
        else:

            start = ele

    return coordinates


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

			coords_purified.extend([None]*(int(coords_without_insertions[i+1])-int(coords_without_insertions[i]-1)))

			NewSeq+=seq[i]
			NewSeq+='-'*(int(coords_without_insertions[i+1])-int(coords_without_insertions[i]-1))

		else:

			coords_purified.append(coords_without_insertions[i])
			NewSeq+=seq[i]


	coords_purified.append(coords_without_insertions[-1])
	NewSeq += seq[-1]

	return coords_purified,NewSeq


def Generate_Alignment_ToPlot(reference_fasta,hap_bam,chromosome,start,end,ref_table,hap_table,label,out):


	if len(hap_bam) == 1:

		bam1_coords,bam1_seq=Get_Alignment_Positions(hap_bam[0],chromosome,start,end)
		bam2_coords,bam2_seq=[],[]


	else:

		bam1_coords,bam1_seq=Get_Alignment_Positions(hap_bam[0],chromosome,start,end)
		bam2_coords,bam2_seq=Get_Alignment_Positions(hap_bam[1],chromosome,start,end)



	if hap_table is None:

		hap1_table = None
		hap2_table = None


	else:

		if len(hap_table[0]) == 1:

			if len(hap_bam) == 1:

				hap1_table = hap_table[0][0]
				hap2_table = None

			else: #len is 2 and cannot be sure which one the BED refers to

				logging.warning('Double input to -bam/--bamfile and single input to -hb/--haplotypebed: cannot say to which BAM BED refers. Skipped repetitions highlighting.')

				hap1_table = None
				hap2_table = None

		else:

			hap1_table = hap_table[0][0]
			hap2_table = hap_table[0][1]



	if len(bam1_seq) == 0 and len(bam2_seq) == 0:


		min_=start-1000
		max_=end+1000


		ref=pyfaidx.Fasta(reference_fasta)
		chrom=ref[chromosome]
		ref_seq=chrom[:len(chrom)].seq
		ref_seq=ref_seq[min_-1:max_]

		Reference_Trace = go.Scatter(
			x = list(range(min_,max_+1)), # works with positions
			y = [0]*len(ref_seq),
			name = 'Reference',
			text = list(ref_seq),
			hoverinfo = 'text+name+x',
			yaxis = 'y',
			mode='markers+lines',
			line=dict(dash='dash',color='rgba(116, 116, 116, 1)'),
			marker=dict(color='rgba(116, 116, 116, 1)')
		)


		data = [Reference_Trace]


		if ref_table is None:

			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'))
			fig = dict(data=data, layout=layout)
			plot(fig,filename=os.path.abspath(out)+'/'+label+'.html',auto_open=False)


		else:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)


			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		
		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(out)+'/'+label+'.html',auto_open=False)


	elif len(bam1_seq) == 0 and len(bam2_seq) != 0:


		bam2_mod_coords,bam2_mod_seq=Modifier(bam2_coords,bam2_seq)

		min_=math.floor(next(item for item in bam2_mod_coords if item is not None))
		max_=math.ceil(next(item for item in reversed(bam2_mod_coords) if item is not None))



		ref=pyfaidx.Fasta(reference_fasta)
		chrom=ref[chromosome]
		ref_seq=chrom[:len(chrom)].seq
		ref_seq=ref_seq[min_-1:max_]

		Reference_Trace = go.Scatter(
			x = list(range(min_,max_+1)), # works with positions
			y = [0]*len(ref_seq),
			name = 'Reference',
			text = list(ref_seq),
			hoverinfo = 'text+name+x',
			yaxis = 'y',
			mode='markers+lines',
			line=dict(dash='dash',color='rgba(116, 116, 116, 1)'),
			marker=dict(color='rgba(116, 116, 116, 1)')
		)


		Haplo2_Trace_Lower = go.Scatter(
			x = bam2_mod_coords, # works with positions
			y = [-1]*len(bam2_mod_seq),
			name = 'Haplotype 2',
			text = list(bam2_mod_seq),
			hoverinfo = 'text+name+x',
			yaxis = 'y',
			mode='markers+lines',
			line=dict(dash='dash',color='rgb(152, 77, 77)'),
			marker=dict(color='rgb(152, 77, 77)')
		)



		data = [Reference_Trace,Haplo2_Trace_Lower]

		if ref_table is None and hap2_table is None:

			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'))
			fig = dict(data=data, layout=layout)
			plot(fig,filename=os.path.abspath(out)+'/'+label+'.html',auto_open=False)


		elif ref_table is not None and hap2_table is None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)


			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)



		elif ref_table is None and hap2_table is not None:


			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]-.1,y0=-.15,x1=hap2_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 2 Repetitions',method = 'relayout',args = ['shapes', cluster_hap2])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)



		elif ref_table is not None and hap2_table is not None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)


			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]-.1,y0=-1.15,x1=hap2_rep['End'][i]+.1,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap2])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)


		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(out)+'/'+label+'.html',auto_open=False)



	elif len(bam1_seq) != 0 and len(bam2_seq) == 0:

		bam1_mod_coords,bam1_mod_seq=Modifier(bam1_coords,bam1_seq)

		min_=math.floor(next(item for item in bam1_mod_coords if item is not None))
		max_=math.ceil(next(item for item in reversed(bam1_mod_coords) if item is not None))


		ref=pyfaidx.Fasta(reference_fasta)
		chrom=ref[chromosome]
		ref_seq=chrom[:len(chrom)].seq
		ref_seq=ref_seq[min_-1:max_]


		Reference_Trace = go.Scatter(
			x = list(range(min_,max_+1)), # works with positions
			y = [0]*len(ref_seq),
			name = 'Reference',
			text = list(ref_seq),
			hoverinfo = 'text+name+x',
			yaxis = 'y',
			mode='markers+lines',
			line=dict(dash='dash',color='rgba(116, 116, 116, 1)'),
			marker=dict(color='rgba(116, 116, 116, 1)')
		)


		Haplo1_Trace_Upper = go.Scatter(
			x = bam1_mod_coords, # works with positions
			y = [1]*len(bam1_mod_seq),
			name = 'Haplotype 1',
			text = list(bam1_mod_seq),
			hoverinfo = 'text+name+x',
			yaxis = 'y',
			mode='markers+lines',
			line=dict(dash='dash',color='rgb(66, 96, 154)'),
			marker=dict(color='rgb(66, 96, 154)')
		)


		data = [Reference_Trace,Haplo1_Trace_Upper]

		if ref_table is None and hap1_table is None:

			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'))
			fig = dict(data=data, layout=layout)
			plot(fig,filename=os.path.abspath(out)+'/'+label+'.html',auto_open=False)


		elif ref_table is not None and hap1_table is None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)


			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)



		elif ref_table is None and hap1_table is not None:


			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]-.1,y0=.85,x1=hap1_rep['End'][i]+.1,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 1 Repetitions',method = 'relayout',args = ['shapes', cluster_hap1])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)




		elif ref_table is not None and hap1_table is not None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)


			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]-.1,y0=.85,x1=hap1_rep['End'][i]+.1,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)


			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 1 Repetitions',method ='relayout',args = ['shapes', cluster_hap1]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap1])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)


		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(out)+'/'+label+'.html',auto_open=False)


	else:


		bam1_mod_coords,bam1_mod_seq=Modifier(bam1_coords,bam1_seq)
		bam2_mod_coords,bam2_mod_seq=Modifier(bam2_coords,bam2_seq)


		#get minimum and maximum for the bams coordinates, so that the reference sequence is rescaled to that range
		mins=[((next(item for item in bam1_mod_coords if item is not None)),(next(item for item in bam2_mod_coords if item is not None)))][0]
		maxs=[((next(item for item in reversed(bam1_mod_coords) if item is not None)),(next(item for item in reversed(bam2_mod_coords) if item is not None)))][0]

		min_=math.floor(min(mins))
		max_=math.ceil(max(maxs))

		ref=pyfaidx.Fasta(reference_fasta)
		chrom=ref[chromosome]
		ref_seq=chrom[:len(chrom)].seq
		ref_seq=ref_seq[min_-1:max_]

		Reference_Trace = go.Scatter(
			x = list(range(min_,max_+1)), # works with positions
			y = [0]*len(ref_seq),
			name = 'Reference',
			text = list(ref_seq),
			hoverinfo = 'text+name+x',
			yaxis = 'y',
			mode='markers+lines',
			line=dict(dash='dash',color='rgba(116, 116, 116, 1)'),
			marker=dict(color='rgba(116, 116, 116, 1)')
		)

		Haplo1_Trace_Upper = go.Scatter(
			x = bam1_mod_coords, # works with positions
			y = [1]*len(bam1_mod_seq),
			name = 'Haplotype 1',
			text = list(bam1_mod_seq),
			hoverinfo = 'text+name+x',
			yaxis = 'y',
			mode='markers+lines',
			line=dict(dash='dash',color='rgb(66, 96, 154)'),
			marker=dict(color='rgb(66, 96, 154)')
		)

		Haplo2_Trace_Lower = go.Scatter(
			x = bam2_mod_coords, # works with positions
			y = [-1]*len(bam2_mod_seq),
			name = 'Haplotype 2',
			text = list(bam2_mod_seq),
			hoverinfo = 'text+name+x',
			yaxis = 'y',
			mode='markers+lines',
			line=dict(dash='dash',color='rgb(152, 77, 77)'),
			marker=dict(color='rgb(152, 77, 77)')
		)


		data = [Reference_Trace,Haplo1_Trace_Upper,Haplo2_Trace_Lower]

		if ref_table is None and hap1_table is None and hap2_table is None:

			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'))
			fig = dict(data=data, layout=layout)
			plot(fig,filename=os.path.abspath(out)+'/'+label+'.html',auto_open=False)


		elif ref_table is not None and hap1_table is None and hap2_table is None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)


			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)


		elif ref_table is None and hap1_table is not None and hap2_table is None:


			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]-.1,y0=.85,x1=hap1_rep['End'][i]+.1,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 1 Repetitions',method = 'relayout',args = ['shapes', cluster_hap1])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)


		elif ref_table is None and hap1_table is None and hap2_table is not None:


			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]-.1,y0=-1.15,x1=hap2_rep['End'][i]+.1,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 2 Repetitions',method = 'relayout',args = ['shapes', cluster_hap2])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)


		elif ref_table is not None and hap1_table is not None and hap2_table is None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)


			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]-.1,y0=.85,x1=hap1_rep['End'][i]+.1,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)


			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 1 Repetitions',method ='relayout',args = ['shapes', cluster_hap1]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap1])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)


		elif ref_table is not None and hap1_table is None and hap2_table is not None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)


			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]-.1,y0=-1.15,x1=hap2_rep['End'][i]+.1,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap2])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)


		elif ref_table is None and hap1_table is not None and hap2_table is not None:

			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]-.1,y0=.85,x1=hap1_rep['End'][i]+.1,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)


			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]-.1,y0=-1.15,x1=hap2_rep['End'][i]+.1,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 1 Repetitions',method = 'relayout',args = ['shapes', cluster_hap1]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_hap1+cluster_hap2])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		else:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]-.1,y0=-.15,x1=ref_rep['End'][i]+.1,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]-.1,y0=.85,x1=hap1_rep['End'][i]+.1,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)


			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]-.1,y0=-1.15,x1=hap2_rep['End'][i]+.1,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 1 Repetitions',method ='relayout',args = ['shapes', cluster_hap1]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'All',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap1+cluster_hap2])]))])
			layout = dict(title='Repetitions in '+chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)


		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(out)+'/'+label+'.html',auto_open=False)


