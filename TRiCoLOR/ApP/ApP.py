#!/usr/bin/python3 env

#Python 3 standard library

import os
import glob
import math
import sys
import re
from collections import defaultdict
from datetime import datetime

# additional modules

import pyfaidx
import pysam
import pandas as pd 
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot


class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	BAM = list()
	OUT = ''
	REF = ''
	chromosome = ''
	start = 0
	end = 0
	label=''
	genomebed=None
	haplotype1bed=None
	haplotype2bed=None


def Get_Alignment_Positions(bam,c):

	'''
	Get coordinates and sequence of interest
	'''

	seq,coords=[],[]

	bamfile=pysam.AlignmentFile(bam,'rb')

	for read in bamfile.fetch(c.chromosome,c.start,c.end):

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary: #we kept only primary consensus alignment

			coords = read.get_reference_positions(full_length=True)
			seq= read.seq

	return coords,seq


def list_duplicates(list_of_seq):

	'''
	Check duplicated entries (insertions) in list of coordinates
	'''

	a_=defaultdict(list)

	for i,item in enumerate(list_of_seq):

		a_[item].append(i)

	return ((key,locs) for key,locs in a_.items() if len(locs) > 1)


def modifier(coordinates):

	'''
	Substitute None with closest not None
	'''

	coordinates=[el+1 if el is not None else el for el in coordinates]
	start=next(ele for ele in coordinates if ele is not None)

	for ind, ele in enumerate(coordinates):
		
		if ele is None:
				
			coordinates[ind] = start
		
		else:

			start = ele

	return coordinates


def Modifier(coords,seq):

	'''
	Modify coordinates. This is a trick-function that modify sequence and coordinates so that they can be plotted properly with respect to reference
	'''

	coords_without_insertions=modifier(coords)

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


def Generate_Alignment_ToPlot(c):

	'''
	Generate plot and store
	'''

	hap1_table=c.haplotype1bed
	hap2_table=c.haplotype2bed
	ref_table=c.genomebed

	bam1_coords,bam1_seq=Get_Alignment_Positions(c.BAM[0],c)
	bam2_coords,bam2_seq=Get_Alignment_Positions(c.BAM[1],c)


	if len(bam1_seq) == 0 and len(bam2_seq) == 0:

		min_=c.start-1000
		max_=c.end+1000

		ref=pyfaidx.Fasta(c.REF)
		chrom=ref[c.chromosome]
		ref_seq=chrom[:len(chrom)].seq
		ref_seq=ref_seq[min_-1:max_]

		Reference_Trace = go.Scatter(
			x = list(range(min_,max_+1)),
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

			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'))
			fig = dict(data=data, layout=layout)
			plot(fig,filename=os.path.abspath(c.OUT+'/'+c.label+'.html'),auto_open=False)

		else:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)
		
		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(c.OUT+'/'+c.label+'.html'),auto_open=False)

	elif len(bam1_seq) == 0 and len(bam2_seq) != 0:

		bam2_mod_coords,bam2_mod_seq=Modifier(bam2_coords,bam2_seq)

		min_=math.floor(next(item for item in bam2_mod_coords if item is not None))
		max_=math.ceil(next(item for item in reversed(bam2_mod_coords) if item is not None))

		ref=pyfaidx.Fasta(c.REF)
		chrom=ref[c.chromosome]
		ref_seq=chrom[:len(chrom)].seq
		ref_seq=ref_seq[min_-1:max_]

		Reference_Trace = go.Scatter(
			x = list(range(min_,max_+1)),
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
			x = bam2_mod_coords,
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

			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'))
			fig = dict(data=data, layout=layout)
			plot(fig,filename=os.path.abspath(c.OUT+'/'+c.label+'.html'),auto_open=False)

		elif ref_table is not None and hap2_table is None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is None and hap2_table is not None:

			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			hap2_rep=hap2_rep.loc[hap2_rep['#Chromosome'] == c.chromosome]
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]+.2,y0=-.15,x1=hap2_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 2 Repetitions',method = 'relayout',args = ['shapes', cluster_hap2])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is not None and hap2_table is not None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			hap2_rep=hap2_rep.loc[hap2_rep['#Chromosome'] == c.chromosome]
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]+.2,y0=-1.15,x1=hap2_rep['End'][i]+.8,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap2])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(c.OUT+'/'+c.label+'.html'),auto_open=False)

	elif len(bam1_seq) != 0 and len(bam2_seq) == 0:

		bam1_mod_coords,bam1_mod_seq=Modifier(bam1_coords,bam1_seq)

		min_=math.floor(next(item for item in bam1_mod_coords if item is not None))
		max_=math.ceil(next(item for item in reversed(bam1_mod_coords) if item is not None))

		ref=pyfaidx.Fasta(c.REF)
		chrom=ref[c.chromosome]
		ref_seq=chrom[:len(chrom)].seq
		ref_seq=ref_seq[min_-1:max_]

		Reference_Trace = go.Scatter(
			x = list(range(min_,max_+1)),
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
			x = bam1_mod_coords,
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

			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'))
			fig = dict(data=data, layout=layout)
			plot(fig,filename=os.path.abspath(c.OUT+'/'+c.label+'.html'),auto_open=False)

		elif ref_table is not None and hap1_table is None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is None and hap1_table is not None:

			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			hap1_rep=hap1_rep.loc[hap1_rep['#Chromosome'] == c.chromosome]
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]+.2,y0=.85,x1=hap1_rep['End'][i]+.8,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 1 Repetitions',method = 'relayout',args = ['shapes', cluster_hap1])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is not None and hap1_table is not None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			hap1_rep=hap1_rep.loc[hap1_rep['#Chromosome'] == c.chromosome]
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]+.2,y0=.85,x1=hap1_rep['End'][i]+.8,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 1 Repetitions',method ='relayout',args = ['shapes', cluster_hap1]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap1])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(c.OUT+'/'+c.label+'.html'),auto_open=False)

	else:

		bam1_mod_coords,bam1_mod_seq=Modifier(bam1_coords,bam1_seq)
		bam2_mod_coords,bam2_mod_seq=Modifier(bam2_coords,bam2_seq)
		mins=[((next(item for item in bam1_mod_coords if item is not None)),(next(item for item in bam2_mod_coords if item is not None)))][0]
		maxs=[((next(item for item in reversed(bam1_mod_coords) if item is not None)),(next(item for item in reversed(bam2_mod_coords) if item is not None)))][0]
		min_=math.floor(min(mins))
		max_=math.ceil(max(maxs))
		ref=pyfaidx.Fasta(c.REF)
		chrom=ref[c.chromosome]
		ref_seq=chrom[:len(chrom)].seq
		ref_seq=ref_seq[min_-1:max_]

		Reference_Trace = go.Scatter(
			x = list(range(min_,max_+1)),
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
			x = bam1_mod_coords,
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
			x = bam2_mod_coords,
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

			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'))
			fig = dict(data=data, layout=layout)
			plot(fig,filename=os.path.abspath(c.OUT+'/'+c.label+'.html'),auto_open=False)

		elif ref_table is not None and hap1_table is None and hap2_table is None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is None and hap1_table is not None and hap2_table is None:

			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			hap1_rep=hap1_rep.loc[hap1_rep['#Chromosome'] == c.chromosome]
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]+.2,y0=.85,x1=hap1_rep['End'][i]+.8,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 1 Repetitions',method = 'relayout',args = ['shapes', cluster_hap1])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is None and hap1_table is None and hap2_table is not None:

			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			hap2_rep=hap2_rep.loc[hap2_rep['#Chromosome'] == c.chromosome]
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]+.2,y0=-1.15,x1=hap2_rep['End'][i]+.8,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 2 Repetitions',method = 'relayout',args = ['shapes', cluster_hap2])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is not None and hap1_table is not None and hap2_table is None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			hap1_rep=hap1_rep.loc[hap1_rep['#Chromosome'] == c.chromosome]
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]+.2,y0=.85,x1=hap1_rep['End'][i]+.8,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 1 Repetitions',method ='relayout',args = ['shapes', cluster_hap1]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap1])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is not None and hap1_table is None and hap2_table is not None:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			hap2_rep=hap2_rep.loc[hap2_rep['#Chromosome'] == c.chromosome]
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]+.2,y0=-1.15,x1=hap2_rep['End'][i]+.8,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap2])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		elif ref_table is None and hap1_table is not None and hap2_table is not None:

			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			hap1_rep=hap1_rep.loc[hap1_rep['#Chromosome'] == c.chromosome]
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]+.2,y0=.85,x1=hap1_rep['End'][i]+.8,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)

			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			hap2_rep=hap2_rep.loc[hap2_rep['#Chromosome'] == c.chromosome]
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]+.2,y0=-1.15,x1=hap2_rep['End'][i]+.8,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 1 Repetitions',method = 'relayout',args = ['shapes', cluster_hap1]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_hap1+cluster_hap2])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		else:

			ref_rep=pd.read_csv(os.path.abspath(ref_table),sep='\t')
			ref_rep=ref_rep.loc[ref_rep['#Chromosome'] == c.chromosome]
			cluster_ref=[]

			for i in range(len(ref_rep['Start'])):

				if ref_rep['Start'][i] >= min_ and ref_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep['Start'][i]+.2,y0=-.15,x1=ref_rep['End'][i]+.8,y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
					cluster_ref.append(cluster_dict)

			hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep='\t')
			hap1_rep=hap1_rep.loc[hap1_rep['#Chromosome'] == c.chromosome]
			cluster_hap1=[]

			for i in range(len(hap1_rep['Start'])):

				if hap1_rep['Start'][i] >= min_ and hap1_rep['End'][i] <= max_: 

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep['Start'][i]+.2,y0=.85,x1=hap1_rep['End'][i]+.8,y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
					cluster_hap1.append(cluster_dict)

			hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep='\t')
			hap2_rep=hap2_rep.loc[hap2_rep['#Chromosome'] == c.chromosome]
			cluster_hap2=[]

			for i in range(len(hap2_rep['Start'])):

				if hap2_rep['Start'][i] >= min_ and hap2_rep['End'][i] <= max_:

					cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep['Start'][i]+.2,y0=-1.15,x1=hap2_rep['End'][i]+.8,y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
					cluster_hap2.append(cluster_dict)

			updatemenus = list([dict(type='buttons',buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 1 Repetitions',method ='relayout',args = ['shapes', cluster_hap1]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'All',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap1+cluster_hap2])]))])
			layout = dict(title='Repetitions in '+c.chromosome+' '+str(min_)+'-'+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title='Genomic Position'),updatemenus=updatemenus)

		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(c.OUT+'/'+c.label+'.html'),auto_open=False)


def run(parser, args):

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] TRiCoLOR ApP v1.1')

	'''
	Check arguments, run functions
	'''

	#fill container

	c.BAM=[os.path.abspath(x) for x in args.bamfile[0]]
	c.OUT=os.path.abspath(args.output)
	c.REF=os.path.abspath(args.genome)
	c.genomebed=args.genomebed
	c.haplotype1bed=args.haplotype1bed
	c.haplotype2bed=args.haplotype2bed
	c.label=args.region

	#fcheck REGION format

	r=re.compile('.*:.*-.*')

	if r.match(args.region) is not None:

		c.chromosome=args.region.split(':')[0]
		c.start=int(args.region.split(':')[1].split('-')[0])
		c.end=int(args.region.split(':')[1].split('-')[1])

	else:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] Invalid region string format')

	
	if not os.path.exists(c.OUT):

		try:

			os.makedirs(c.OUT)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Cannot create the output folder')
			sys.exit(1)

	else:

		if not os.access(os.path.abspath(c.OUT),os.W_OK):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Missing write permissions on the output folder')
			sys.exit(1)
			
		elif os.path.exists(os.path.abspath(c.OUT+'/'+c.label+'.html')):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] The output folder already contains this plot: specify another output folder or remove the previous one')
			sys.exit(1)

	if len(c.BAM) != 2:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] ApP strictly requires a couple of splitted consensus BAM haplotypes generated by REFER')
		sys.exit(1)

	#this is probably redundat, as BAM are generated by REFER. 

	for bam in c.BAM:

		try:

			pysam.quickcheck(bam)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Errror] BAM ' + bam + ' does not exist, is not readable or is not a valid BAM')
			sys.exit(1)

		if not os.path.exists(bam + '.bai'):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Warning] Missing ' + bam + ' index. Creating')

			try:

				pysam.index(bam)

			except:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Errror] BAM ' + bam + ' could not be indexed')
				sys.exit(1)	

	try:

		ref=pyfaidx.Fasta(c.REF)

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Errror] Reference file does not exist, is not readable or is not a valid FASTA')
		sys.exit(1)


	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Plotting')

	try:

		Generate_Alignment_ToPlot(c)

	except Exception as e:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] Unexpected error while plotting: ' + str(e))
		sys.exit(1)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Done')
	sys.exit(0)