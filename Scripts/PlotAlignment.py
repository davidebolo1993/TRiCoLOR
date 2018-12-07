#!/usr/bin/python

import os
import glob
import pyfaidx
import pysam
from collections import defaultdict
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot
import argparse



def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Interactive Viewer for main TRiCoLOR results''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute( EMBL/EBI)''')
	parser.add_argument("-g", "--genome", metavar='', help="reference fasta")
	parser.add_argument("-mb1", "--merged_bam1",metavar='', help=".srt.bam file created by merging all the  consensus .srt.bam file generated for haplotype 1")
	parser.add_argument("-mb2","--merged_bam2", metavar='', help=".srt.bam file created by merging all the  consensus .srt.bam file generated for haplotype 2")
	parser.add_argument("-chr", "--chromosome",metavar='', help="chromosome number")
	parser.add_argument("-s", "--start", type=int, metavar='',help="start coordinate of the repetition you are interested to look at. Even the neighbors repetitions are plotted")
	parser.add_argument("-e","--end", type=int, metavar='',help="end coordinate of the repetition you are interested to look at. Even the neighbors repetitions are plotted")
	parser.add_argument("-rt", "--ref_table", metavar='',default=None, help=".tsv file containing repetitions found in reference")
	parser.add_argument("-b1t", "--bam1_table", metavar='', default=None, help=".tsv file containing repetitions found in haplotype 1")
	parser.add_argument("-b2t", "--bam2_table", metavar='', default=None, help=".tsv file containing repetitions found in haplotype 2")
	parser.add_argument("-l", "--label", metavar='', help="label to identify the plot")
	parser.add_argument("-o", "--out", metavar='', help="where to save the .html file")
	args = parser.parse_args()

	Generate_Alignment_ToPlot(args.genome,args.merged_bam1,args.merged_bam2,args.chromosome,args.start,args.end,args.ref_table,args.bam1_table,args.bam2_table,args.label,args.out)



def Get_Alignment_Positions(bamfile,chromosome,start,end):

	BamFile=pysam.AlignmentFile(bamfile,"rb")

	for read in BamFile.fetch(chromosome, start, end):

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
	start = next(ele for ele in coordinates if ele is not None)

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
			coords_purified.extend([None]*(coords_without_insertions[i+1]-coords_without_insertions[i]-1))
			NewSeq+=seq[i]
			NewSeq+="-"*(coords_without_insertions[i+1]-coords_without_insertions[i]-1)

		else:

			coords_purified.append(coords_without_insertions[i])
			NewSeq+=seq[i]

	return coords_purified,NewSeq



def Generate_Alignment_ToPlot(reference_fasta,hap1_bam,hap2_bam,chromosome,start,end,ref_table,hap1_table,hap2_table,label,out):

	# deal with Haplo1 and Haplo2
	bam1_coords,bam1_seq=Get_Alignment_Positions(hap1_bam,chromosome,start,end)
	bam1_mod_coords,bam1_mod_seq=Modifier(bam1_coords,bam1_seq)

	bam2_coords,bam2_seq=Get_Alignment_Positions(hap2_bam,chromosome,start,end)
	bam2_mod_coords,bam2_mod_seq=Modifier(bam2_coords,bam2_seq)


	#get minimum and maximum for the bams coordinates, so that the reference sequence is rescaled to that range

	mins=[((next(item for item in bam1_mod_coords if item is not None)),(next(item for item in bam2_mod_coords if item is not None)))][0]
	maxs=[((next(item for item in reversed(bam1_mod_coords) if item is not None)),(next(item for item in reversed(bam2_mod_coords) if item is not None)))][0]

	min_=min(mins)
	max_=max(maxs)

	ref=pyfaidx.Fasta(reference_fasta)
	chrom=ref[chromosome]
	ref_seq=chrom[:len(chrom)].seq
	ref_seq=ref_seq[min_-1:max_]

	Reference_Trace = go.Scatter(
		x = list(range(min_,max_+1)), # works with positions
		y = [0]*len(ref_seq),
		name = "Reference",
		text = list(ref_seq),
		hoverinfo = 'text+name+x',
		yaxis = "y",
		mode="markers+lines",
		line=dict(dash="dash",color='rgba(116, 116, 116, 1)'),
		marker=dict(color='rgba(116, 116, 116, 1)')
	)

	Haplo1_Trace_Upper = go.Scatter(
		x = bam1_mod_coords, # works with positions
		y = [1]*len(bam1_mod_seq),
		name = "Haplotype 1",
		text = list(bam1_mod_seq),
		hoverinfo = 'text+name+x',
		yaxis = "y",
		mode="markers+lines",
		line=dict(dash="dash",color='rgb(66, 96, 154)'),
		marker=dict(color='rgb(66, 96, 154)')
	)

	Haplo2_Trace_Lower = go.Scatter(
		x = bam2_mod_coords, # works with positions
		y = [-1]*len(bam2_mod_seq),
		name = "Haplotype 2",
		text = list(bam2_mod_seq),
		hoverinfo = 'text+name+x',
		yaxis = "y",
		mode="markers+lines",
		line=dict(dash="dash",color='rgb(152, 77, 77)'),
		marker=dict(color='rgb(152, 77, 77)')
	)


	data = [Reference_Trace,Haplo1_Trace_Upper,Haplo2_Trace_Lower]

	if ref_table is None and hap1_table is None and hap2_table is None:

		layout = dict(title='Repetitions in '+chromosome+" "+str(min_)+"-"+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title="Genomic Position"))
		fig = dict(data=data, layout=layout)
		plot(fig,filename=os.path.abspath(out)+"/"+label+".html",auto_open=False)


	elif ref_table is not None and hap1_table is None and hap2_table is None:

		ref_rep=pd.read_csv(os.path.abspath(ref_table),sep="\t")
		cluster_ref=[]

		for i in range(len(ref_rep["Start"])):

			if ref_rep["Start"][i] >= min_ and ref_rep["End"][i] <= max_:

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep["Start"][i],y0=-.15,x1=ref_rep["End"][i],y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
				cluster_ref.append(cluster_dict)


		updatemenus = list([dict(type="buttons",buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref])]))])
		layout = dict(title='Repetitions in '+chromosome+" "+str(min_)+"-"+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title="Genomic Position"),updatemenus=updatemenus)


	elif ref_table is None and hap1_table is not None and hap2_table is None:


		hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep="\t")
		cluster_hap1=[]

		for i in range(len(hap1_rep["Start"])):

			if hap1_rep["Start"][i] >= min_ and hap1_rep["End"][i] <= max_: 

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep["Start"][i],y0=.85,x1=hap1_rep["End"][i],y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
				cluster_hap1.append(cluster_dict)

		updatemenus = list([dict(type="buttons",buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 1 Repetitions',method = 'relayout',args = ['shapes', cluster_hap1])]))])
		layout = dict(title='Repetitions in '+chromosome+" "+str(min_)+"-"+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title="Genomic Position"),updatemenus=updatemenus)


	elif ref_table is None and hap1_table is None and hap2_table is not None:


		hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep="\t")
		cluster_hap2=[]

		for i in range(len(hap2_rep["Start"])):

			if hap2_rep["Start"][i] >= min_ and hap2_rep["End"][i] <= max_:

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep["Start"][i],y0=-.15,x1=hap2_rep["End"][i],y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
				cluster_ref.append(cluster_dict)

		updatemenus = list([dict(type="buttons",buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 2 Repetitions',method = 'relayout',args = ['shapes', cluster_hap2])]))])
		layout = dict(title='Repetitions in '+chromosome+" "+str(min_)+"-"+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title="Genomic Position"),updatemenus=updatemenus)


	elif ref_table is not None and hap1_table is not None and hap2_table is None:

		ref_rep=pd.read_csv(os.path.abspath(ref_table),sep="\t")
		cluster_ref=[]

		for i in range(len(ref_rep["Start"])):

			if ref_rep["Start"][i] >= min_ and ref_rep["End"][i] <= max_:

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep["Start"][i],y0=-.15,x1=ref_rep["End"][i],y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
				cluster_ref.append(cluster_dict)


		hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep="\t")
		cluster_hap1=[]

		for i in range(len(hap1_rep["Start"])):

			if hap1_rep["Start"][i] >= min_ and hap1_rep["End"][i] <= max_: 

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep["Start"][i],y0=.85,x1=hap1_rep["End"][i],y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
				cluster_hap1.append(cluster_dict)


		updatemenus = list([dict(type="buttons",buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 1 Repetitions',method ='relayout',args = ['shapes', cluster_hap1]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap1])]))])
		layout = dict(title='Repetitions in '+chromosome+" "+str(min_)+"-"+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title="Genomic Position"),updatemenus=updatemenus)


	elif ref_table is not None and hap1_table is None and hap2_table is not None:

		ref_rep=pd.read_csv(os.path.abspath(ref_table),sep="\t")
		cluster_ref=[]

		for i in range(len(ref_rep["Start"])):

			if ref_rep["Start"][i] >= min_ and ref_rep["End"][i] <= max_:

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep["Start"][i],y0=-.15,x1=ref_rep["End"][i],y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
				cluster_ref.append(cluster_dict)


		hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep="\t")
		cluster_hap2=[]

		for i in range(len(hap2_rep["Start"])):

			if hap2_rep["Start"][i] >= min_ and hap2_rep["End"][i] <= max_:

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep["Start"][i],y0=-1.15,x1=hap2_rep["End"][i],y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
				cluster_hap2.append(cluster_dict)

		updatemenus = list([dict(type="buttons",buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap2])]))])
		layout = dict(title='Repetitions in '+chromosome+" "+str(min_)+"-"+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title="Genomic Position"),updatemenus=updatemenus)


	elif ref_table is None and hap1_table is not None and hap2_table is not None:

		hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep="\t")
		cluster_hap1=[]

		for i in range(len(hap1_rep["Start"])):

			if hap1_rep["Start"][i] >= min_ and hap1_rep["End"][i] <= max_: 

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep["Start"][i],y0=.85,x1=hap1_rep["End"][i],y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
				cluster_hap1.append(cluster_dict)


		hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep="\t")
		cluster_hap2=[]

		for i in range(len(hap2_rep["Start"])):

			if hap2_rep["Start"][i] >= min_ and hap2_rep["End"][i] <= max_:

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep["Start"][i],y0=-1.15,x1=hap2_rep["End"][i],y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
				cluster_hap2.append(cluster_dict)

		updatemenus = list([dict(type="buttons",buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Haplotype 1 Repetitions',method = 'relayout',args = ['shapes', cluster_hap1]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'Both',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap2])]))])
		layout = dict(title='Repetitions in '+chromosome+" "+str(min_)+"-"+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title="Genomic Position"),updatemenus=updatemenus)

	else:

		ref_rep=pd.read_csv(os.path.abspath(ref_table),sep="\t")
		cluster_ref=[]

		for i in range(len(ref_rep["Start"])):

			if ref_rep["Start"][i] >= min_ and ref_rep["End"][i] <= max_:

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=ref_rep["Start"][i],y0=-.15,x1=ref_rep["End"][i],y1=.15, opacity=.25,line=dict(color='rgb(0, 190, 110)'),fillcolor='rgb(0, 190, 110)')
				cluster_ref.append(cluster_dict)

		hap1_rep=pd.read_csv(os.path.abspath(hap1_table),sep="\t")
		cluster_hap1=[]

		for i in range(len(hap1_rep["Start"])):

			if hap1_rep["Start"][i] >= min_ and hap1_rep["End"][i] <= max_: 

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap1_rep["Start"][i],y0=.85,x1=hap1_rep["End"][i],y1=1.15, opacity=.25,line=dict(color='#d10cf1'),fillcolor='#d10cf1')
				cluster_hap1.append(cluster_dict)


		hap2_rep=pd.read_csv(os.path.abspath(hap2_table),sep="\t")
		cluster_hap2=[]

		for i in range(len(hap2_rep["Start"])):

			if hap2_rep["Start"][i] >= min_ and hap2_rep["End"][i] <= max_:

				cluster_dict=dict(type='rectangule',xref='x',yref='y',x0=hap2_rep["Start"][i],y0=-1.15,x1=hap2_rep["End"][i],y1=-.85, opacity=.25,line=dict(color='#f1c31f'),fillcolor='#f1c31f')
				cluster_hap2.append(cluster_dict)

		updatemenus = list([dict(type="buttons",buttons=list([dict(label = 'None',method = 'relayout',args = ['shapes', []]),dict(label = 'Reference Repetitions',method = 'relayout',args = ['shapes', cluster_ref]),dict(label = 'Haplotype 1 Repetitions',method ='relayout',args = ['shapes', cluster_hap1]),dict(label = 'Haplotype 2 Repetitions',method ='relayout',args = ['shapes', cluster_hap2]),dict(label = 'All',method = 'relayout',args = ['shapes', cluster_ref+cluster_hap1+cluster_hap2])]))])
		layout = dict(title='Repetitions in '+chromosome+" "+str(min_)+"-"+str(max_),yaxis=dict(range=[-4,+4],showticklabels=False,ticks=''),xaxis=dict(title="Genomic Position"),updatemenus=updatemenus)


	fig = dict(data=data, layout=layout)
	plot(fig,filename=os.path.abspath(out)+"/"+label+".html",auto_open=False)



if __name__ == main():

	main()
