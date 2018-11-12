#!/usr/bin/python

import sys
import os
import re
import pyfaidx
import glob
import subprocess
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot
from RepFinder import *
from BamParser import *
import gzip
import shutil


def Fasta_Writer(header,seq,out):

	dirname=os.path.abspath(out)

	if not os.path.exists(dirname):

		os.makedirs(dirname)

	fname=os.path.abspath(dirname+"/ReferenceRegion.fa")
	f=open(fname,"w")
	f.write(header + "\n" + seq)
	f.close()


def Ref_Repeats(reference, chromosome, start, end, kmer, times, strict, out):

	out_=out+"/Reference"

	if not os.path.exists(out_):

		os.makedirs(out_)

	ref=pyfaidx.Fasta(reference)
	chrom=ref[chromosome]
	seq=chrom[:len(chrom)].seq
	wanted=seq[start-1:end]
	header=">"+chromosome+":"+str(start)+"-"+str(end)
	Fasta_Writer(header,wanted,out)
	rep=list(RepeatsFinder(wanted,kmer,times,strict))
	ResultsWriter_FromRef(rep,out_,"Reference",start)



def Haplo1_Repeats(bamfile1, chromosome, start, end, size, out, kmer, times, strict,mmi_ref):

	out_=out+"/haplotype1"

	if not os.path.exists(out_):

		os.makedirs(out_)

	seq,coord = Bamfile_Analyzer(bamfile1,chromosome,start,end,size)
	filseq,filcoord=Com_Seq(seq,coord)
	Fasta_Generator(filseq,filcoord,out_)
	MA(out_,mmi_ref)
	consensus_bams=glob.glob(os.path.abspath(out_)+"/*.srt.bam")
	consensus_labels=[name.split(".")[0] for name in consensus_bams]

	for bam,label in  zip(consensus_bams,consensus_labels):

		coords,seq=Get_Alignment_Position(bam)
		bam_start=Consensus_Real_Start(coords)
		repetitions=list(RepeatsFinder(seq,kmer,times,strict))
		alignment_file=Consensus_Alignment_Generator(bam_start, seq, repetitions)
		rep_coord_list=Get_Alignment_Coordinates(coords,alignment_file,repetitions)
		ResultsWriter_FromCons(rep_coord_list,out_,label)



def Haplo2_Repeats(bamfile2, chromosome, start, end, size, out, kmer, times, strict,mmi_ref):

	out_=out+"/haplotype2"

	if not os.path.exists(out_):

		os.makedirs(out_)

	seq,coord = Bamfile_Analyzer(bamfile2,chromosome,start,end,size)
	filseq,filcoord=Com_Seq(seq,coord)
	Fasta_Generator(filseq,filcoord,out_)
	MA(out_,mmi_ref)
	consensus_bams=glob.glob(os.path.abspath(out_)+"/*.srt.bam")
	consensus_labels=[name.split(".")[0] for name in consensus_bams]

	for bam,label in  zip(consensus_bams,consensus_labels):

		coords,seq=Get_Alignment_Position(bam)
		bam_start=Consensus_Real_Start(coords)
		repetitions=list(RepeatsFinder(seq,kmer,times,strict))
		alignment_file=Consensus_Alignment_Generator(bam_start, seq, repetitions)
		rep_coord_list=Get_Alignment_Coordinates(coords,alignment_file,repetitions)
		ResultsWriter_FromCons(rep_coord_list,out_,label)



def Concat_Tables(list_of_paths):

	List_of_tables=[]

	for tab in list_of_paths:

		Tab=pd.read_csv(tab, sep="\t")
		List_of_tables.append(Tab)

	By_Row=pd.concat(List_of_tables, axis=0,ignore_index=True)

	return By_Row


def isNaN(value):

    return value != value


def Compare_Results(out):

	os.chdir(os.path.abspath(out))

	out_=[os.path.abspath(out+j) for j in ["/Reference", "/haplotype1", "/haplotype2"]]

	Reference_Tsv=glob.glob(os.path.abspath(out_[0])+"/*.tsv")
	Haplo_One_Tsvs=glob.glob(os.path.abspath(out_[1])+"/*.tsv")
	Haplo_Two_Tsvs=glob.glob(os.path.abspath(out_[2])+"/*.tsv")

	Reference_Tab=pd.read_csv(Reference_Tsv[0], sep="\t")
	Haplo_One_Tab=Concat_Tables(Haplo_One_Tsvs)
	Haplo_Two_Tab=Concat_Tables(Haplo_Two_Tsvs)


	#Compare Haplotyped Bam Repetitions

	MergeHaploTab=Haplo_One_Tab.merge(Haplo_Two_Tab, how='outer', indicator=True)
	MergeHaploTab.sort_values(by=["Start"], inplace=True)
	MergeHaploTab.replace(to_replace={"_merge":{"left_only":"haplotype_1", "right_only":"haplotype_2"}}, inplace=True)
	MergeHaploTab.rename({"_merge": "haplo_differences"}, axis='columns', inplace=True) #rename to avoid problem with the next indicator


	MergedRefHaplo=MergeHaploTab=Reference_Tab.merge(MergeHaploTab, how='outer', indicator=True)
	MergedRefHaplo.sort_values(by=["Start"], inplace=True) #don't know if it's always needed


	MergedRefHaplo.fillna(0,inplace=True)


	FinalCol=[]

	for i in range(len(MergedRefHaplo["_merge"])):

		if not isNaN(MergedRefHaplo["haplo_differences"][i]):

			if MergedRefHaplo["haplo_differences"][i] == "both" and  MergedRefHaplo["_merge"][i] == "both":

				FinalCol.append("all")

			elif MergedRefHaplo["haplo_differences"][i] != "both" and  MergedRefHaplo["_merge"][i] == "both":

				FinalCol.append("reference and "+str(MergedRefHaplo["haplo_differences"][i]))

			elif MergedRefHaplo["haplo_differences"][i] != "both" and  MergedRefHaplo["_merge"][i] != "both":

				if  MergedRefHaplo["_merge"][i] == "left_only":

					FinalCol.append("reference")

				else:

					FinalCol.append(str(MergedRefHaplo["haplo_differences"][i]))

			elif MergedRefHaplo["haplo_differences"][i] == "both" and  MergedRefHaplo["_merge"][i] != "both":

				FinalCol.append("haplotype_1 and haplotype_2")

		else:

				FinalCol.append("reference")

	MergedRefHaplo.drop(columns=["_merge","haplo_differences"], inplace=True)
	MergedRefHaplo["Where"]=FinalCol
	MergedRefHaplo.to_csv("Comparison_RepetitionsTable.tsv",sep="\t",index=False)


	#Reference

	BarWidth=Reference_Tab["End"]-Reference_Tab["Start"]
	AnnotationsText=Reference_Tab["Core_Repetition"].tolist()
	AnnotationXPos=Reference_Tab["Start"].tolist()
	AnnotationYPos=Reference_Tab["Repetitions"].tolist()
	HoverInfoStart=list(map(str,Reference_Tab["Start"].tolist()))
	HoverInfoEnd=list(map(str,Reference_Tab["End"].tolist()))
	HoverInfo=["Reference" + "\n""Repetition Interval: "+a+"-"+b+"\n"+"Repetition: " + c for a,b,c in zip(HoverInfoStart,HoverInfoEnd,AnnotationsText)]
	Trace0=go.Bar(x=Reference_Tab["Start"], y=Reference_Tab["Repetitions"],text=HoverInfo,hoverinfo="text",hoverlabel=dict(bgcolor="white",bordercolor="black"),opacity=0.8,name="Reference",marker=dict(color='rgb(180, 180, 181)'))

	width=[]

	for i in range(len(BarWidth)):

	    if BarWidth[i] != 0:

	        width.append(BarWidth[i])

	    else:

	        width.append(Table["Repetitions"][i])

	Trace0["width"]=width


	#Haplotype 1

	BarWidth=Haplo_One_Tab["End"]-Haplo_One_Tab["Start"]
	AnnotationsText=Haplo_One_Tab["Core_Repetition"].tolist()
	AnnotationXPos=Haplo_One_Tab["Start"].tolist()
	AnnotationYPos=Haplo_One_Tab["Repetitions"].tolist()
	HoverInfoStart=list(map(str,Haplo_One_Tab["Start"].tolist()))
	HoverInfoEnd=list(map(str,Haplo_One_Tab["End"].tolist()))
	HoverInfo=["Haplotype1" + "\n" +"Repetition Interval: "+a+"-"+b+"\n"+"Repetition: " + c for a,b,c in zip(HoverInfoStart,HoverInfoEnd,AnnotationsText)]
	Trace1=go.Bar(x=Haplo_One_Tab["Start"], y=Haplo_One_Tab["Repetitions"],text=HoverInfo,hoverinfo="text",hoverlabel=dict(bgcolor="white",bordercolor="black"),opacity=0.8,name="Haplotype 1",marker=dict(color='rgb(9, 30, 218)'))

	width=[]

	for i in range(len(BarWidth)):

	    if BarWidth[i] != 0:

	        width.append(BarWidth[i])

	    else:

	        width.append(Table["Repetitions"][i])

	Trace1["width"]=width

	#Haplotype 2

	BarWidth=Haplo_Two_Tab["End"]-Haplo_Two_Tab["Start"]
	AnnotationsText=Haplo_Two_Tab["Core_Repetition"].tolist()
	AnnotationXPos=Haplo_Two_Tab["Start"].tolist()
	AnnotationYPos=Haplo_Two_Tab["Repetitions"].tolist()
	HoverInfoStart=list(map(str,Haplo_Two_Tab["Start"].tolist()))
	HoverInfoEnd=list(map(str,Haplo_Two_Tab["End"].tolist()))
	HoverInfo=["Haplotype2" + "\n" + "Repetition Interval: "+a+"-"+b+"\n"+"Repetition: " + c for a,b,c in zip(HoverInfoStart,HoverInfoEnd,AnnotationsText)]
	Trace2=go.Bar(x=Haplo_Two_Tab["Start"], y=Haplo_Two_Tab["Repetitions"],text=HoverInfo,hoverinfo="text",hoverlabel=dict(bgcolor="white",bordercolor="black"),opacity=0.8,name="Haplotype 2",marker=dict(color='rgb(189, 10, 10)'))

	width=[]

	for i in range(len(BarWidth)):

	    if BarWidth[i] != 0:

	        width.append(BarWidth[i])

	    else:

	        width.append(Table["Repetitions"][i])

	Trace2["width"]=width

	#Combine

	data=[Trace0,Trace1,Trace2]
	layout = dict(title = 'Repetitions in Reference, Haplotype 1 and Haplotype 2', yaxis=dict(range=[0,max(MergedRefHaplo["Repetitions"]+1)]), xaxis=dict(range=[min(MergedRefHaplo['Start'])-200, max(MergedRefHaplo['End'])+200]))
	fig = go.Figure(data=data, layout=layout)
	plot(fig,filename="Comparison_RepetitionsBarChart.html", auto_open=False)
