#!/usr/bin/python

import os
import re
import pysam
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot


def RepeatsFinder(string,kmer,times,strict):

    if kmer != 0 and times == 0:

        my_regex = r"(.{"  + re.escape(str(kmer)) + r"})\1+"

    elif kmer != 0 and times != 0:

        my_regex = r"(.{"  + re.escape(str(kmer)) + r"})\1{" + re.escape(str(times)) + r"}" if strict==0 else r"(.{"  + re.escape(str(kmer)) + r"})\1{" + re.escape(str(times-1)) + r",}"

    elif kmer == 0 and times != 0:

        my_regex = r"(.+?)\1{" + re.escape(str(kmer)) + r"}" if strict==0 else r"(.+?)\1{" + re.escape(str(times-1)) + r",}"

    elif kmer == 0 and times == 0:

        my_regex = r"(.+?)\1+"

    else:

        sys.exit('Allowed kmers/times arguments are any integer')

    r=re.compile(my_regex)

    for match in r.finditer(string):

        yield (match.group(1), match.start(1), match.start(1)+len(match.group(1))*int(len(match.group(0))/len(match.group(1)))-1,int(len(match.group(0))/len(match.group(1))))



def ResultsWriter_FromRef(repetitions,out,label,realstart):

    seq=[el[0] for el in repetitions]
    start=[realstart+el[1] for el in repetitions]
    end=[realstart+el[2] for el in repetitions]
    rep=[el[3] for el in repetitions]

    os.chdir(os.path.abspath(out))

    Table=pd.DataFrame({"Core_Repetition":seq, "Start":start,"End":end, "Repetitions":rep},columns=['Core_Repetition', 'Start', 'End', 'Repetitions'])
    Table.to_csv(label + "_RepetitionsTable.tsv",sep="\t",index=False)
    BarWidth=Table["End"]-Table["Start"]
    AnnotationsText=Table["Core_Repetition"].tolist()
    AnnotationXPos=Table["Start"].tolist()
    AnnotationYPos=Table["Repetitions"].tolist()
    HoverInfoStart=list(map(str,Table["Start"].tolist()))
    HoverInfoEnd=list(map(str,Table["End"].tolist()))
    HoverInfo=["Repetition Interval: "+a+"-"+b for a,b in zip(HoverInfoStart,HoverInfoEnd)]
    Trace0=go.Bar(x=Table["Start"], y=Table["Repetitions"], width=BarWidth,marker=dict(color='rgba(204,204,204,1)'),text=HoverInfo,hoverinfo="text",hoverlabel=dict(bgcolor="white",bordercolor="black"))
    TraceLayout=go.Layout(title= "Repeated Sequences", xaxis=dict(range=[start[0]-1000,start[-1]+1000],title='Base number in consensus'),yaxis=dict(title='Number of repetitions'), showlegend=False)
    annotations=[]

    for i in range(len(AnnotationsText)):

        annotations.append(dict(x=AnnotationXPos[i],y=AnnotationYPos[i], text=AnnotationsText[i], showarrow=True, arrowhead=5))

    TraceLayout["annotations"]=annotations
    fig=go.Figure([Trace0],TraceLayout)
    plot(fig,filename=label+"_RepetitionsBarchart.html",auto_open=False)



def Get_Alignment_Position(bamfile):

	BamFile=pysam.AlignmentFile(bamfile,"rb")

	for read in BamFile.fetch():

		coords = read.get_reference_positions(full_length=True)

		seq=read.seq

	return coords,seq


def Get_Alignment_Coordinates(coord_list,alignment_file,repetitions): #works for inserted and matched repetitions: think it works even for deleted repetitions; the following coordinated are reporetd 1-based

	Start=[]
	End=[]

	for i in range(len(alignment_file)):

		if i == 0: #for soft clipped at the left edge

			if alignment_file[i] == "-":

				pass

			else:

				if coord_list[i] is not None:

					Start.append(coord_list[i]+1)

				else:

					l=i

					while coord_list[l] is None:

						l+=1

					Start.append(coord_list[l]+1-l)

		else:

			if alignment_file[i] == "-":

				pass

			else:

				if alignment_file[i-1] == "-" and alignment_file[i+1] != "-":

					if coord_list[i] is not None:

						Start.append(coord_list[i]+1)

					else: #this bases are soft clipped or inserted

						l=i

						while coord_list[l] is None:

							l -= 1

						Start.append(coord_list[l]+1)


				elif alignment_file[i-1] != "-" and alignment_file [i+1] != "-": #it's not start or end

					pass

				else:

					if coord_list[i] is not None:

						End.append(coord_list[i]+1)

					else: #this bases are soft clipped or inserted

						l=i

						while coord_list[l] is None:

							l-= 1

						End.append(coord_list[l]+1)

	rep_file=[(a,b,c,d) for a,b,c,d in zip([repetition[0] for repetition in repetitions],Start,End,[repetition[3] for repetition in repetitions])]

	return rep_file


def Consensus_Real_Start(coord):

    if coord[0] is not None:

        return coord[0]

    else:

        counter=0

        while coord[counter] is None:

            counter +=1

        return coord[counter]+1-counter



def Consensus_Alignment_Generator(realstart, string, repetitions):


    dim=realstart+len(string)

    seq=[el[0] for el in repetitions]
    start=[realstart+el[1] for el in repetitions]
    end=[realstart+el[2] for el in repetitions]
    rep=[el[3] for el in repetitions]

    i=realstart
    j=0

    NewSeq=''

    while i < dim:

        if i in start:

            PartialSeq=seq[j]*rep[j]
            NewSeq+=PartialSeq
            i+=len(PartialSeq)
            j+=1

        else:

            NewSeq+="-"
            i+=1

    return NewSeq



def ResultsWriter_FromCons(repetitions_with_coord,out,label): #can be softclipped due to a repeat insertion

    seq=[el[0] for el in repetitions_with_coord]
    start=[el[1] for el in repetitions_with_coord]
    end=[el[2] for el in repetitions_with_coord]
    rep=[el[3] for el in repetitions_with_coord]

    os.chdir(os.path.abspath(out))

    Table=pd.DataFrame({"Core_Repetition":seq, "Start":start,"End":end, "Repetitions":rep},columns=['Core_Repetition', 'Start', 'End', 'Repetitions'])
    Table.to_csv(label + "_RepetitionsTable.tsv",sep="\t",index=False)
    BarWidth=Table["End"]-Table["Start"]
    AnnotationsText=Table["Core_Repetition"].tolist()
    AnnotationXPos=Table["Start"].tolist()
    AnnotationYPos=Table["Repetitions"].tolist()
    HoverInfoStart=list(map(str,Table["Start"].tolist()))
    HoverInfoEnd=list(map(str,Table["End"].tolist()))
    HoverInfo=["Repetition Interval: "+a+"-"+b for a,b in zip(HoverInfoStart,HoverInfoEnd)]
    Trace0=go.Bar(x=Table["Start"], y=Table["Repetitions"],text=HoverInfo,hoverinfo="text",hoverlabel=dict(bgcolor="white",bordercolor="black"))

    width=[]
    colors=[]
    for i in range(len(BarWidth)):

        if BarWidth[i] != 0:

            width.append(BarWidth[i])
            colors.append('rgba(204,204,204,1)')

        else:

            width.append(Table["Repetitions"][i])
            colors.append('rgba(222,45,38,0.8)')


    Trace0["width"]=width
    Trace0["marker"]=dict(color=colors)


    TraceLayout=go.Layout(title= "Repeated Sequences", xaxis=dict(range=[start[0]-1000,start[-1]+1000],title='Base number in consensus'),yaxis=dict(title='Number of repetitions'), showlegend=False)
    annotations=[]

    for i in range(len(AnnotationsText)):

        annotations.append(dict(x=AnnotationXPos[i],y=AnnotationYPos[i], text=AnnotationsText[i], showarrow=True, arrowhead=5))

    TraceLayout["annotations"]=annotations
    fig=go.Figure([Trace0],TraceLayout)
    plot(fig,filename=label+"_RepetitionsBarchart.html",auto_open=False)
