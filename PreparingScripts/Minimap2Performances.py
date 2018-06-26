#!/usr/bin/python


import os
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly import tools


File = "/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/PyInfo.txt"
Names = ["File_Name", "Number_of_Sequences", "Minimap_Time", "Mean_Length", "Min_Length", "Max_Length"]


Table = pd.read_table(File, delimiter='\t',names=Names)
PyTable = Table.sort_values('Number_of_Sequences')


MinimapTime = go.Scatter(
    x=PyTable['Number_of_Sequences'],
    y=PyTable['Minimap_Time'],
    mode = 'lines+markers',
    text=PyTable['File_Name'],
    fill='tozeroy'
)

MinimapTimeLayout = go.Layout(
    title= '#Alignment time',
    xaxis= dict(title='#Time (minutes)'),
    yaxis=dict(title='#Number of Reads', side='right'),
)

MiniPlot = [MinimapTime]
MinimapFigure = go.Figure(data=MiniPlot, layout=MinimapTimeLayout)
plot(MinimapFigure,filename="Time.html")


MeanLength = go.Scatter(	
    x=PyTable['Number_of_Sequences'],
    y=PyTable['Mean_Length'],
    mode = 'lines+markers',
    text=PyTable['File_Name'],
    fill='tozeroy',
    name='Mean'

)
MinLength = go.Scatter(
    x=PyTable['Number_of_Sequences'],
    y=PyTable['Min_Length'],
    mode = 'lines+markers',
    text=PyTable['File_Name'],
    fill='tozeroy',
    name='Min'
)
MaxLength = go.Scatter(
    x=PyTable['Number_of_Sequences'],
    y=PyTable['Max_Length'],
    mode = 'lines+markers',
    text=PyTable['File_Name'],
    fill='tozeroy',
    name='Max'

)

LengthPlot = [MeanLength, MinLength, MaxLength]
LengthPlotLayout = go.Layout(
    title= '#Sequence min/mean/max length',
    xaxis= dict(title='#Number of Reads'),
    yaxis=dict(title='#Length'),
)
LengthFigure = go.Figure(data=LengthPlot, layout=LengthPlotLayout)
plot(LengthFigure,filename="Length.html" )


MyDir = "/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2Mapping"

BamList=[]
ID=[]

for file in os.listdir(MyDir):
    if file.endswith(".srt.bam"):
        BamList.append(os.path.join(MyDir, file))
        IDLab = os.path.basename(os.path.join(MyDir, file))
        NewIdLab=IDLab.split(".")[0]
        ID.append(NewIdLab)

Total = []
MapList = []
UnMapList = []


for Bam in BamList:
	TotalNumber = pysam.view("-c",Bam) 
	bamfile = pysam.AlignmentFile(Bam, "rb")
	Total.append(int(TotalNumber))
	MapList.append(int(bamfile.mapped)) # == samtools view -F 0x40 "srt.bam" | cut -f1 | sort | wc -l
	UnMapList.append(int(bamfile.unmapped)) # == samtools view -f 4 | cut -f1 | sort | wc -l



InfoTable = pd.DataFrame(
    {'ID': ID,
     'Number of sequeces in bam': Total,
     'Mapped': MapList,
     'Unmapped':UnMapList
    })


Mapped = go.Bar(
    x=InfoTable['ID'],
    y=InfoTable['Mapped'],
    name='Mapped',
    text=InfoTable['Number of sequeces in bam']
)

UnMapped = go.Bar(
    x=InfoTable['ID'],
    y=InfoTable['Unmapped'],
    name='Unmapped'
)

data = [Mapped, UnMapped]

BarLay = go.Layout(
    barmode='stack',
    title= 'Mapped/Unmapped sequences',
    xaxis= dict(title='#Read Identifier'),
    yaxis=dict(title='#Mapped/Unmapped Reads')
)


MapUnmap = go.Figure(data=data, layout=BarLay)

plot(MapUnmap, filename="BarChart.html")


MeanLength=[3000,8000,13000,3000,8000,13000,3000,8000,13000,10000,100000] #data used for create syntetic fastq with pbsim

MeanAccuracy=[0.85,0.85,0.85,0.90,0.90,0.90,0.95,0.95,0.95,0.90,0.98] #data used for create syntetic fastq with pbsim


TableInfoNew = go.Table(
    header=dict(values=['#ID','#Mean Length','#Mean Accuracy', '#Mapped Reads', '#Unmapped Reads', '#Total Reads'],
                line = dict(color='#7D7F80'),
                fill = dict(color='#a1c3d1'),
                align = ['left'] * 5),
    cells=dict(values=[InfoTable['ID'],
    	               MeanLength,
                       MeanAccuracy,
                       InfoTable['Mapped'],
                       InfoTable['Unmapped'],
                       InfoTable['Number of sequeces in bam']
                       ],
               line = dict(color='#7D7F80'),
               fill = dict(color='#EDFAFF'),
               align = ['left'] * 5))

PlotTable = [TableInfoNew]

plot(PlotTable, filename = 'InfoTable.html')

