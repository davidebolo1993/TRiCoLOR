#!/usr/bin/python


import os
import pysam
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly import tools


MyDir = "/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2Mapping"

BamList=[]
ID=[]

for file in os.listdir(MyDir):
    if file.endswith(".srt.bam"):
        BamList.append(os.path.join(MyDir, file))
        ID.append(os.path.basename(os.path.join(MyDir, file)))


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
    barmode='stack'
)

MapUnmap = go.Figure(data=data, layout=BarLay)

plot(MapUnmap, filename="BarChart.html")

MeanLength=[3000,8000,13000,3000,8000,13000,3000,8000,13000,10000,100000]

MeanAccuracy=[0.85,0.85,0.85,0.90,0.90,0.90,0.95,0.95,0.95,0.90,0.98]


TableInfoNew = go.Table(
    header=dict(values=['ID','MeanLength','MeanAccuracy', 'Mapped', 'Unmapped', 'Total'],
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

