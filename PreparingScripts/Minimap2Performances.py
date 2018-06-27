#!/usr/bin/python


import os
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly import tools
import pysam

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
    title= 'Alignment time',
    xaxis= dict(title='# Reads'),
    yaxis=dict(title='Time (minutes)', side='left'),
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
    title= 'Sequence min/mean/max length',
    xaxis= dict(title='# Reads'),
    yaxis=dict(title='Length'),
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
TotalMap = []
TotalUnmap = []
MapInRegionList = []
UnmapInRegionList = []


for Bam in BamList:
	TotalAligned = pysam.view("-c",Bam) 
	Total.append(int(TotalAligned))
	bamfile = pysam.AlignmentFile(Bam, 'rb')
	TotalMap.append(int(bamfile.mapped)) # == samtools view -F 0x40 "srt.bam" | cut -f1 | sort | wc -l
	TotalUnmap.append(int(bamfile.unmapped)) # == samtools view -f 4 | cut -f1 | sort | wc -l
	mapinregion = 0
	unmapinregion = 0
	for read in bamfile.fetch('chr1',1000001,40000000):
		if not read.is_unmapped:
			mapinregion+=1
		else:
			unmapinregion+=1
	MapInRegionList.append(mapinregion)
	UnmapInRegionList.append(unmapinregion)


InfoTable = pd.DataFrame(
    {'ID': ID,
     '# Mapped in Bam': TotalMap,
     '# Unmapped in Bam':TotalUnmap,
     '# Mapped in region': MapInRegionList,
     '# Unmapped in region': UnmapInRegionList
    })

InfoTable = InfoTable [['ID','# Mapped in Bam','# Unmapped in Bam','# Mapped in region','# Unmapped in region']]

TotMapped = go.Bar(
    x=InfoTable['ID'],
    y=InfoTable['# Mapped in Bam'],
    name='# Mapped in Bam',
)

TotUnMapped = go.Bar(
    x=InfoTable['ID'],
    y=InfoTable['# Unmapped in Bam'],
    name='# Unmapped in Bam'
)

MappedInRegion = go.Bar(
    x=InfoTable['ID'],
    y=InfoTable['# Mapped in region'],
    name='# Mapped in region'
)


UnMappedInRegion = go.Bar(
    x=InfoTable['ID'],
    y=InfoTable['# Unmapped in region'],
    name='# Unmapped in region'
)


data = [TotMapped, TotUnMapped, MappedInRegion, UnMappedInRegion]

BarLay = go.Layout(
    barmode='base-bar',
    title= '# Mapped/Unmapped reads',
    xaxis= dict(title='Read Identifier'),
    yaxis=dict(title='# Mapped and Unmapped Reads')
)


MapUnmap = go.Figure(data=data, layout=BarLay)

plot(MapUnmap, filename="BarChart.html")


MeanLength=[3000,8000,13000,3000,8000,13000,3000,8000,13000,10000,100000] #data used for create syntetic fastq with pbsim

MeanAccuracy=[0.85,0.85,0.85,0.90,0.90,0.90,0.95,0.95,0.95,0.90,0.98] #data used for create syntetic fastq with pbsim

MapFraction = InfoTable['# Mapped in Bam']/(InfoTable['# Mapped in Bam']+InfoTable['# Unmapped in Bam'])
RightMapped = InfoTable['# Mapped in region']/InfoTable['# Mapped in Bam']


TableInfoNew = go.Table(
    header=dict(values=['ID','Mean Length','Mean Accuracy', '# Mapped Reads in Bam', '# Unmapped Reads in Bam', '# Mapped Reads in Region', '# UnMapped Reads in Region', '# Mapped Reads in Bam / #Total', '# Mapped Reads in Region / # Mapped Reads in Bam'],
                line = dict(color='#7D7F80'),
                fill = dict(color='#a1c3d1'),
                align = ['left'] * 5),
    cells=dict(values=[InfoTable['ID'],
    	               MeanLength,
                       MeanAccuracy,
                       InfoTable['# Mapped in Bam'],
                       InfoTable['# Unmapped in Bam'],
                       InfoTable['# Mapped in region'],
                       InfoTable['# Unmapped in region'],
                       ["{0:.4f}".format(i) for i in MapFraction],
                       ["{0:.4f}".format(i) for i in RightMapped]
                       ],
               line = dict(color='#7D7F80'),
               fill = dict(color='#EDFAFF'),
               align = ['left'] * 5))

PlotTable = [TableInfoNew]

plot(PlotTable, filename = 'InfoTable.html')

