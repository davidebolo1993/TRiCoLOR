#!/usr/bin/python

import os
import numpy as np
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly import tools

MyDir = "/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/ChromosomesCoverage"

Coverage=[]
ID=[]

for file in os.listdir(MyDir):
    if file.endswith(".coverage"):
        Coverage.append(os.path.join(MyDir, file))
        IDLab = os.path.basename(os.path.join(MyDir, file))
        NewIdLab=IDLab.split(".")[0]
        ID.append(NewIdLab)


NewDir = "Python_Coverage_Plots"

os.mkdir(NewDir)
os.chdir(NewDir)


def rep(x, repeat):
	return np.repeat(x, repeat).tolist()


for i in range(len(Coverage)):
	
	Table= np.genfromtxt(Coverage[i], delimiter='\t',dtype=['S5',"i8",'i8'],names=['Chromosome','Locus','Depth'])
	Gap = np.diff(Table['Locus'])
	Indexes = np.where(Gap>1)
	Pos=[]
	MedianCov=[]
	AltLen=[]

	Posit=int(round(np.median(Table['Locus'][0:Indexes[0][0]]).tolist())) # median position for first interval
	Cover=int(round(np.median(Table['Depth'][0:Indexes[0][0]]).tolist())) # median depth for first interval
	Pos.append(Posit)
	MedianCov.append(Cover)
	AltLen.append(len(Table['Locus'][0:Indexes[0][0]]))

	for j in range(len(Indexes[0])-1):
		Posit = int(round(np.median(Table['Locus'][(Indexes[0][j]+1):Indexes[0][j+1]]).tolist())) 
		Cover = int(round(np.median(Table['Depth'][(Indexes[0][j]+1):Indexes[0][j+1]]).tolist()))
		Pos.append(Posit)
		MedianCov.append(Cover)
		AltLen.append(len(Table['Locus'][(Indexes[0][j]+1):Indexes[0][j+1]]))

	Posit = int(round(np.median(Table['Locus'][(Indexes[0][-1]+1):len(Table['Locus'])]).tolist())) # median position for last interval
	Cover = int(round(np.median(Table['Depth'][(Indexes[0][-1]+1):len(Table['Depth'])]).tolist())) # median depth for last interval
	
	Pos.append(Posit)
	MedianCov.append(Cover)
	AltLen.append(len(Table['Locus'][(Indexes[0][-1]+1):len(Table['Locus'])]))


	Title=["Coverage " + ID[i]]
	FileName=["Coverage_" + ID[i] + ".html"]


	CoveragePlot = go.Bar(x=Pos,y=MedianCov,name="Coverage in event",text=AltLen)
	MeanLine = go.Scatter(x=Pos,y=rep(np.median(MedianCov).tolist(),len(MedianCov)),mode = 'lines',name="Median coverage")

	CoverageLay = go.Layout(title=''.join(Title),xaxis= dict(title='Mean interval position'),yaxis=dict(title='Mean interval coverage'),barmode='base-bar')

	CovPlot = [CoveragePlot,MeanLine]
	Cov = go.Figure(data=CovPlot, layout=CoverageLay)
    
	plot(Cov,filename=''.join(FileName),auto_open=False)


