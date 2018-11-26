#!/usr/bin/python

############################################ Plot Minimap2 Time and Accuracy Performances in Alignment ############################################


### Plot Minimap2, Ngmlr and BWA time performances in alignment ###


import os
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly import tools
import pysam


def TimeLengthPlotter(FilesList):
 
    Header1 = ["File_Name", "Number_of_Sequences", "Minimap2_Time", "Mean_Length", "Min_Length", "Max_Length"]
    Header2 = ["File_Name", "Number_of_Sequences", "Ngmlr_Time", "Mean_Length", "Min_Length", "Max_Length"]
    Header3 = ["File_Name", "Number_of_Sequences", "BWA_Time", "Mean_Length", "Min_Length", "Max_Length"]


    Table1 = pd.read_table(FilesList[0], delimiter='\t',names=Header1)
    PyTable1 = Table1.sort_values('Number_of_Sequences')
    Time1 = go.Scatter(x=PyTable1['Number_of_Sequences'],y=PyTable1['Minimap2_Time'],mode = 'lines+markers',text=PyTable1['File_Name'],fill='tozeroy',name='Minimap2')
    Table2 = pd.read_table(FilesList[1], delimiter='\t',names=Header2)
    PyTable2 = Table2.sort_values('Number_of_Sequences')
    Time2 = go.Scatter(x=PyTable2['Number_of_Sequences'],y=PyTable2['Ngmlr_Time'],mode = 'lines+markers',text=PyTable2['File_Name'],fill='tozeroy',name='Ngmlr')
    Table3 = pd.read_table(FilesList[2], delimiter='\t',names=Header3)
    PyTable3 = Table3.sort_values('Number_of_Sequences')
    Time3 = go.Scatter(x=PyTable3['Number_of_Sequences'],y=PyTable3['BWA_Time'],mode = 'lines+markers',text=PyTable3['File_Name'],fill='tozeroy',name='BWA')
    TimeLay = go.Layout(title= 'Minimap2, Ngmlr and BWA alignment time', xaxis= dict(title='# Reads'),yaxis=dict(title='Time (minutes)', side='left'))
    Time=[Time1,Time2,Time3]
    TimePlot=go.Figure(data=Time, layout=TimeLay)

    plot(TimePlot,filename="TimePlot.html",auto_open=False)

    MeanLength = go.Scatter(x=PyTable1['Number_of_Sequences'],y=PyTable1['Mean_Length'],mode='lines+markers',text=PyTable1['File_Name'],fill='tozeroy',name='Mean')
    MinLength = go.Scatter(x=PyTable1['Number_of_Sequences'],y=PyTable1['Min_Length'],mode='lines+markers',text=PyTable1['File_Name'],fill='tozeroy',name='Min')
    MaxLength = go.Scatter(x=PyTable1['Number_of_Sequences'],y=PyTable1['Max_Length'],mode='lines+markers',text=PyTable1['File_Name'],fill='tozeroy',name='Max')
    LengthPlot = [MeanLength, MinLength, MaxLength]
    LengthPlotLayout = go.Layout(title= 'Sequences min/mean/max length',xaxis= dict(title='# Reads'),yaxis=dict(title='Length'))
    LengthFigure = go.Figure(data=LengthPlot, layout=LengthPlotLayout)
    
    plot(LengthFigure,filename="LengthPlot.html",auto_open=False)


## Executing codes ... ##

Files = ["/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2_Time_Performances/Minimap2_PyInfo.txt","/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Ngmlr_Time_Performances/Ngmlr_PyInfo.txt","/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/BWA_Time_Performances/BWA_PyInfo.txt"]

TimeLengthPlotter(Files)




### Plot Minimap2, Ngmlr and BWA accuracy performances in alignment ###


import os
import pysam
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly import tools

def Bam_and_ID_List_Generator(ListOfDirs):

    Bam = []
    ID = []

    for Dirs in ListOfDirs:

        BamInDir = []
        IDInDir = []

        for file in os.listdir(Dirs):

            if file.endswith(".srt.bam"):

                BamInDir.append(os.path.join(Dirs, file))
                IDLab = os.path.basename(os.path.join(Dirs, file))
                IDInDir.append(IDLab)

        Bam.append(BamInDir)
        ID.append(IDInDir)

    return Bam,ID

def BAM_Analyzer(BamList):

    TMap = []
    TUnmap = []
    RMap = []
    RUnmap = []


    for i in range(len(BamList)):

        TotalMapped = []
        TotalUnmapped = []
        RegionMapped = []
        RegionUnmapped = []

        for Bam in BamList[i]:

            BamFile = pysam.AlignmentFile(Bam, "rb")
            TotalMapped.append(int(BamFile.mapped))
            TotalUnmapped.append(int(BamFile.unmapped))
            MappedInRegion = 0
            UnMappedInRegion = 0

            for reads in BamFile.fetch("chr1",1000001,40000000):

                if not reads.is_unmapped:

                    MappedInRegion+=1
                
                else:

                    UnMappedInRegion+=1

            RegionMapped.append(MappedInRegion)
            RegionUnmapped.append(UnMappedInRegion)

        TMap.append(TotalMapped)
        TUnmap.append(TotalUnmapped)
        RMap.append(RegionMapped)
        RUnmap.append(RegionUnmapped)


    return TMap,TUnmap,RMap,RUnmap

def Results_Plotter(TMap,TUnmap,RMap,RUnmap,ID,Simul_Length,Simul_Accuracy):

    Methods1 = ["Minimap2_Accuracy_Plots","Ngmlr_Accuracy_Plots","BWA_Accuracy_Plots"]
    Methods2 = ["Minimap2_Accuracy_Table","Ngmlr_Accuracy_Table","BWA_Accuracy_Table"]

    for i in range(len(ID)):

        Id = ID[i]
        TotalMapped = TMap[i]
        TotalUnmapped = TUnmap[i]
        RegionMapped = RMap[i]
        RegionUnmapped = RUnmap[i]

        InfoTable = pd.DataFrame({"ID":[i.split(".")[0] for i in Id], "# Mapped Reads in Bam": TotalMapped, "# Unmapped Reads in Bam":TotalUnmapped,"# Mapped Reads in Region": RegionMapped, "# Unmapped Reads in Region":RegionUnmapped})

        MapFraction = InfoTable["# Mapped Reads in Bam"]/(InfoTable["# Mapped Reads in Bam"]+InfoTable["# Unmapped Reads in Bam"])
        RightMapped = InfoTable["# Mapped Reads in Region"]/InfoTable["# Mapped Reads in Bam"]


        TotMapped = go.Bar(x = InfoTable["ID"], y = InfoTable["# Mapped Reads in Bam"], name = "# Mapped Reads in Bam")
        TotUnMapped = go.Bar(x = InfoTable["ID"], y = InfoTable["# Unmapped Reads in Bam"], name = "# Unmapped Reads in Bam")
        MappedInRegion = go.Bar(x = InfoTable["ID"], y = InfoTable["# Mapped Reads in Region"], name = "# Mapped Reads in Region")
        UnMappedInRegion = go.Bar(x = InfoTable["ID"], y = InfoTable["# Unmapped Reads in Region"], name = "# Unmapped Reads in Region")

        
        TableInfoNew = go.Table(header=dict(values=["ID","Mean Length","Mean Accuracy", "# Mapped Reads in Bam", "# Unmapped Reads in Bam", "# Mapped Reads in Region", "# UnMapped Reads in Region", "# Mapped Reads in Bam / #Total", "# Mapped Reads in Region / # Mapped Reads in Bam"],line = dict(color="#7D7F80"),fill = dict(color="#a1c3d1"),align = ["left"] * 5), cells=dict(values=[InfoTable["ID"],Simul_Length,Simul_Accuracy,InfoTable["# Mapped Reads in Bam"],InfoTable["# Unmapped Reads in Bam"],InfoTable["# Mapped Reads in Region"],InfoTable["# Unmapped Reads in Region"],["{0:.4f}".format(j) for j in MapFraction],["{0:.4f}".format(j) for j in RightMapped]],line = dict(color="#7D7F80"),fill = dict(color="#EDFAFF"),align = ["left"] * 5))
        TableLay = go.Layout(title = Methods2[i])

        Fig1 = tools.make_subplots(rows=2, cols=2, specs=[[{}, {}], [{}, {}]],subplot_titles=("# Mapped Reads in Bam","# Unmapped Reads in Bam","# Mapped Reads in Region","# Unmapped Reads in Region"))
        Fig2 = go.Figure(data = [TableInfoNew], layout=TableLay)


        Fig1.append_trace(TotMapped, 1, 1)
        Fig1.append_trace(TotUnMapped, 1, 2)
        Fig1.append_trace(MappedInRegion, 2, 1)
        Fig1.append_trace(UnMappedInRegion, 2, 2)
        Fig1['layout'].update(title=Methods1[i])

        
        plot(Fig1, filename=Methods1[i]+".html", auto_open=False)
        plot(Fig2, filename=Methods2[i]+".html", auto_open=False )



## Executing codes ... ##

MyDirs = ["/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2_Accuracy_Performances","/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Ngmlr_Accuracy_Performances","/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/BWA_Accuracy_Performances"]
Bam,Id = Bam_and_ID_List_Generator(MyDirs)
TMap,TUnmap,RMap,RUnmap = BAM_Analyzer(Bam)
SimulatedReadsLength = [3000,8000,13000,3000,8000,13000,3000,8000,13000,10000,100000] #length (bp) used to create synthetic fastq with pbsim
SimulatedReadsAccuracy = [0.85,0.85,0.85,0.90,0.90,0.90,0.95,0.95,0.95,0.90,0.98] #accuracy used to create synthetic fastq with pbsim
Results_Plotter(TMap,TUnmap,RMap,RUnmap,Id,SimulatedReadsLength,SimulatedReadsAccuracy)





