#!/usr/bin/python

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
    x=PyTable['Minimap_Time'],
    y=PyTable['Number_of_Sequences'],
    mode = 'lines+markers',
    text=PyTable['File_Name'],
    fill='tozeroy'
)

MinimapTimeLayout = go.Layout(
    title= 'Alignment time',
    xaxis= dict(title='Time (minutes)'),
    yaxis=dict(title='Number of reads', side='right'),
)

MiniPlot = [MinimapTime]
MinimapFigure = go.Figure(data=MiniPlot, layout=MinimapTimeLayout)
plot(MinimapFigure,filename="Time.html" )


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
    xaxis= dict(title='Number of reads'),
    yaxis=dict(title='Length'),
)
LengthFigure = go.Figure(data=LengthPlot, layout=LengthPlotLayout)
plot(LengthFigure,filename="Length.html" )

