#!/usr/bin/python

import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly import tools


File = "/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/PyInfo.txt"
Names = ["File_Name", "Number_of_Sequences", "Minimap_Time", "Mean_Length", "Median_Length"]


Table = pd.read_table(File, delimiter='\t',names=Names)
PyTable = Table.sort_values('Number_of_Sequences')


MinimapTime = go.Scatter(
    x=PyTable['Number_of_Sequences'],
    y=PyTable['Minimap_Time'],
    mode = 'lines+markers',
    text=PyTable['File_Name'],
    fill='tozeroy'
)

MeanLength = go.Scatter(	
    x=PyTable['Number_of_Sequences'],
    y=PyTable['Mean_Length'],
    mode = 'lines+markers',
    text=PyTable['File_Name'],
    fill='tozeroy'
)

MedianLength = go.Scatter(
    x=PyTable['Number_of_Sequences'],
    y=PyTable['Median_Length'],
    mode = 'lines+markers',
    text=PyTable['File_Name'],
    fill='tozeroy'
)

Satistics = tools.make_subplots(rows=3, cols=1,subplot_titles=('Minimap2 alignment time for different number of reads','Mean length of reads', 'Median length of reads'))
Satistics.append_trace(MinimapTime, 1, 1)
Satistics.append_trace(MeanLength, 2, 1)
Satistics.append_trace(MedianLength, 3, 1)
Satistics['layout']['xaxis1'].update(title='Number of reads')
Satistics['layout']['xaxis2'].update(title='Number of reads')
Satistics['layout']['xaxis3'].update(title='Number of reads')
Satistics['layout']['yaxis1'].update(title='Time (minutes)')
Satistics['layout']['yaxis2'].update(title='Mean length')
Satistics['layout']['yaxis3'].update(title='Median length')
Satistics['layout'].update(height=1000, width=1200, title='Statistics',showlegend=False)
plot(Satistics, filename="Statistics.html")

