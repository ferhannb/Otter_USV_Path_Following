#!/usr/bin/env python
# -*- coding: utf-8 -*-


import plotly.graph_objects as go
from plotly.subplots import make_subplots


class Graphs():

    def __init__(self,rows=4,cols=1):
        self.figure = go.Figure()
        self.rows = rows
        self.cols = cols 
        self.figure = make_subplots(rows=self.rows,cols=self.rows)


    def plotly_results(self,time,speed,commad_speed,):


    fig = make_subplots(rows=4, cols=1)
    # Add traces
    fig.add_trace(go.Scatter(x=time, y=speed,
                        mode='lines',
                        name='araç hiz'),row=1, col=1)
    fig.add_trace(go.Scatter(x=timeotter, y=command_speed,
                        mode='lines',
                        name='cmd hiz'),row=1, col=1)

    

    fig.add_trace(go.Scatter(x=timeotter, y=speed_error_list,
                        mode='lines',
                        name='speed error'),row=1, col=1)

    fig.add_trace(
        go.Scatter(x=timeotter, y=u_actual_list_sancak,name=' Sancak pervane'),row=2, col=1
        )
    
    fig.add_trace(
        go.Scatter(x=timeotter, y=u_actual_list_iskele,name='İskele Pervane'),
        row=2, col=1)

    fig.add_trace(
        go.Scatter(x=timeotter, y=controlsignal_sancak,name='control signali_sancak'),
        row=2, col=1)
    
    fig.add_trace(
        go.Scatter(x=timeotter, y=controlsignal_iskele,name='control signali_iskele'),
        row=2, col=1)

    

    fig.add_trace(
        go.Scatter(x=timeotter, y=heading_signal_list,name='heading referans signali'),
        row=3, col=1)
    
    fig.add_trace(
        go.Scatter(x=timeotter, y=heading_error_list,name='heading error'),
        row=3, col=1)

    fig.add_trace(
    go.Scatter(x=timeotter, y=heading_list,name='heading otter'),
    row=3, col=1)

    fig.add_trace(
        go.Scatter(x=x_list, y=y_list,name='x-y pervane'),row=4, col=1
        )
    fig.update_layout( title_text="Sinus referansi")
    # plotly.offline.plot(fig, filename="Kare1.5-2.2-2.html")

    fig.show()