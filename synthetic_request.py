#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:52:27 2020

@author: admin
"""

import os, glob
import obspy
import numpy as np
from obspy import read
from obspy.clients.syngine import Client
client = Client()
from obspy.taup import TauPyModel
model = TauPyModel(model="prem")

with open('events.dat') as f:
    lines = f.readlines()
    
dt = 1.0
work_dir = './data/'   
for line in lines:
    o, evla, evlo, evdp, mag, gcmt = line.split()   
    evdir = glob.glob(work_dir+o+'*')[0]
    os.chdir(evdir)
    
    for sac in glob.glob("*.T"):
        T = read(sac)[0]
        net = T.stats.sac['knetwk']
        sta = T.stats.sac['kstnm']
        evdp = T.stats.sac['evdp'] / 1000.
        dist = T.stats.sac['gcarc']
        print(net,sta)
        st = client.get_waveforms(model="prem_i_2s",
                              network=net,
                              station=sta,
                              eventid="GCMT:"+gcmt,
                              components="ZRT",
                              starttime=0,
                              endtime=3500)
        st.filter("lowpass", freq=dt/2, corners=2, zerophase=True)
        st.interpolate(sampling_rate=dt)
        st.filter("bandpass", freqmin=0.013, freqmax=0.067, corners=2, zerophase=True)
        tr = st.select(component="T")[0]
        tr.detrend('linear')
        tr.detrend('demean')
        tr.normalize()
        
        # reverse polarity if necessary
        arrivals = model.get_travel_times(source_depth_in_km=evdp,
                                          distance_in_degree=dist,
                                          phase_list=["SS"])
        ss = arrivals[0].time
        ts1, ts2 = ss - 50, ss + 50        
        is1 = int(round( (ts1)/dt ))
        is2 = int(round( (ts2)/dt ))
        win_SS = tr.data[is1:is2]
        max_ind = np.argmax(abs(win_SS))
        if win_SS[max_ind] < 0:
            tr.data = tr.data * -1
        
        o = tr.stats.starttime
        tr.trim(o+ss-900, o+ss+400)
        # save to file
        tr.stats.sac = obspy.core.AttribDict()
        tr.stats.sac['stla'] = T.stats.sac['stla']
        tr.stats.sac['stlo'] = T.stats.sac['stlo']
        tr.stats.sac['stel'] = T.stats.sac['stel']
        tr.stats.sac['stdp'] = T.stats.sac['stdp']
        tr.stats.sac['evla'] = T.stats.sac['evla']
        tr.stats.sac['evlo'] = T.stats.sac['evlo']
        tr.stats.sac['evel'] = T.stats.sac['evel']
        tr.stats.sac['evdp'] = T.stats.sac['evdp']
        tr.stats.sac['mag'] = T.stats.sac['mag']
        
        syn_dir = evdir.replace('./data', './synthetics')
        if not os.path.exists(syn_dir):
            os.makedirs(syn_dir)
        tr.write(os.path.join(syn_dir,sac), format="SAC")
        
    os.chdir("../..")
