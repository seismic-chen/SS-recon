#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 09:33:52 2021

@author: admin
"""

import os
import glob
import shutil
from obspy import read
from numpy import deg2rad, rad2deg, cos, sin, arctan2, sqrt

def gc_midpoint(lat1, lon1, lat2, lon2):
    lat1r = deg2rad(lat1)
    lat2r = deg2rad(lat2)
    lon1r = deg2rad(lon1)
    lon2r = deg2rad(lon2)
    
    bx = cos(lat2r) * cos(lon2r-lon1r)
    by = cos(lat2r) * sin(lon2r-lon1r)
    
    bplatr = arctan2( sin(lat1r)+sin(lat2r), sqrt((cos(lat1r)+bx)**2+by**2 ))
    bplonr = lon1r + arctan2(by, cos(lat1r)+bx)
    
    bplat = rad2deg(bplatr)
    bplon = rad2deg(bplonr)
    
    return bplat, bplon


trlist = glob.glob("./database/*/*.T")
outdir = '/media/admin/Data/western_pacific/data'
for f in trlist:
    dirname = f.split('/')[2]
    st = read(f)
    tr = st[0]
    
    evla = tr.stats.sac.evla
    evlo = tr.stats.sac.evlo
    stla = tr.stats.sac.stla
    stlo = tr.stats.sac.stlo
    
    bpla, bplo = gc_midpoint(evla, evlo, stla, stlo)
    
    if 20 <= bpla <= 60 and 110 <= bplo <= 160:
        destination = os.path.join(outdir, dirname)
        if not os.path.exists(destination):
            os.makedirs(destination)
        shutil.copy(f, destination)
        

## output bounce points
trlist = glob.glob("./data/*/*.T")
n = len(trlist)
bp = np.zeros((n,2))
k = 0
for f in trlist:
    
    st = read(f)
    tr = st[0]
    
    evla = tr.stats.sac.evla
    evlo = tr.stats.sac.evlo
    stla = tr.stats.sac.stla
    stlo = tr.stats.sac.stlo
    
    bpla, bplo = gc_midpoint(evla, evlo, stla, stlo)
    bp[k,:] = bpla, bplo
    k += 1
 
with open("Bounce_Point.txt", "w") as txt_file:
    for line in bp:
        txt_file.write("%10.3f %10.3f\n" % (line[0],line[1]))
txt_file.close()