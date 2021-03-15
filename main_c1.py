#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Constellation Revisit Time"""
import helperfunctions as hf
from revisitcalc_py import revisitcalc
from percov import percentcoverage
import math
"""Python 3.7
   Simpson Aerospace (c) 2021
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
#CASE 1 Shandong 
#sat range
Hsat = [i for i in range(250,500,50)]#km
ecc  = 0 #[]
#point parameters
lon  = 360-(175 + (6/60) + (58/3600))
lat  = (21 + (42/60) + (30/3600))
El   = 10

#Walker delta constellation
sats = 8 #No of Sats
plan = [p for p in range(1,sats+1)] #planes
spac = [f for p in plan for f in range(p)] #spacing
lim  = 1  #Loop control [days]

#determine inc to evaluate
test = [math.radians(i) for i in range(1,179,1)]
# inc  = [90, 88, 96] 
inc  = [] #inclination sample of largest percentage coverage
incs = [] #inclinations providing some coverage
cper = [] #percent coverage 
for h in Hsat:
    for i in test:
        a, b = percentcoverage(math.radians(lat), i, h, math.radians(El))
        if (a!=0):
            incs.append(i)
            cper.append(b)
for i in range(len(plan)):
    inc.append(incs.pop(cper.index(max(cper))))
    print(cper.pop(cper.index(max(cper))))

OE  = []
mrt = []
wal = []
#height, planes & inclination, spacing
for h in Hsat:
    #find inclination for minimum time
    inc.append(hf.minInc(4, math.radians(60), lat, El))
    for p in plan:
        for f in spac:
            if(p>1):
                OE.append([hf.apo(h,ecc), ecc, inc[:p], math.radians(lon), 0])
            else:
                #5 orbital elements: A, E, I, RAAN, TA
                OE.append([hf.apo(h,ecc), ecc, inc[-1], math.radians(lon), 0])
            a, b = revisitcalc(OE[-1],
                               math.radians(lat),
                               math.radians(lon),
                               math.radians(El),
                               sats,
                               p,
                               f,
                               lim)
            mrt.append(a)
            wal.append([4,p,f])
with open('case1start.txt','w') as f:
    for i in range(len(mrt)):
        f.write(f"{mrt[i]},{wal[i]},{OE[i]}\n")
    