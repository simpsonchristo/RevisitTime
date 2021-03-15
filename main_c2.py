#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Constellation Revisit Time"""
import helperfunctions as hf
from revisitcalc_py import revisitcalc
from percov import pdensity
import math
"""Python 3.7
   Simpson Aerospace (c) 2021
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
#CASE 2 Mischief Reef
#sat range
Hsat = [i for i in range(250,500,50)]#km
ecc  = 0 #[]
#point parameters
lon  = (115 + (32/60) + (0/3600))
lat  = (9 + (54/60) + (0/3600))
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
cper = [] #coverage density
for h in Hsat:
    for i in test:
        if(i>math.radians(lat) and i<math.radians(lat + 90)):
            cper.append(pdensity(math.radians(lat), i))
            incs.append(i)
for i in range(len(plan)):
    inc.append(incs.pop(cper.index(max(cper))))
    print(cper.pop(cper.index(max(cper))))

OE  = []
mrt = []
wal = []
gmti = []
#height, planes & inclination, spacing
for h in Hsat:
    #find inclination for minimum time
    inc.append(hf.minInc(ecc, hf.apo(h,ecc), 0, lat, El))
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
            gmti.append(b)
            wal.append([4,p,f])
with open('case2start.txt','w') as f:
    for i in range(len(mrt)):
        f.write(f"{mrt[i]},{gmti[i]},{wal[i]},{OE[i]}\n")
    