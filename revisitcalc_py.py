#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Constellation Revisit Time"""
import math
from passes import lon1sat, satsInPlan, planPass, tview
"""Python 3.7
   Simpson Aerospace (c) 2021
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
#constants
Re = 6378.15 #km
Rp = 6.3567523e+6/1000 #km
We = 7.292115856e-5 #rad/sec
mu = 3.986012e+5 #km3/sec2
ecliptic = (23*math.pi)/180 #rad
J2 = -1.082626683e-03



def revisitcalc(oe, phi, lam, El, NoSats, NoPlan, RelSpace, Limit):
    '''
    REVISITCALC, Calculate Mean Revisit Time (MRT) and Maximum Gap (Mgap) of
    target point given latitude, phi, and Longitude, lam. 

    Parameters
    ----------
    oe : list
        Orbital elements [apo (km), ecc (-), inc (rad), RAAN (rad), AoP (rad)].
    phi : float
        Target latitude (rad).
    lam : float
        Target longitude (rad).
    El : float
        Elevation constraint on viewing angle (rad).
    NoSats : list
        Number of satellites.
    NoPlan : list
        Number of equispaced orbital planes.
    RelSpace : list
        Relative spacing of satellites in planes.
    Limit : int
        Time in days for limit of computation.

    Returns
    -------
    None.

    '''
    #unpack orbital elements
    a, ecc, inc, RAAN, AoP = oe
    #no. of satellites per plane
    NoSatPerPlane = NoSats/NoPlan
    
    #reference satellite: longitude and time of passes at latitude phi
    if type(inc)==list:
        lon = []
        tpass = []
        theta = []
        dRAAN = []
        for i in inc:
            x, y, z, node = lon1sat([a, ecc, i, RAAN, AoP], phi, El, Limit)
            lon.append(x)
            tpass.append(y)
            theta.append(z) 
            dRAAN.append(node)
    else:
        lon, tpass, theta, dRAAN = lon1sat(oe, phi, El, Limit)
    lon_all = []
    tp_all = []
    #planar passes
    if NoPlan>1:
        lon_p, tp_p= planPass(oe, phi, lon, tpass, [NoSats, NoPlan, RelSpace])
        lon = lon_p
        tpass = tp_p
    #satellites in plane passes
    if NoSatPerPlane>1:
        if(type(lon[0])==list):
            if(type(inc)==list):
                for i in range(len(lon)):
                    lon_s, tp_s= satsInPlan([a, ecc, inc[i], RAAN, AoP], phi, lon[i], tpass[i], [NoSats, NoPlan, RelSpace])
                    lon_all.append(lon_s)
                    tp_all.append(tp_s)
            else:
                for i in range(len(lon)):
                    lon_s, tp_s= satsInPlan([a, ecc, inc, RAAN, AoP], phi, lon[i], tpass[i], [NoSats, NoPlan, RelSpace])
                    lon_all.append(lon_s)
                    tp_all.append(tp_s)
        else:
            lon_s, tp_s= satsInPlan(oe, phi, lon, tpass, [NoSats, NoPlan, RelSpace])
            lon_all.append(lon_s)
            tp_all.append(tp_s)
    #check if point in view of pass
    if type(inc)==list:
        inview = lon_all
        for k in range(len(lon_all)):
            for j in range(len(lon_all[k])):
                for i in range(len(lon_all[k][j])):
                    inview[k][j][i] = abs(lon_all[k][j][i] - lam) <= theta[k]     
    else:
        inview = [[[abs(sat[i] - lam) <= theta for i in range(len(sat))] for sat in plane] for plane in lon_all]     
    
    #return time and longitude of passes
    lonpass = []
    whenview = []  
    view_time = []
    mrt = [tp_all[j][i][k] for j in range(len(inview)) for i in range(len(inview[j])) for k in range(1,len(inview[j][i]))]
    for j in range(len(inview)):
        for i in range(len(inview[j])):
            for k in range(len(inview[j][i])):
                if(inview[j][i][k]):
                    whenview.append(tp_all[j][i][k])
                    lonpass.append(lon_all[j][i][k])
                    if type(inc)==list:
                        oem = [a, ecc, inc[j], RAAN - dRAAN[j][k], AoP]
                        view_time.append(tview(oem, phi, El))
                    else:
                        oem = [a, ecc, inc, RAAN - dRAAN[k], AoP]
                        view_time.append(tview(oem, phi, El))
    
    t_start = sorted(whenview)
    t_end   = [t_start[i] + view_time[whenview.index(t_start[i])] for i in range(len(view_time))]
    t_end.insert(0,0)
    mrt = [t_end[i] - t_end[i-1] if (t_end[i] - t_end[i-1])>0 else 0 for i in range(1,len(t_end))]
    if(len(view_time)!=0):
        mrt  = sum(mrt)/len(view_time)
        gmti = [view_time[whenview.index(t_start[i])] for i in range(len(view_time))]
        gmti = sum(gmti)/len(view_time)
    else:
        mrt  = float('nan')
        gmti = float('nan')
    return mrt, gmti