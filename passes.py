#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Helper Functions"""
import math
import helperfunctions as hf
import warnings
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
def tview(oe, phi, el):
    #time point p in view
    #unpack orbital elements
    a, ecc, inc, RAAN, AoP = oe
    #drift
    Pn = hf.nodalperiod(oe, phi)
    #position
    p   = hf.semilatus(ecc, a)
    Rlat= hf.geodeticRadius(phi)
    Rs  = hf.satgeodeticradius(p, ecc, phi)
    #lambda max
    CeA = math.pi/2 - (math.asin(math.cos(el)*(Rlat/(Rs)))) + el 
    #off-nadir
    eta = math.asin(math.sin(CeA)*math.sin((math.pi - inc)%(math.pi/2)))
    return Pn*math.acos(math.cos(CeA)*math.cos(eta))/math.pi

#Given target latitude, height above sphere, eccentricity, and grazing angle return longitude coverage
def loncover(lat, oe, el, descend):
    a, ecc, inc, RAAN, AoP = oe 
    if descend==0:
        Ta  = hf.trueanom(lat, AoP, inc)#rad, true anomaly
    else:
        Ta  = hf.trueanom(lat, AoP, inc)#rad, true anomaly
        u = (Ta + AoP) % 2*math.pi
        Ta  = (math.pi - (u + AoP)) % (2*math.pi)
    
    p   = hf.semilatus(ecc, a)
    Rlat= hf.geodeticRadius(lat)
    Rs  = hf.satgeodeticradius(p, ecc, Ta)
    Hsat= Rs - Rlat
    # 30 x (y x 5,400)
    if Hsat<0:
        warnings.warn('Hsat is below oblate spheriod')
    
    CeA = math.pi/2 - (math.asin(math.cos(el)*(Rlat/(Rs)))) + el 
    # 30 x [(y x 5,400) OR len(h) x 30 OR 0]
    lon = math.atan(math.tan(CeA)/math.sin((math.pi/2)-lat))
    # 30 x [(y x 5,400) OR len(h) x 30 OR 0]
    return Ta, CeA, lon

def lon1sat(oe, phi, el, lim):
    #unpack orbital elements
    a, ecc, inc, RAAN, AoP = oe
    #check longitudinal coverage
    TA, theta, lon = loncover(phi, oe, el, 0)
    TA2, theta2, lon2 = loncover(phi, oe, el, 1)
    #drift
    Pn = hf.nodalperiod(oe, phi)
    dRAAN = hf.dRAAN(oe, phi)
    dlon = Pn*(dRAAN - We)
    
    
    #Nodal Passes
    NoPass = [i for i in range(math.ceil(lim*3600*24/Pn))]
    if(TA<TA2):
        #initial ascend
        AoPTA = AoP + TA
        lon0 = math.atan2((math.cos(AoPTA)*math.sin(RAAN)) 
                          + (math.sin(AoPTA)*math.cos(RAAN)*math.cos(inc)),
                          (math.cos(AoPTA)*math.cos(RAAN)) 
                          - (math.sin(AoPTA)*math.sin(RAAN)*math.cos(inc)))
        lon0 += ((AoPTA % 2*math.pi)/(2*math.pi))*dlon
        #nodal passes
        tpass0 = Pn*(AoPTA % (math.pi*2))/(math.pi*2)
        tpasssingle = [(Npass * Pn) + tpass0 for Npass in NoPass]
        dRAAN = [(Npass * dRAAN) for Npass in NoPass]
        #given the pass time for each node, determine longitude pass at latitude
        #single satellite
        lonsingle = [(lon0 + (dlon*Npass))% (2*math.pi) for Npass in NoPass]
        theta = theta
    if(TA>TA2):
        #initial descend
        AoPTA = AoP + TA2
        lon0 = math.atan2((math.cos(AoPTA)*math.sin(RAAN)) 
                          + (math.sin(AoPTA)*math.cos(RAAN)*math.cos(inc)),
                          (math.cos(AoPTA)*math.cos(RAAN)) 
                          - (math.sin(AoPTA)*math.sin(RAAN)*math.cos(inc)))
        lon0 += ((AoPTA % 2*math.pi)/(2*math.pi))*dlon
        #nodal passes
        tpass0 = Pn*(AoPTA % (math.pi*2))/(math.pi*2)
        tpasssingle = [(Npass * Pn) + tpass0 for Npass in NoPass]
        dRAAN = [(Npass * dRAAN) for Npass in NoPass]
        #given the pass time for each node, determine longitude pass at latitude
        #single satellite
        lonsingle = [(lon0 + (dlon*Npass))% (2*math.pi) for Npass in NoPass]
        theta = theta2
    return lonsingle, tpasssingle, theta, dRAAN

def satsInPlan(oe, phi, lon, tpass, walker):
    #passes by satellites in plane
    #unpack orbital elements
    a, ecc, inc, RAAN, AoP = oe
    #unpack walker
    t, p, f = walker #no sats, no planes, rel spacing
    #drift
    Pn = hf.nodalperiod(oe, phi)
    dRAAN = hf.dRAAN(oe, phi)
    dlon = Pn*(dRAAN - We)
    #no. of satellites per plane
    s = t/p
    
    lon_s = [[(lj + ((l/s)*dlon))% (2*math.pi) for lj in lon] for l in range(int(s))]
    tp_s = [[tj + (Pn - ((l/s)*Pn)) for tj in tpass] for l in range(int(s))]
    return lon_s, tp_s

def planPass(oe, phi, lon, tpass, walker):
    #passes of reference satellite in each plane
    #unpack orbital elements
    a, ecc, inc, RAAN, AoP = oe
    #unpack walker
    t, p, f = walker #no sats, no planes, rel spacing
    #drift
    if type(inc)==list:
        Pn = []
        for i in inc:
            Pn.append(hf.nodalperiod([a, ecc, i, RAAN, AoP], phi))
        lon_p = [[(lj[m] + math.pi*2*m*((1/p) + (f/t)))% (2*math.pi) for lj in lon] for m in range(p)]
        tp_p = [[tj[m] + Pn[m]*math.pi*2*m*f/t for tj in tpass] for m in range(p)]
    else:
        Pn = hf.nodalperiod(oe, phi)
        lon_p = [[(lj + math.pi*2*m*((1/p) + (f/t)))% (2*math.pi) for lj in lon] for m in range(p)]
        tp_p = [[tj + Pn*math.pi*2*m*f/t for tj in tpass] for m in range(p)]
    return lon_p, tp_p

    
    
    