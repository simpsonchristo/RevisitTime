#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Constellation Revisit Time"""
import helperfunctions as hf
import math
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

#Given target latitude, height above sphere, eccentricity, and grazing angle return longitude coverage
def loncover(lat, oe, el, descend):
    a, ecc, inc, RAAN, AoP = oe 
    if descend==0:
        Ta  = hf.trueanom(lat, AoP, inc)#rad, true anomaly
    else:
        Ta  = hf.trueanom(lat, AoP, inc)#rad, true anomaly
        u = (Ta + AoP) % 2*math.pi
        Ta  = (math.pi - (u + AoP)) % (2*math.pi)
    # 5,400
    p   = hf.semilatus(ecc, a)
    # y
    Rlat= hf.geodeticRadius(lat)
    # 30
    Rs  = [[p[i]/(1 + (ecc[i]*math.cos(nu[j]))) for i in range(len(p))] for nu in Ta for j in range(len(nu))]
    # y x 5,400
    Hsat= [[Rsat - Rphi for sublist in Rs for Rsat in sublist] for Rphi in Rlat]
    # 30 x (y x 5,400)
    if min([Hcalc for sublist in Hsat for Hcalc in sublist])<0:
        warnings.warn('Hsat is below oblate spheriod')
    
    CeA = [math.pi/2 - (math.asin(math.cos(el)*(Rlat[i]/(Rs[i])))) + el for i in range(len(Rlat))]
    # 30 x [(y x 5,400) OR len(h) x 30 OR 0]
    lon = [math.atan(math.tan(CeA[i])/math.sin((math.pi/2)-lat[i])) for i in range(len(lat))]
    # 30 x [(y x 5,400) OR len(h) x 30 OR 0]
    return Ta, CeA, lon

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
    #check lat against inc
    if phi > any(inc):
        warnings.warn('Inclination may be too small.')
    
    #check longitudinal coverage
    TA, theta, lon = loncover(phi, oe, El, 0)
    TA2, theta2, lon2 = loncover(phi, oe, El, 1)
    #drift
    Pn = hf.nodalperiod(oe, phi)
    dRAAN = hf.dRAAN(oe, phi)
    dlon = Pn*(dRAAN - We)
    
    #Nodal Passes
    NoPass = [i for i in range(Limit[0]*3600*24/Pn)]
    #initial ascend
    AoPTA = AoP + TA
    lon0 = math.atan2((math.cos(AoPTA)*math.sin(RAAN)) 
                      + (math.sin(AoPTA)*math.cos(RAAN)*math.cos(inc)),
                      (math.cos(AoPTA)*math.cos(RAAN)) 
                      - (math.sin(AoPTA)*math.sin(RAAN)*math.cos(inc)))
    lon0 += ((AoPTA % 2*math.pi)/(2*math.pi))*dlon
    #initial descending
    AoPTA2 = AoP + TA2
    lon1 = math.atan2((math.cos(AoPTA2)*math.sin(RAAN)) 
                      + (math.sin(AoPTA2)*math.cos(RAAN)*math.cos(inc)),
                      (math.cos(AoPTA2)*math.cos(RAAN)) 
                      - (math.sin(AoPTA2)*math.sin(RAAN)*math.cos(inc)))
    lon1 += ((AoPTA2 % 2*math.pi)/(2*math.pi))*dlon
    #nodal passes
    tpass0 = Pn*(AoPTA % (math.pi*2))/(math.pi*2)
    tpassall = [(Npass * Pn) + tpass0 for Npass in NoPass]
    #given the pass time for each node, determine longitude pass at latitude
    #single satellite
    lonsingle = [(lon0 + (dlon*Npass))% (2*math.pi) for Npass in NoPass]
    
    if NoPlan>1:
        dphase = RelSpace*math.pi*2/NoSats
        lonp = []
        for m in range(NoPlan):
            lonp.append([lonsingle + math.pi*2*m*((1/NoPlan)+(dphase/NoSats))])
    if NoSatPerPlane > 1:
        lons = []
        for l in range(NoSatPerPlane):
            lons.append([lonsingle + ((l/NoSatPerPlane)*dlon)])
    