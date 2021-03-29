#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Helper Functions"""
import math
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

def apo(Hsat, ecc):
    if  (type(Hsat) == list and type(ecc) == list):
        return [[((1+((1-e)/(1+e)))/2)*(Re + alt) for alt in Hsat] for e in ecc]
    elif(type(Hsat) == list):
        return [((1+((1-ecc)/(1+ecc)))/2)*(Re + alt) for alt in Hsat]
    elif(type(ecc) == list):
        return [((1 + ((1 - e)/(1 + e)))/2)*(Re + Hsat) for e in ecc]
    else:
        return ((1+((1-ecc)/(1+ecc)))/2)*(Re + Hsat)
def trueanom(lat, CeA, aop, inc):
    if(inc>math.pi/2):
        inc = math.pi - inc
    if(inc>lat):
        return math.asin(math.sin(lat)/math.sin(inc)) - aop 
    else:
        ang = lat - CeA if(lat-CeA>0) else 0
        return math.asin(math.sin(ang)/math.sin(inc)) - aop 
    
def period(apoapsis):
    if(type(apoapsis) == list):
        return [2*math.pi*math.sqrt(a**3/mu) for a in apoapsis]
    else:
        return 2*math.pi*math.sqrt(apoapsis**3/mu)
def semilatus(ecc, apo):
    if (type(apo) == list and type(ecc) == list):
        return [apo[i] * (1 - ecc[i]**2) for i in range(len(ecc))]
    else:
        return apo * (1 - ecc**2)
def geodeticRadius(lat):
    if(type(lat)==list):
        return [math.sqrt(((math.cos(phi)*Re**2)**2 + (math.cos(phi)*Rp**2)**2)
                          / ((math.cos(phi)*Re)**2 + (math.cos(phi)*Rp)**2)) for phi in lat]
    else:
        return math.sqrt(((math.cos(lat)*Re**2)**2 + (math.cos(lat)*Rp**2)**2)
                         / ((math.cos(lat)*Re)**2 + (math.cos(lat)*Rp)**2))
def satgeodeticradius(p, ecc, nu):
    if(type(p)==list and type(nu)==list):
        return [[p[i]/(1 + (ecc[i]*math.cos(nu[j]))) for i in range(len(p))] for j in range(len(nu))]
    elif(type(p)==list):
        return [p[i]/(1 + (ecc[i]*math.cos(nu))) for i in range(len(p))]
    elif(type(nu)==list):
        return [p/(1 + (ecc*math.cos(nu[j]))) for j in range(len(nu))]
    else:
        return p/(1 + (ecc*math.cos(nu)))
#nodal period for J2 perturbed orbit
def nodalperiod(oe, lat):
    a, ecc, inc, RAAN, AoP = oe
    Pk = period(a)
    # y
    p  = semilatus(ecc, a)
    # y
    if(type(Pk)==list and type(inc)==list):
        return [[Pk[j]/(1 + ((3*J2*Re/(4*p[j]))*((math.sqrt(1 - ecc[j]**2)*(2 - 3*math.sin(i)**2)) + (4 - 5*math.sin(i)**2)))) for j in range(len(Pk))]for i in inc]
    elif(type(Pk)==list):
        return [Pk[j]/(1 + ((3*J2*Re/(4*p[j]))*((math.sqrt(1 - ecc[j]**2)*(2 - 3*math.sin(inc)**2)) + (4 - 5*math.sin(inc)**2)))) for j in range(len(Pk))]
    elif(type(inc)==list):
        return [Pk/(1 + ((3*J2*Re/(4*p))*((math.sqrt(1 - ecc**2)*(2 - 3*math.sin(i)**2)) + (4 - 5*math.sin(i)**2)))) for i in inc]
    else:
        return Pk/(1 + ((3*J2*Re/(4*p))*((math.sqrt(1 - ecc**2)*(2 - 3*math.sin(inc)**2)) + (4 - 5*math.sin(inc)**2))))
def dRAAN(oe, lat):
    a, ecc, inc, RAAN, AoP = oe
    p  = semilatus(ecc, a)
    n  = math.sqrt(mu/a**3)
    return ((-3/2)*n*J2*((Re/p)**2)*math.cos(inc)) + ((3/32)*n*(J2**2)*((Re/p)**4)*math.cos(inc)*(12 - 4*ecc**2 - (80 + 5*ecc**2)*math.sin(inc)**2))

def minInc(lon, lat_sat, lat_P, el):
    lat = lat_sat
    # p   = semilatus(ecc, a)
    # Rlat= geodeticRadius(lat)
    # Rs  = satgeodeticradius(p, ecc, 0)
    
    # CeA = (math.pi/2 - (math.asin(math.cos(el)*(Rlat/(Rs)))) + el)%(math.pi*2)
    CeA = math.acos(math.cos(math.pi/2 - lat)*math.cos(math.pi/2 - lat_P) 
                    + (math.sin(math.pi/2 - lat)*math.sin(math.pi/2 - lat_P)*math.cos(lon)))
    # lon = math.atan2(math.tan(CeA),math.sin((math.pi/2)-lat))
    if((lat_P - lat_sat)>0):
        return math.asin(math.sin(lat_P)*math.sin(CeA))
    else:
        return math.asin(math.sin(lat_P)*math.sin(CeA))%math.pi