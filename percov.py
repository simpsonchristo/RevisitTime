#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Constellation Revisit Time"""
import helperfunctions as hf
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

def percentcoverage(latP, inc, hsat, elconstraint):
    '''PERCENTCOVERAGE Analytical percent coverage and No. of coverage regions
    '''
    Rphi = hf.geodeticRadius(latP)
    Slam = Rphi/(Rphi+hsat)
    ups = math.asin(math.cos(elconstraint)*Slam)
    #max earth central angle
    lammax = math.pi/2 - (ups + elconstraint)
    #angular sensor range
    inc = inc if inc<math.pi/2 else math.pi - inc
    Ctheta1 = (-math.sin(lammax) + (math.cos(inc)*math.sin(latP)))/(math.sin(inc)*math.cos(latP))
    Ctheta2 = (math.sin(lammax) + (math.cos(inc)*math.sin(latP)))/(math.sin(inc)*math.cos(latP))
    if (latP>(lammax + inc)):
        return [0, 0]
    elif(latP<(lammax + inc) and latP>(lammax - inc)):
        theta1 = math.acos(Ctheta1)
        return [1, theta1/math.pi]
    elif(latP<(lammax - inc)):
        theta1 = math.acos(Ctheta1)
        theta2 = math.acos(Ctheta2)
        return [2, (theta1 - theta2)/math.pi]

def pdensity(latP, inc):
    return (math.cos(latP)/(math.pi*math.sqrt(math.sin(inc)**2 - math.sin(latP)**2)))