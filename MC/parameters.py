#!/usr/bin/env python

## LIBRARIES

import sys
import numpy as np


##########################################################################################
##########################################################################################
# FUNCTIONS

# This function takes in input the parameters of the partilces and returns
# a dictionary where the parameters are stored
# Note that we also set the center of mass in the origin and the patches aligned
# along the z axis.
def set_IPP(params):

    a, sigma_c, sigma_p, delta_c, delta_p, K = params

    IPP = {}

    IPP['sigma_c'] = sigma_c
    IPP['a'] = a

    ### Quantities that are specific to the Overalp-Of-Spheres model
    IPP['sigma_p'] = sigma_p
    IPP['delta_c'] = delta_c
    IPP['delta_p'] = delta_p

    ### Quantities that are specific to the Exponential model
    IPP['K'] = K


    ### Center of mass:
    IPP['r'] = np.array([0, 0, 0])

    ### Patches:
    p1 = np.array([0, 0, a])
    p2 = np.array([0, 0, -a])

    IPP['p1'] = {}
    IPP['p2'] = {}

    IPP['p1']['r'] = p1
    IPP['p2']['r'] = p2

    return IPP



# This function takes in input two IPPs, the u vector and the two parameters a and sigma_c.
# It computes the matrix W, inverts it and return the vector epsilon
# To compute the matrix W, we set the particles in each of the three reference configurations
# And we compute the three matrix elements associated to the three pair of interactions site (center-center), (center-off center), (off center-off center)
def compute_e_vec(IPP1, IPP2, u, a, sigma_c, k, weight):

    W = np.zeros((3,3))                                                          # Weight matrix

    confs = ['EE', 'EP', 'PP']
    for i, conf in enumerate(confs):

        IPP1, IPP2 = reference_config(IPP1, IPP2, a, sigma_c, conf)             # Put in reference configuration

        if weight == 'os':
            w = overlap_vector(IPP1, IPP2)
        elif weight == 'exp':
            w = exponential_weights(IPP1, IPP2, k)

        W[i, :] = w

    e = e_from_u(W, u)                                                           # e_EE, e_EP, e_PP

    return e


# This function puts the particles in one of the reference configurations, to be given in input by conf
def reference_config(IPP1, IPP2, a, sigma_c, conf):

    # Particle 1
    IPP1['r'] = np.array([0, 0, 0])
    IPP1['p1']['r'] = np.array([0, 0, a])
    IPP1['p2']['r'] = np.array([0, 0, -a])

    # Particle 2
    IPP2['r'] = np.array([0, 2*sigma_c, 0])
    IPP2['p1']['r'] = IPP2['r'] + np.array([0, 0, a])
    IPP2['p2']['r'] = IPP2['r'] + np.array([0, 0, -a])

    # They are now in EE config, so:

    # (i) return EE
    if conf == 'EE':
        return IPP1, IPP2


    # Define the rotations
    T = np.pi/2
    v = np.array([1, 0, 0])

    # Rotate particle 2 so that EE -> EP
    if conf == 'EP':
        IPP2 = rotate_IPP_around_its_center(v, T, IPP2)
        return IPP1, IPP2

    # Rotate both particles so that EE -> PP
    if conf == 'PP':
        IPP1 = rotate_IPP_around_its_center(v, -T, IPP1)
        IPP2 = rotate_IPP_around_its_center(v, T, IPP2)
        return IPP1, IPP2

## Rotation matrix of a rotation by an angle theta around a versor v
def rotation_matrix(v, theta):

    ct = np.cos(theta)
    omct = 1 - ct
    st = np.sin(theta)
    omst = 1 - st

    R = np.array((
        [ct + v[0]*v[0]*omct, v[0]*v[1]*omct - v[2]*st, v[0]*v[2]*omct + v[1]*st],
        [v[0]*v[1]*omct + v[2]*st, ct + v[1]*v[1]*omct, v[1]*v[2]*omct - v[0]*st],
        [v[0]*v[2]*omct - v[1]*st, v[1]*v[2]*omct + v[0]*st, ct + v[2]*v[2]*omct] ))

    return R

## Generate the rotation matrix and rotate the patches of one IPP
def rotate_IPP_around_its_center(v, dtheta, IPP):

    R = rotation_matrix(v, dtheta)

    ## Get the absolute vector of the first patch, rotate it and reattach it to the COM
    dr1 = IPP['p1']['r'] - IPP['r']
    IPP['p1']['r'] = IPP['r'] + R.dot(dr1)

    ## Get the absolute vector of the second patch, rotate it and reattach it to the COM
    dr2 = IPP['p2']['r'] - IPP['r']
    IPP['p2']['r'] = IPP['r'] + R.dot(dr2)


    return IPP

'''
Overlap of spheres model
'''

## Compute the three matrix elements of W for a given configuration
def overlap_vector(IPP1, IPP2):

    R = IPP1['sigma_c'] + IPP1['delta_c']               # Radius of the interaction sphere
    r = IPP1['sigma_p']                                 # Radius of the hard sphere

    rab = IPP1['r'] - IPP2['r']
    rab = np.sqrt(np.dot(rab, rab))                     # Distance between the center of mass

    if rab < 2*IPP1['sigma_c']:                         # OVERLAP!
        return np.array([1e12, 1e12, 1e12])

    # EE
    O_EE = single_overlap_value(R, R, IPP1['r'], IPP2['r'])

    # EP
    O_EP = 0
    O_EP += single_overlap_value(R, r, IPP1['r'], IPP2['p1']['r'])    # Equator 1, patch 2 A
    O_EP += single_overlap_value(R, r, IPP1['r'], IPP2['p2']['r'])    # Equator 1, patch 2 B
    O_EP += single_overlap_value(r, R, IPP1['p1']['r'], IPP2['r'])    # Equator 2, patch 1 A
    O_EP += single_overlap_value(r, R, IPP1['p2']['r'], IPP2['r'])    # Equator 2, patch 1 B

    # PP
    O_PP = 0
    O_PP += single_overlap_value(r, r, IPP1['p1']['r'], IPP2['p1']['r'])  # patch 1 A, patch 2 A
    O_PP += single_overlap_value(r, r, IPP1['p1']['r'], IPP2['p2']['r'])  # patch 1 A, patch 2 B
    O_PP += single_overlap_value(r, r, IPP1['p2']['r'], IPP2['p1']['r'])  # patch 1 B, patch 2 A
    O_PP += single_overlap_value(r, r, IPP1['p2']['r'], IPP2['p2']['r'])  # patch 1 B, patch 2 B


    # Normalization
    O = (4/3) * np.pi * (IPP1['sigma_c']**3)

    return np.array([O_EE, O_EP, O_PP]) / O


## This fuction computes the overlap volume of two spheres with radii Ra and Rb
## whose centers are in coordinates specified by the vectors ra and rb
def single_overlap_value(Ra, Rb, ra, rb):

    rab = ra - rb
    rab = np.sqrt(np.dot(rab, rab))

    Rm = min(Ra, Rb)

    rmax = Ra + Rb
    rmin = np.fabs(Ra - Rb)

    if rab >= rmax:
        return 0
    elif rab <= rmin:
        return (4/3)*np.pi*(Rm**3)
    elif rmin <= rab and rab <= rmax:

        f1 = 2*Ra + (Ra**2 - Rb**2 + rab**2) / (2*rab)
        f2 = Ra - (Ra**2 - Rb**2 + rab**2) / (2*rab)

        f3 = 2*Rb - (Ra**2 - Rb**2 - rab**2) / (2*rab)
        f4 = Rb + (Ra**2 - Rb**2 - rab**2) / (2*rab)

        f = f1*f2*f2 + f3*f4*f4

        return f * np.pi/3

    else:
        print('Something wierd is going on')
        return np.nan


'''
Exponential model
'''

def exponential_weights(IPP1, IPP2, k):

    pp = 2

    rcc = 2.0 * IPP1['sigma_c']
    rcp = 2.0 * IPP1['sigma_c'] - IPP1['a']
    rpp = 2.0 * (IPP1['sigma_c'] - IPP1['a'])

    # EE
    O_EE = single_expo_weight(IPP1['r'], IPP2['r'], rcc, k)

    # EP
    O_EP = 0
    for i in range(1, pp+1):
        key = 'p%d' %i
        O_EP += single_expo_weight(IPP1['r'], IPP2[key]['r'], rcp, k)
        O_EP += single_expo_weight(IPP1[key]['r'], IPP2['r'], rcp, k)

    # EP
    O_PP = 0
    for i in range(1, pp+1):
        key1 = 'p%d' %i
        for j in range(1, pp+1):
            key2 = 'p%d' %j
            O_PP += single_expo_weight(IPP1[key1]['r'], IPP2[key2]['r'], rpp, k)

    return np.array([O_EE, O_EP, O_PP])


def single_expo_weight(ra, rb, rc, k):

    rab = np.linalg.norm(ra-rb)

    return np.exp(-k * (rab - rc))


# These functions switch representation between energy vectors: invert the matrix W and compute epsilon

def e_from_u(W, u):
    WI = np.linalg.inv(W)
    e = WI.dot(u)
    return e

def exp_cutoff(IPP1, u, k):

    rcc = 2.0 * IPP1['sigma_c']
    rcp = 2.0 * IPP1['sigma_c'] - IPP1['a']
    rpp = 2.0 * (IPP1['sigma_c'] - IPP1['a'])


    rc_ee = rcc - np.log(1e-3 / u[0]) / k
    rc_ep = rcp - np.log(1e-3 / np.fabs(u[1])) / k
    rc_pp = rpp - np.log(1e-3 / u[2]) / k

    rc = max([rc_ee, rc_ep, rc_pp])

    return rc

##########################################################################################
##########################################################################################
# MAIN

interaction = sys.argv[1]
if interaction == '0':
    weight = 'os'
elif interaction == '1':
    weight = 'exp'
else:
    print('UNKNOWN INTERACTION')
    raise SystemExit

## Basic paramters
sigma_c = 0.5                                    # Radius

'''
OVERLAP OF SPHERES
'''

VALUE_OF_DELTA = 0.2                             # This is delta MEASURED IN THE UNITS OF THE PARTICLE DIAMETER (NOT IN UNITS OF SIGMA_C!)
delta = VALUE_OF_DELTA / 2

## Option 1: a is set and sigma_p is computed -> IPC constraint is on

a = 0.22
sigma_p = sigma_c + delta - a
delta_c = delta
delta_p = delta

## Option 2: gamma is set and a, sigma_p are computed -> IPC constraint is on

# VALUE_OF_GAMMA = 40
# gamma = VALUE_OF_GAMMA * np.pi / 180
#
# cg = np.cos(gamma)
# a = delta*(2*sigma_c + delta) / (2*(sigma_c+delta-sigma_c*cg))
# sigma_p = sigma_c + delta - delta*(2*sigma_c + delta) / (2*(sigma_c+delta-sigma_c*cg))
# delta_c = delta
# delta_p = delta

## Option 3: IPC constrain is off, all parameters are independent

# a = 0.22
# sigma_p = 0.4
# delta_c = delta
# delta_p = sigma_p + a - sigma_c


'''
EXPONENTIAL MODEL
'''

a = 0.22
k = 13

'''
COMPUTING EPSILON VECTOR
'''

## All the input parameters
params = (a, sigma_c, sigma_p, delta_c, delta_p, k)

## Particle setting
IPP1 = set_IPP(params)
IPP2 = set_IPP(params)


# u vector in input, epsilon vector in output
u = np.array([0.1, -1, 4.0])
e = compute_e_vec(IPP1, IPP2, u, IPP1['a'], IPP1['sigma_c'], k, weight)

# Cutoff of the exponential model
r_cut = exp_cutoff(IPP1, u, k)

# Output file

if weight == 'os':
    model_name = 'EE%.1f_PP%.1f_D%.1f_A%.2f' %(u[0], u[2], VALUE_OF_DELTA, a)
elif weight == 'exp':
    model_name = 'EE%.1f_PP%.1f_D%.1f_K%d' %(u[0], u[2], VALUE_OF_DELTA, k)

F = open('./systems/%s.txt' %model_name, 'w')
F.write('a = %.12f\n' %a)
F.write('sigma_c = %f\n' %sigma_c)
F.write('sigma_p = %.12f\n' %sigma_p)
F.write('delta_c = %.12f\n' %delta_c)
F.write('delta_p = %.12f\n' %delta_p)
F.write('e_EE = %.12f\n' %e[0])
F.write('e_EP = %.12f\n' %e[1])
F.write('e_PP = %.12f\n' %e[2])
F.write('k_EE = %.12f\n' %k)
F.write('k_EP = %.12f\n' %k)
F.write('k_PP = %.12f\n' %k)
F.write('r_cut = %.12f\n' %r_cut)
F.close()
