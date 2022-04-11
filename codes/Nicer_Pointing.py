from astropy.io import fits
from datetime import datetime
import geopack
import geopack.geopack as gp
import numpy as np
from dateutil import parser
import math
import pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

##### Import of data
# Data location is just a folder on my computer
RootPath = '/Users/emilatz/Dropbox/IDL_Lib/NICER_Cusp/NICER_obs_orbit/'
DataName = 'ni4202470103'
# fits.open is used to read the fits files
orb_fits = fits.open(RootPath+DataName+'.orb')
att_fits = fits.open(RootPath+DataName+'.att')
# the header contains relevant variabl information
orb_header = orb_fits[1].header
att_header = att_fits[1].header
# the data file contains the relevant data
orb_data = orb_fits[1].data
att_data = att_fits[1].data
# If you are searching for a variable, this is one way you can find it
# cols = att_data.columns
# print(cols)
# If you are searching for info on the header, this is how you print it
# print(repr(att_header))
# apple
# Getting the orbit position data
Times_Pos = orb_data['TIME']
X = orb_data['X']
Y = orb_data['Y']
Z = orb_data['Z']
Times_Q = att_data['TIME']
Quats = att_data['QPARAM']
# The TIME data is set as mission elapsed time since the start
# This needs to be calibrated to UTC time
NICERStart = datetime(year=2014,month=1,day=1,hour=0,minute=0,second=0)
NICERStart = datetime.timestamp(NICERStart)
LeapSec = 2
Times_Pos = Times_Pos + NICERStart - LeapSec
Times_Q = Times_Q + NICERStart - LeapSec

##### Data manipulation
# Unfortunatley the quaternion data is higher resolution than the position data
# This means that we need to interpolate one into the other to match data points
Q1_Func = interpolate.interp1d(Times_Q, Quats[:,0], kind='cubic')
Q2_Func = interpolate.interp1d(Times_Q, Quats[:,1], kind='cubic')
Q3_Func = interpolate.interp1d(Times_Q, Quats[:,2], kind='cubic')
Q4_Func = interpolate.interp1d(Times_Q, Quats[:,3], kind='cubic')

# Apply the times of the position data to the quaternion function to generate
# quaternions at every time of the position data
Q1_at_P = Q1_Func(Times_Pos)
Q2_at_P = Q2_Func(Times_Pos)
Q3_at_P = Q3_Func(Times_Pos)
Q4_at_P = Q4_Func(Times_Pos)

# Initialize some new variables for the calculations
RX = np.zeros((len(Times_Pos ),3))
RY = np.zeros((len(Times_Pos ),3))
RZ = np.zeros((len(Times_Pos ),3))
RA = np.zeros((len(Times_Pos ),1))
Dec = np.zeros((len(Times_Pos ),1))
Roll = np.zeros((len(Times_Pos ),1))
# Now convert the quaternions into angles
for i in range(0,len(Times_Pos)):
    q1 = Q1_at_P[i] #x
    q2 = Q2_at_P[i] #y
    q3 = Q3_at_P[i] #z
    q4 = Q4_at_P[i] #w This is the scalar

    # This is the calculation of the rotation matrix. This rotation matrix rotates
    # The coordinate system (like ECI) into the coordinates pointing of the subject.
    # Therefore, multiplying the rotation matrix by [[1],[0],[0]] gives the pointing
    # of the X axis, multiply by [[0],[1],[0]] gives pointing of Y axis.... etc
    # because each row of the rotation matrix is the unit vector by which the
    # the axis points
    RX[i,0] = +q1*q1 - q2*q2 - q3*q3 + q4*q4
    RX[i,1] = 2*(q1*q2 + q3*q4)
    RX[i,2] = 2*(q1*q3 - q2*q4)
    RY[i,0] = 2*(q1*q2 - q3*q4)
    RY[i,1] = -q1*q1 + q2*q2 - q3*q3 + q4*q4
    RY[i,2] = 2*(q2*q3 + q1*q4)
    RZ[i,0] = 2*(q1*q3 + q2*q4)
    RZ[i,1] = 2*(q2*q3 - q1*q4)
    RZ[i,2] = -q1*q1 - q2*q2 + q3*q3 + q4*q4
    # Calculation of the pointing angles from quaternions
    # This comes from code by Fred Eckert and Aspire - specifically Quaternion.cpp
    #####
    # Aspire - High-performance, lightweight IP-based messaging middleware.
    # *   Copyright (C) 2010-2011 Ken Center, Fred Eckert & Rod Green
    #####
    #THIS IS ROLL
    denom = q1*q3 - q2*q4 #Prep denominator
    psi = math.atan2(-(q2*q3 + q1*q4),denom) #Calculate arc tan with quadrant dependency
    psi = psi + 0
    # These if statements bring the roll into the positive rotation by wrapping 2pi
    if (psi < -np.pi):
        n = (np.ceil(-psi / (2*np.pi)))
        Roll[i] = psi+(n*2*np.pi)
    elif (psi >= np.pi):
        n = (psi / (2*np.pi));
        Roll[i] = psi-(n*2*np.pi)
    else:
        Roll[i] = psi
    #####
    # THIS IS RIGHT ASCENSION
    denom = q3*q1 + q2*q4 #Prep denominator
    phi = math.atan2(q3*q2 - q1*q4,denom) #Calculate arc tan with quadrant dependency
    # These if statements bring the RA into the positive rotation by wrapping 1pi
    if (phi < 0):
        n = (np.ceil(-phi / (2*np.pi)))
        RA[i] = phi+(n*2*np.pi)
    elif (phi >= 2*np.pi):
        n = (phi / (2*np.pi));
        RA[i] = phi-(n*2*np.pi)
    else:
        RA[i] = phi
    #####
    # THIS IS DECLINATION
    acos = q3*q3 + q4*q4 - q2*q2 - q1*q1 #This is actually an element of a directional cosine matrix
    if (acos > 1):
        acos = 1
    if (acos < -1):
        acos = -1
    theta = math.acos(acos) #get the angle
    # These if statements bring the Declination to be defined from the equator
    if (theta >=0):
        Dec[i] = np.pi/2 - theta
    else:
        Dec[i] = -np.pi/2 - theta

# Convert the RA and Dec to a point in ECI
x_eci = np.cos(RA)*np.cos(Dec)
y_eci = np.sin(RA)*np.cos(Dec)
z_eci = np.sin(Dec)
# Summarize it in a column
RaDec = np.c_[x_eci,y_eci,z_eci]
# Initialize some variable matricies
P_gsm = np.zeros((len(Times_Pos),3))
RaDec_gsm = np.zeros((len(RaDec[:,0] ),3))

###### Conversion to GSM
# Iterate through the matricies to convert to GSM with geopack.
for i in range(len(Times_Pos)):
    gp.recalc(Times_Pos[i]) # set the time in geopack to get proper location of ECI
    tmp1,tmp2,tmp3 = gp.geigeo(X[i],Y[i],Z[i], 1) # the 1 means go from GEI to GEO
    P_gsm[i,0],P_gsm[i,1],P_gsm[i,2] = gp.geogsm(tmp1,tmp2,tmp3, 1) # from GEO to GSM
    # P_gse[i,0],P_gse[i,1],P_gse[i,2] = gp.gsmgse(P_gsm[i,0],P_gsm[i,1],P_gsm[i,2], 1)

    tmp1,tmp2,tmp3 = gp.geigeo(RaDec[i,0],RaDec[i,1],RaDec[i,2], 1)
    RaDec_gsm[i,0],RaDec_gsm[i,1],RaDec_gsm[i,2] = gp.geogsm(tmp1,tmp2,tmp3, 1)
    # RaDec_gse[i,0],RaDec_gse[i,1],RaDec_gse[i,2] = gp.gsmgse(RaDec_gsm[i,0],RaDec_gsm[i,1],RaDec_gsm[i,2], 1)

P_gsm = P_gsm/1000 #change from meters to km

##### Plotting

fig = plt.figure()
ax =  plt.axes(projection = '3d')
ax.set_title('NICER Pointing GSM')
ax.set_xlabel('GSM X (km)')
ax.set_ylabel('GSM Y (km)')
ax.set_zlabel('GSM Z (km)')
ax.set_xlim(-7000.0, 7000.0)
ax.set_ylim(-7000.0, 7000.0)
ax.set_zlim(-7000.0, 7000.0)
ax.view_init(20, 45)
ax.plot([0,7000], [0,0], [0,0], color="y") # plot the Earth-sun lin
ax.plot3D(P_gsm[:,0], P_gsm[:,1], P_gsm[:,2], color = 'r',label='ISS orbit') #Scatter the points of the sc position

#Plot an Earth sphere (thank you StackOverflow)
# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# Sx = np.cos(u)*np.sin(v)*6371
# Sy = np.sin(u)*np.sin(v)*6371
# Sz = np.cos(v)*6371
# ax.plot_surface(Sx, Sy, Sz, color="b",linewidth=0.0)
# ax.plot_surface(x, y, z, linewidth=0.0)

# Plot the first vector so you get the label
ax.quiver(P_gsm[0,0],P_gsm[0,1],P_gsm[0,2],RaDec_gsm[0,0],RaDec_gsm[0,1],RaDec_gsm[0,2],length=1600,color="k",label='NICER view')

# For loop to go over the pointing vectors
for i in range(0,len(P_gsm[:,0])):
    if i % 20 == 0: # use only every n'th value. It gets confusing with too many
        ax.quiver(P_gsm[i,0],P_gsm[i,1],P_gsm[i,2],RaDec_gsm[i,0],RaDec_gsm[i,1],RaDec_gsm[i,2],length=1600,color="k")

ax.legend()
plt.show()
fig.savefig(RootPath+'OrbitPlot_GSM.png')
plt.close()

print("Doin Science Yo!")
