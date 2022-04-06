#!/usr/bin/env python
# coding: utf-8

# # Test of UWB_initialization functions
# 
# Sid Mahapatra

from UWB_Initialization import *
import numpy as np

#%% Define a probem with 4 robot poses and one anchor
RobotPoses = np.array([[0,0,0],[2,0,0],[0,2,0],[1,1,2]])
AnchorPose = np.array([3,3,3])
Ranges = calculate_distance_matrix(RobotPoses,AnchorPose,squared=False)
print('Ranges')
print(Ranges)

#%% Solve the initialization problem using the MDS and Bancroft algorithm
solutionMDS  = MDS_loc_engine_3D(RobotPoses, Ranges)
solutionBanc = Bancroft_loc_engine_3D(RobotPoses, Ranges)
print('Solution of MDS')
print(solutionMDS)
print('Solution of Bancroft')
print(solutionBanc)

#%%Add disturbance to the ranges
sigma = 0.03
DisturbedRanges = Ranges + np.random.normal(0, sigma, Ranges.size).reshape((-1,1))
print('Disturbed Ranges')
print(DisturbedRanges)

#%%Solve the disturbed initialization problem using the MDS and Bancroft algorithm
solutionMDS_dist  = MDS_loc_engine_3D(RobotPoses, DisturbedRanges)
solutionBanc_dist = Bancroft_loc_engine_3D(RobotPoses, DisturbedRanges)
print('Solution of MDS')
print(solutionMDS_dist)
print('Solution of Bancroft')
print(solutionBanc_dist)

#%%Check the initialization function

tf = checkInitialization(RobotPoses,DisturbedRanges,solutionMDS_dist)
