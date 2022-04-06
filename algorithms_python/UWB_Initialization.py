# -*- coding: utf-8 -*-
"""
Algorithms for the UWB landmark initialization

@author: Sid Mahapatra
"""

import numpy as np

def MDS_loc_engine_3D(anchorCoordinates,Ranges):
    n,p = anchorCoordinates.shape
    k,m = Ranges.shape
    assert n==k, f"The number of anchors does not match the number of ranges"
    D_squared = calculate_distance_matrix(anchorCoordinates)
    D_full = np.block([[D_squared, Ranges**2],[Ranges.T**2,0]])
    C = np.eye(n+m)-np.ones((n+m,n+m))/(n+m)
    B = -0.5*C@D_full@C
    U,S,Vt = np.linalg.svd(B,full_matrices=True)
    P = U[:, :p]@np.sqrt(np.diag(S[:p]))
    
    R,t = Kabsch(P[:-1,:], anchorCoordinates)
    out = P[-1,:]@R+t
    
    return out
    
def Bancroft_loc_engine_3D(anchorCoordinates,Ranges):
    B = np.block([anchorCoordinates,-Ranges])
    B_ = np.linalg.pinv(B)
    a = 0.5*dotLorentz(B,B)
    c2 = dotLorentz(np.sum(B_,axis=1),np.sum(B_,axis=1))
    c1 = 2*(dotLorentz(B_@a,np.sum(B_,axis=1))-1)
    c0 = dotLorentz(B_@a,B_@a)
    if (c1*c1-4*c0*c2)>=0:
        root1 = (-c1+np.sqrt(c1*c1-4*c0*c2))/(2*c2)
        root2 = (-c1-np.sqrt(c1*c1-4*c0*c2))/(2*c2)
    else:
        root1 =-c1/(2*c2)
        root2 = root1  
        
    u1 = B_@(a+root1)
    u2 = B_@(a+root2)
    U =np.array([u1,u2])
    U = U[np.abs(U[:, 3]).argsort()]
    return U[0,:3]

def calculate_distance_matrix(A,B=None,squared=True):
    if B is None:
        B=A
    if B.ndim<2:
        B = B.reshape((1,-1)) 
    A=A
    B=B
    m = A.shape[0]
    n = B.shape[0]
    assert A.shape[1] == B.shape[1], f"The number of components for vectors in A \
        {A.shape[1]} does not match that of B {B.shape[1]}!"

    Adot = np.sum(A*A,axis=1).reshape((m,1))*np.ones(shape=(1,n))
    Bdot = np.sum(B*B,axis=1)*np.ones(shape=(m,1))
    ABdot = np.dot(A,B.T)
    D_squared =  Adot + Bdot -2*ABdot
    
    if squared:
        return D_squared
    else:
        D = np.sqrt(D_squared)
        return D
    
def Kabsch(P,Q):
    assert P.shape == Q.shape, f"shape of P and Q do not match"
    p0 = np.average(P,axis=0)
    q0 = np.average(Q,axis=0)
    P_centered = P-p0
    Q_centered = Q-q0
    C = P_centered.T@Q_centered
    U,S,Vt = np.linalg.svd(C)
    R = U@Vt
    t = q0-p0@R
    return R,t

def dotLorentz(u,v,dim=1):
    if u.ndim>1:
        if dim ==0:
            out = np.inner(u.T[:,:-1],v.T[:,:-1]).diagonal()-u[-1,:]*v[-1,:]
        elif dim ==1:
            out = np.inner(u[:,:-1],v[:,:-1]).diagonal()-u[:,-1]*v[:,-1]
        else:
            out = None
    else:
        out = np.inner(u[:-1],v[:-1])-u[-1]*v[-1]
    return out
    
def checkInitialization(anchorCoordinates,Ranges,poseEstimate,PdopThs=20,sigmaThs=0.1,planeRatioThs=1):
    n,p = anchorCoordinates.shape
    PdopThs = PdopThs/np.sqrt(n)
    
    A = (anchorCoordinates-poseEstimate)/Ranges
    Q = np.linalg.inv(A.T@A)
    PDOP = np.sqrt(np.trace(Q))
    
    RangesEstimate = np.linalg.norm(anchorCoordinates-poseEstimate,axis=1).reshape(-1,1)
    sigma = np.linalg.norm(RangesEstimate-Ranges)/np.sqrt(n-1)
    
    meanAnchorCoordinates = np.average(anchorCoordinates,axis=0)
    U,S,Vt = np.linalg.svd(anchorCoordinates-meanAnchorCoordinates)
    V = Vt.T
    q = poseEstimate-meanAnchorCoordinates
    angle = np.arcsin(np.abs(np.dot(V[:,-1],q))/np.linalg.norm(q))
    b = -np.linalg.norm(q)*np.sin(angle)+np.sqrt(np.sin(angle)**2*np.linalg.norm(q)**2+sigma*(2*np.linalg.norm(q)+sigma))
    planeRatio = S[-1]/np.sqrt(n)/b
    out = (planeRatio>planeRatioThs and PDOP<PdopThs and sigma<sigmaThs)
    return out
    
    
    
    
    
    
