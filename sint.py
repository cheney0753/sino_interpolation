# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 11:12:17 2017

@author: zhong

This is the module for SINT algorithm
"""
import numpy as np

class sint(S_known, thetas, dif_th):
    """ 
    Main algorithm of SINT
    """
    global L, N
    
    # N is the number of the pixels in one projection image
    N = S.shape[0]
    
    # L is the number of warps
    L = 0
    S_est=np.zeros(thetas, N)
    
    # Begin the main algorithm
    ka = np.ceil( diff_th / ( 2 * np.arcsin((1/(N-1))) ) )
    
    for h in range (1, N-2):
        for z in range(1 , k-1):
            
            th_h = th_h_m1 + dif_th
            
            [WM, L] = warps_generator(S, th_h)
            
            [A, IM] = relation_generator(S, WM, L)
            
            G = eignenvector_generator(A, WM)
            
            R = rule_generator(WM, IM)
            
            S_est[h, z] = sinocol_generator(A, G, R, S_k, eps)

    return S_est 


"""
-------------------------------------------------------------------------------
Algorithm 1: generation of the warp matrix
"""

def warps_generator(S, thetas, h):
    """
    Algorithm 1
    Generation of the warp matrix \PHI
    
    VARIABLES:
    S:          the known sinogram
    h:          column index of the missing column in the estimated sinogram
    thetas:     angles of the estimated sinogram
    """   
    
    # Define an empty warp matrix
    W_list = []    
    
    k = 0
        
    th_1 = thetas(h-1 )  # theta_(h-1), the angle at h-1 column
    th_3 = thetas(h+1)  # theta_(h+1)
    theta = thetas(h)   # angle of the missing column h
    
    for i_1 in range(0, N-1):
        
        for i_3 in range(0, N-1):
            
            # compute a vector of gray values (s_wp) corresponding to a warp
            
            s_wp = np.zeros([1, H])
            for j in range(0, H-1):
                theta_j = thetas(j)
                id_wp = np.round(warp_calculator(theta_j, i_1, i_3, th_1, th_3))
                s_wp = S[id_wp, j]
                
            s_bool = s_wp > 0   
            
            # Calculate the k_th warp using EQ 1.9
            w_k = warp_calculator( theta, i_1, i_3, th_1, th_3)
            if s_bool.all():
                i_2 = np.ceil(w_k)
                W_list.append([i_1, i_2, i_3])
                k += 1
                    
    WM = np.array(W_list)
    L = WM.shape[0] # L is the number of the warps
    
    return [WM , L]
    
def warp_calculator(theta, thetas, i_1, i_3, th_1, th_3):
    """
    Inplementation of EQ 1.9 for calculating a warp, which determine a warp 
    based on i_1 and i_3 and the corresponding angles th_1, th_3
    
    VARIABLES:
    i_1, i_3 : row index in the sinogram at -1 and +1 columns

    """
    q_1 = (N+1)/2 - i_1
    q_2 = (N+1)/2 - i_3
        
    # note that arccot is implemented as arctan here
    phi = np.arctan( ( q_2 * np.sin(th_1) - q_1*np.sin(th_3) ) 
                        / ( q_2 * np.cos(th_1) - q_1 * np.cos(th_3) ) )
                                                
    amp = q_2 / np.sin( th_3 - phi )
    
    w_k = warp_equation(theta, amp, phi)
    
    return w_k
    
def warp_equation(theta, amp, phi ):
    """
    Inplementation of EQ 1.9 for calculating a warp
    """
    w_k = amp * np.sin( theta - phi) + ( N + 1 ) / 2
    
    return w_k


"""
-------------------------------------------------------------------------------
Algorithm 2: generation of the relation matrix A and \PSi
"""    
def relation_generator(S_nb, WM, L):
    """
    Generation of the matrices A and \PSI (IM)
    VARIABLES:
    S_nb:   a N-by-2 matrix containing two columns S_(h-1) and S_(h+1)
    WM:     warp matrix
    L:      number of warps
    """
    c = 0
          
    A = np.zeros([N, L])
    
    IM = np.zeros([N, L], dtype = int)
    
    for j in range(1, N):
        
        for k in range(1, L):
            
            if WM[k,1] == j:
                IM[j, c:c+1] = [WM[k, 0], WM[k, 2]]
                A[j, c:c+1] = [S_nb[WM[k,0], 0], S_nb[WM[k,2], 1]]
                c += 2

    return [A, IM]


"""
-------------------------------------------------------------------------------
Algorithm 3: generation of the eignvector matrix G
"""   
def eignenvector_generator(S_nb, A, WM, L):
    """
    Generation of the eignvector matrix G
    VARIABLES:
    S_nb:   a N-by-2 matrix containing two columns S_(h-1) and S_(h+1)
    A:      relation matrix
    WM:     warp matrix
    L:      number of warps
    """    
    AA = np.zeros(2*L, 2*L)
    
    for i in range(0, L-1):
        d_1 = WM[i, 0]
        d_2 = WM[i, 2]
        AA[d_1, d_2] = S_nb[d_1,0] / S_nb[d_2, 1]
        AA[d_2, d_1] = S_nb[d_2,0] / S_nb[d_1, 1]
        
    # Calculate all the positive eigenvectors v of the matrix AA
    # G = [v(1), v(2) ... v(L)]  
    [w, EG] = np.linalg.eig(AA)
    
    G = np.array(2*L, L) #eigen vector matrix 
    k = 0 # Number of nonnegative eigenvectors for AA
    for i in range(0, np.size(w) -1 ):
        if all( EG[:, i] > 0):
            G[:, k] = EG[:, i] / w[i]
            k += 1
    
    if k != L:
        raise Exception('The number of positive eigen vectors is not L!')
                   
    return G
    
"""
-------------------------------------------------------------------------------
Algorithm 4: generation of the rule matrix R
"""      
def rule_generator(WM, IM):
    """
    Generation of the rule matrix R
    """
    k = 0
    R = 
    N_1 = size( S[h-1, :] != 0) 
    N_2 = size( S[h+1, :] != 0)
    
    R = np.zeros(N_1*N_2, L)
    
    for h in [1 ,3]:
        
        for d in range (0, N-1):
            
            for i in range(0, N-1):
                
                for j in range(0, 2L-1):
                    
                    if IM[i, j] == WM[d, h]:
                        R[k, j] = 1
                        
        k += 1

    return R

def sinocol_generator(A, G, R, S_k, eps):
    """
    Genration of the sinogram column S_h. 
    """
    beta = 0
    r = 0
    r_bar = eps
    
    while r < r_bar:

        beta = beta + eps
        RG = R*G
        temp = RG.transpose() * RG + beta* np.identity( RG.shape[0])
        
        S_hb =  A* np.linalg.inv(temp)* RG
        
        r_bar = r
        
        r = abs(S_hb.sum() - S_k.sum())
        
        if r < r_bar:
            S_h = S_hb
            
    return S_h
    