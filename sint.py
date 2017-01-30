# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 11:12:17 2017

@author: zhong

This is the module for SINT algorithm
"""
import numpy as np

def sint(S_known, thetas_known):
    """ 
    Main algorithm of SINT
    """   
    # N is the number of the pixels in one projection image
    N = S_known.shape[0]
    
    incr_th_ori = thetas_known[1] - thetas_known[0]
    
    # Begin the main algorithm
    ka = int( np.ceil( incr_th_ori / ( 2 * np.arcsin((1/(N-1)) ) ) ) )

    S_est = np.zeros([ka*(S_known.shape[1]-1), S_known.shape[0]])    
    eps =0.01    
    
    for h in range (1, N-1):
        
        for z in range(1 , ka-1):
            
            th_1 = thetas_known[h-1]
            th_3 = thetas_known[h]
            
            th_int = th_1 + ( z / ka ) * incr_th_ori
                        
            WM = warps_generator(S_known,thetas_known, th_1, th_3, th_int)
            
            [A, IM] = relation_generator(S_known[:, h-1:h+1], WM)
            
            G = eignenvector_generator(S_known[:, h-1:h+1], A, WM)
            
            R = rule_generator(S_known[:, h-1:h+1],WM, IM)
            
            S_est[h, z] = sinocol_generator(A, G, R, S_known, eps)

    return S_est 


"""
-------------------------------------------------------------------------------
Algorithm 1: generation of the warp matrix
"""

def warps_generator(S, thetas, th_1, th_3, th_2):
    """
    Algorithm 1
    Generation of the warp matrix \PHI
    
    Parameters:
    S:          the known sinogram
    h:          column index of the missing column in the estimated sinogram
    thetas:     angles of the estimated sinogram
    """   
    
    # Define an empty warp matrix
    W_list = []    
    
    k = 0
       
    N =  S.shape[0]     # number of pixels  
    
    H = S.shape[1]  # number of projection images
    
    for i_1 in range(0, N):
        
        for i_3 in range(0, N):
            
            # determine a warp sino function given i_1 and i_3
            [amp, phi] = warp_equation(i_1+1 , i_3+1 , th_1, th_3, N)
            
            # use the warp only if the warp does not exceed the sinogram
            if abs(amp) > (N-1)/2 :            
                continue
            
            # testify if all the pixels passed by the warp function are positive
            id_wp = []            
            s_wp = []
            for j in range(0, H):
                id_wp.append( int( np.round( warp_calculator(thetas[j], amp, phi, N) ))-1 )
                s_wp.append( S[id_wp[j], j])

            s_wp_arr = np.array(s_wp)                    
            # id_wp_arr = np.array(id_wp)
            

            # id_wp_bool = (id_wp_arr > 0) & (id_wp_arr < N-1)
                        
            # compute a vector of gray values (s_wp) corresponding to a warp
            s_bool = s_wp_arr > 0
 
            # Calculate the k_th warp using EQ 1.9
            w_k = warp_calculator( th_2, amp, phi, N)
            if s_bool.all():
                i_2 = np.round(w_k) - 1
                W_list.append([i_1, i_2, i_3])
                k += 1
                #DEBUG
                        
    WM = np.array(W_list, dtype = int)
    L = WM.shape[0] # L is the number of the warps
    
    print('{} warps found'.format(L))   

    return WM 
    
def warp_equation(i_1, i_3, th_1, th_3, N):
    """
    Inplementation of EQ 1.9 for calculating a warp, which determine a warp 
    based on i_1 and i_3 and the corresponding angles th_1, th_3
    
    VARIABLES:
    i_1, i_3 : row index in the sinogram at -1 and +1 columns

    """
    q_1 = (N+1)/2 - i_1
    q_2 = (N+1)/2 - i_3
        
    # note that arccot is implemented as arctan here
    phi = np.arctan( ( q_2 * np.sin(th_1) - q_1 * np.sin(th_3) ) 
                        / ( q_2 * np.cos(th_1) - q_1 * np.cos(th_3) ) )
                                                
    amp = q_2 / np.sin( th_3 - phi )
        
    return [amp, phi]
    
def warp_calculator(theta, amp, phi, N ):
    """
    Inplementation of EQ 1.9 for calculating a warp
    
    return: 
    pixel index of the warp
    """
    w_k = -amp * np.sin( theta - phi) + ( N + 1 ) / 2
    
    return w_k


"""
-------------------------------------------------------------------------------
Algorithm 2: generation of the relation matrix A and \PSi
"""    
def relation_generator(S_nb, WM):
    """
    Generation of the matrices A and \PSI (IM)
    VARIABLES:
    S_nb:   a N-by-2 matrix containing two columns S_(h-1) and S_(h+1)
    WM:     warp matrix
    L:      number of warps
    """
    L = WM.shape[0]    
    N = S_nb.shape[0]
    
    c = 0
          
    A = np.zeros([N, 2*L])
    
    IM = np.zeros([N, 2*L], dtype = int)
    
    for j in range(0, N):
        
        for k in range(0, L):
            
            if WM[k,1] == j:
                IM[j, c:c+2] = [WM[k, 0], WM[k, 2]]
                A[j, c:c+2] = [S_nb[WM[k,0], 0], S_nb[WM[k,2], 1]]
                c += 2

    return [A, IM]


"""
-------------------------------------------------------------------------------
Algorithm 3: generation of the eignvector matrix G
"""   
def eignenvector_generator(S_nb, A, WM):
    """
    Generation of the eignvector matrix G
    VARIABLES:
    S_nb:   a N-by-2 matrix containing two columns S_(h-1) and S_(h+1)
    A:      relation matrix
    WM:     warp matrix
    L:      number of warps
    """   
    WM = WM[0:5, :]    
    
    L = WM.shape[0]    
    
    AA = np.zeros([2*L, 2*L])
    
    for i in range(0, L):
        d_1 = WM[i, 0]
        d_2 = WM[i, 2]
        AA[d_1, d_2] = S_nb[d_1,0] / S_nb[d_2, 1]
        AA[d_2, d_1] = S_nb[d_2,1] / S_nb[d_1, 0]
        
    # Calculate all the positive eigenvectors v of the matrix AA
    # G = [v(1), v(2) ... v(L)]  
    [w, EG] = np.linalg.eig(AA)
    
    G = np.zeros([2*L, L]) #eigen vector matrix 
    k = 0 # Number of nonnegative eigenvectors for AA
    for i in range(0, np.size(w)):
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
def rule_generator(S_nb, WM, IM):
    """
    Generation of the rule matrix R
    """
    L = WM.shape[0]    
    N = S_nb.shape[0]

    k = 0
    N_1 = np.sum(S_nb[:, 0] != 0) 
    N_2 = np.sum(S_nb[:, 1] != 0)
    
    R = np.zeros(N_1*N_2, L)
    
    for h in [1 ,3]:
        
        for d in range (0, N):
            
            for i in range(0, N):
                
                for j in range(0, 2*L):
                    
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
        temp = np.matmul( RG.transpose(), RG) + beta* np.identity( RG.shape[0])
        
        S_hb =  np.matmul(A.matmul( np.linalg.inv(temp)), RG)
        
        r_bar = r
        
        r = abs(S_hb.sum() - S_k.sum())
        
        if r < r_bar:
            S_h = S_hb
            
    return S_h
    