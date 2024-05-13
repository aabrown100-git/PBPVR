"""
Physics-based pressure-volume relationship (PBPVR) model for heart research.
Date: 2023-02-22
Author: Yunxiao Zhang, Moritz Kalhöfer-Köchling, Eberhard Bodenschatz, and Yong Wang
Email: yunxiao.zhang@ds.mpg.de; yong.wang@ds.mpg.de
Cite it: 
"""

from pbpvr_help_functions import *

# Basis
# Normalized R
def Rn(d): 
    '''
    R_hat (Zhang_2023 Eq. 13). From Eq. 13, can show R_hat = 1+delta
    Input:
        d: normalized thickness variable; delta = (R - R_endo)/R_endo in Eq. 13
    '''
    return 1+d

# Normalized r
def rn(d,Vn): 
    '''
    r_hat (Zhang_2023 Eq. 14)
    Input:
        d: normalized thickness variable; delta = (R - R_endo)/R_endo in Eq. 13
        Vn: normalized volume variable; V = V/V0 in Eq. 13
    '''
    return (Rn(d)**3+Vn-1)**(1/3)

# radial strain lambda_rho
def lamd(d,Vn): 
    '''
    lambda_rho (Zhang_2023 Eq. 8)
    Input:
        d: normalized thickness variable; delta = (R - R_endo)/R_endo in Eq. 13
        Vn: normalized volume variable; V = V/V0 in Eq. 13
    '''
    return Rn(d)**2/rn(d,Vn)**2

# the first invariant of the right Cauchy-Green deformation tensor
def I1(lamd): 
    '''
    I_1 (Zhang_2023 Eq. 11)
    Input:
        lamd: lambda_rho in Eq. 8
    '''
    return (lamd**2+2/lamd)

# Passive Part
# P_ED part to integrate
def dW(d,Vn,a,b): 
    '''
    Integrand of passive part of PBPVR model (Zhang_2023 Eq. 19)
    Input:
        d: normalized thickness variable; delta = (R - R_endo)/R_endo in Eq. 13
        Vn: normalized volume variable; V = V/V0 in Eq. 13
        a: scaling material parameter a in Eq. 18
        b: exponential material parameter b in Eq. 18
    '''
    return  2*a * (1-lamd(d,Vn)**3)/rn(d,Vn) * exp(b*(I1(lamd(d,Vn))-3))

# P_ED part for single input
def moritz_curve(Vn,D,a,b):
    '''
    Performs integral over shell thickness in Eq. 19
    Input:
        Vn: normalized volume variable; V = V/V0 in Eq. 13
        D: normalized thickness; Delta = (R_epi - R_endo)/R_endo in Eq. 13
        a: scaling material parameter a in Eq. 18
        b: exponential material parameter b in Eq. 18
    '''
    # Vn = np.array(Vn)
    return quad(dW,0,D,args=(Vn,a,b))[0]

# P_ED part for array input
def vmoritz_curve(Vn_list:Vn2,D,a,b):
    '''
    Performs integral over shell thickness in Eq. 19 for array input
    Input:
        Vn_list: normalized volume variable; V = V/V0 in Eq. 13
        D: normalized thickness; Delta = (R_epi - R_endo)/R_endo in Eq. 13
        a: scaling material parameter a in Eq. 18
        b: exponential material parameter b in Eq. 18
    '''
    # Vn = np.array(Vn)
    return np.array([quad(dW,0,D,args=(Vn,a,b))[0] for Vn in Vn_list])

# P_ED part for array input
def moritz_curve_for_fitting(Vn1_list:Vn1,D,a:mmHg,b):
    Vn2_list=Vn1_to_Vn2_2(Vn1_list,V30_V0=V30_V0(D,a,b))
    return np.array([quad(dW,0,D,args=(Vn2,a,b))[0] for Vn2 in Vn2_list]) #kPa

# P_ED part for array: convert p to v
def vmoritz_curve_ptov(p:mmHg,D,a:mmHg,b)->Vn2:
    f = lambda Vn2: (vmoritz_curve(Vn2,D,a,b)-p)**2
    return minimize(f,1.6).x[0]

# Find Normalized V30 of 2nd method
def V30_V0(D,a:mmHg,b)->Vn2:
    f = lambda Vn2: (vmoritz_curve(Vn2,D,a,b)-30)**2
    return minimize(f,1.6).x[0]


# Active Part
# P_a part to integrate
def dW_a(d,Ta,Vn,lamd_0=1.58/1.85):
    '''
    Integrand of active part of PBPVR model (Zhang_2023 Eq. 25)
    Input:
        d: normalized thickness variable; delta = (R - R_endo)/R_endo in Eq. 13
        Ta: maximum active stress Ta in Eq. 20
        Vn: normalized volume variable; V = V/V0 in Eq. 13
        lamd_0: sarcomere length ratio in Eq. 20
    '''
    return 2*Ta * (1-lamd_0*sqrt(lamd(d,Vn)/2))/rn(d,Vn)

# P_a part for single input
def yunxiao_active(Vn,D,Ta,lamd_0):
    '''
    Performs integral over shell thickness in Eq. 25
    Input:
        Vn: normalized volume variable; V = V/V0 in Eq. 13
        D: normalized thickness; Delta = (R_epi - R_endo)/R_endo in Eq. 13
        Ta: maximum active stress Ta in Eq. 20
        lamd_0: sarcomere length ratio in Eq. 20
    '''
    return quad(dW_a,0,D,args=(Ta,Vn,lamd_0))[0]

# Full PBPVR model with scaling factor, for for single input
def PBPVR(Vn,D,a,b,Ta,p_passive,p_active,lamd_0=1.58/1.85):
    '''
    Physical-based pressure-volume relationship model for heart research.
    Given normalized volume (and material parameters), calculate pressure.
    Input:
        Vn: normalized volume variable; V = V/V0 in Eq. 13
        D: normalized thickness; Delta = (R_epi - R_endo)/R_endo in Eq. 13
        a: scaling material parameter a in Eq. 18
        b: exponential material parameter b in Eq. 18
        Ta: maximum active stress Ta in Eq. 20
        p_passive: passive pressure coefficient [0,1]
        p_active: active pressure coefficient [0,1]
        lamd_0: sarcomere length ratio in Eq. 20
    '''
    p = p_passive * moritz_curve(Vn,D,a,b) + p_active * yunxiao_active(Vn,D,Ta,lamd_0)
    return p

# Test
if __name__ == "__main__":
    x = [1,2,3]
    y = [2,4,7]
    a = to_equalx_data([1,2,3,4,5,4,3,32,2],x,y)
    print(a)
