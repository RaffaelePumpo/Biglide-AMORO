import math
import numpy as np

# Geometric parameters
l = 0.2828427
d = 0.4
# Dynamic parameters
mp = 3.0
mf = 1.0


def dgm(q11, q21, assembly_mode):
    a2a1 = np.array([[-d],
                     [q11-q21]])
    a2h = 1/2*a2a1
    a = np.linalg.norm(a2h)
    h = math.sqrt(l**2 - a**2)
    u  = np.array([[0, -1],
                     [1, 0]])
    hc = assembly_mode*(h/a) * (np.matmul(u, a2h))
    p = np.array([[d/2],
                  [q21]]) + a2h + hc                           
    x = p[0, 0]
    y = p[1, 0]
    return x,y


def igm(x, y, gamma1, gamma2):
    q11 = y + gamma1* math.sqrt(-(x+d/2)**2 + l**2)
    q21 = y + gamma2* math.sqrt(-(x-d/2)**2 + l**2)
    return q11, q21


def dgm_passive(q11, q21, assembly_mode):
    x,y = dgm(q11, q21, assembly_mode)
    q12 = math.atan2(y-q11,x+d/2)
    q22 = math.atan2(y-q21,x-d/2)
    return q12, q22


# You can create intermediate functions to avoid redundant code
def compute_A_B(q11,q21):
	x,y = dgm(q11, q21, -1)
	A = np.array([[2*(x+d/2), 2*(y-q11)],
				[2*(x-d/2), 2*(y-q21)]])
	B = np.array([[-2*(y-q11), 0],
				[0, -2*(y-q21)]])
	return A,B


def dkm(q11,q12,q21,q22,q11D, q21D):
	qD = np.array([q11D,q21D])
	A,B = compute_A_B(q11, q21)
	pD = (np.matmul(-np.linalg.inv(A),B))@qD
	xD = pD[0]
	yD = pD[1]
	return xD,yD


def ikm(q11, q12, q21, q22, xD, yD):
	pD = np.array([xD,yD])
	A,B = compute_A_B(q11, q21)
	qD = (np.matmul(-np.linalg.inv(B),A))@pD
	q11D = qD[0]
	q21D = qD[1]
	return q11D, q21D


def compute_Ap_Bp(q11,q12,q21,q22):
	Ap = np.array([[-l*math.sin(q12), l*math.sin(q22)],
				[l*math.cos(q12), -l*math.cos(q22)]])
	Bp = np.array([[0, 0],
				[1, -1]])
	return Ap,Bp

def dkm_passive(q11, q12, q21, q22, q11D, q21D, xD, yD):
	qD = np.array([q11D,q21D])
	Ap,Bp = compute_Ap_Bp(q11,q12,q21,q22)
	pD = (np.matmul(-np.linalg.inv(Ap),Bp))@qD
	q12D = pD[0]
	q22D = pD[1]
	return q12D, q22D


def compute_AD_BD(q11,q12,q21,q22,q11D, q21D):
	xD,yD = dkm(q11,q12,q21,q22,q11D, q21D)
	A = np.array([[2*xD, 2*(yD-q11D)],
				[2*xD, 2*(yD-q21D)]])
	B = np.array([[-2*(yD-q11D), 0],
				[0, -2*(yD-q21D)]])
	return A,B

def dkm2(q11, q12, q21, q22, q11D, q12D, q21D, q22D, q11DD, q21DD):
	A,B = compute_A_B(q11, q21)
	qD =np.array([q11D,q21D])
	qDD =np.array([q11DD,q21DD])
	AD,BD = compute_AD_BD(q11,q12,q21,q22,q11D, q21D)
	xD,yD = dkm(q11,q12,q21,q22,q11D, q21D)
	pD =np.array([xD,yD])
	BqDD = B@qDD
	ADpD = AD@pD
	BDqD = BD@qD
	PDD = (np.matmul(-np.linalg.inv(A),BqDD+ADpD+BDqD))
	xDD = PDD[0]
	yDD = PDD[1]
	return xDD, yDD


def ikm2(q11, q12, q21, q22, q11D, q12D, q21D, q22D, xDD, yDD):
	A,B = compute_A_B(q11, q21)
	qD =np.array([q11D,q21D])
	pDD =np.array([xDD,yDD])
	AD,BD = compute_AD_BD(q11,q12,q21,q22,q11D, q21D)
	xD,yD = dkm(q11,q12,q21,q22,q11D, q21D)
	pD =np.array([xD,yD])
	ApDD = A@pDD
	ADpD = AD@pD
	BDqD = BD@qD
	qDD = (np.matmul(-np.linalg.inv(B),ApDD+BDqD+ADpD))
	q11DD = qDD[0]
	q21DD = qDD[1]
	return q11DD, q21DD

def compute_ApD_BpD(q11,q12,q21,q22,q11D, q21D,xD,yD):				
	q12D,q22D = dkm_passive(q11, q12, q21, q22, q11D, q21D, xD, yD)
	ApD = np.array([[-l*math.cos(q12)*q12D,l*math.cos(q22)*q22D],
				[-l*math.sin(q12)*q12D,l*math.sin(q22)*q22D]])
	BpD = np.array([[0,0],
				[0,0]])
	return ApD,BpD	
	
def dkm2_passive(q11, q12, q21, q22, q11D, q12D, q21D, q22D, q11DD, q21DD, xDD, yDD):
    Ap,Bp = compute_Ap_Bp(q11,q12,q21,q22)
    qD =np.array([q11D,q21D])
    qDD =np.array([q11DD,q21DD])
    xD,yD = dkm(q11,q12,q21,q22,q11D, q21D)
    ApD,BpD = compute_ApD_BpD(q11,q12,q21,q22,q11D, q21D,xD,yD)
    qpD =np.array([q12D,q22D])
    ApDqpD = ApD@qpD
    BpDqD = BpD@qD
    BpqDD = Bp@qDD
    qpDD = (np.matmul(-np.linalg.inv(Ap),ApDqpD+BpDqD+BpqDD))
    q12DD = qpDD[0]
    q22DD = qpDD[1]
    return q12DD, q22DD


def dynamic_model(q11, q12, q21, q22, q11D, q12D, q21D, q22D):
	A,B = compute_A_B(q11,q21)
	J = np.matmul(-np.linalg.inv(A),B)
	Jt = J.transpose()
	xD,yD = dkm(q11,q12,q21,q22,q11D, q21D)
	pD =np.array([xD,yD])
	AD,BD = compute_AD_BD(q11,q12,q21,q22,q11D, q21D)
	qD =np.array([q11D,q21D])
	ADpD = AD@pD
	BDqD = BD@qD
	b = np.matmul(-np.linalg.inv(A),ADpD+BDqD)
	M = np.zeros((2, 2))
	M1 = mp*np.matmul(Jt,J)
	Mf = np.array([[mf,0],
				[0,mf]])
	M = M1+Mf
	c = mp*Jt@b
	return M, c
