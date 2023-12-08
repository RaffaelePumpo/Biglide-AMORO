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
def compute_A_B(q11, q12, q21, q22):
	u1 = np.array([math.cos(q12),math.sin(q12)])
	u2 = np.array([math.cos(q22),math.sin(q22)])
	y0 = np.array([0.0, 1.0])
	A = np.array([u1,u2])
	B = np.array([[u1.dot(y0), 0.0], [0.0, u2.dot(y0)]])
	return A, B

def dkm(q11, q12, q21, q22, q11D, q21D):
	qd = np.array([q11D, q21D])
	[A, B] = compute_A_B(q11, q12, q21, q22)
	Ainv = np.linalg.inv(A)
	J = np.matmul(Ainv, B)
	posD = J.dot(qd)
	xD = posD[0]
	yD = posD[1]
	return xD, yD

def ikm(q11, q12, q21, q22, xD, yD):
	posD = np.array([xD, yD])
	[A, B] = compute_A_B(q11, q12, q21, q22)
	Binv = np.linalg.inv(B)
	qda = np.matmul(Binv, A)
	qd = np.matmul(qda, posD)
	q11D = qd[0]
	q21D = qd[1]
	return q11D, q21D

def dkm_passive(q11, q12, q21, q22, q11D, q21D, xD, yD):
	posD = np.array([xD, yD])
	v1 = np.array([-math.sin(q12),math.cos(q12)])
	v2 = np.array([-math.sin(q22),math.cos(q22)])
	y0 = np.array([0.0, 1.0])
	A = np.array([v1, v2])
	B = np.array([[v1.dot(y0), 0.0], [0.0, v2.dot(y0)]])
	A_posD = np.matmul(A, posD)
	qaD = np.array([q11D, q21D])
	B_qa = np.matmul(B,qaD)
	
	qdD = 1/l*(A_posD-B_qa)
	
	q12D = qdD[0]
	q22D = qdD[1]
	return q12D, q22D

def dkm2(q11, q12, q21, q22, q11D, q12D, q21D, q22D, q11DD, q21DD):
	qaDD = np.array([q11DD, q21DD])
	[A, B] = compute_A_B(q11, q12, q21, q22)
	Ainv = np.linalg.inv(A)
	u1 = np.array([math.cos(q12),math.sin(q12)])
	u2 = np.array([math.cos(q22),math.sin(q22)])

	d = np.array([-l*(q12D)**2, -l*(q22D)**2])
	
	B_qaDD = np.matmul(B,qaDD)
	F = B_qaDD + d
	posDD = np.matmul(Ainv,F)
	
	xDD = posDD[0]
	yDD = posDD[1]
	return xDD, yDD

def ikm2(q11, q12, q21, q22, q11D, q12D, q21D, q22D, xDD, yDD):
	posDD = np.array([xDD, yDD])
	[A, B] = compute_A_B(q11, q12, q21, q22)
	Binv = np.linalg.inv(B)
	u1 = np.array([math.cos(q12),math.sin(q12)])
	u2 = np.array([math.cos(q22),math.sin(q22)])

	d = np.array([-l*(q12D)**2, -l*(q22D)**2])
	
	A_posDD = np.matmul(A, posDD)
	F = A_posDD - d
	
	qaDD = np.matmul(Binv, F)
	q11DD = qaDD[0]
	q21DD = qaDD[1]
	return q11DD, q21DD

def dkm2_passive(q11, q12, q21, q22, q11D, q12D, q21D, q22D, q11DD, q21DD, xDD, yDD):
	posDD = np.array([xDD, yDD])
	qaDD = np.array([q11DD, q21DD])

	u1 = np.array([math.cos(q12),math.sin(q12)])
	u2 = np.array([math.cos(q22),math.sin(q22)])
	v1 = np.array([-math.sin(q12),math.cos(q12)])
	v2 = np.array([-math.sin(q22),math.cos(q22)])
	y0 = np.array([0.0, 1.0])

	A = np.array([v1, v2])
	B = np.array([[v1.dot(y0), 0.0], [0.0, v2.dot(y0)]])

	A_posDD = np.matmul(A, posDD)
	B_qaDD = np.matmul(B,qaDD)

	qdDD = 1/l*(A_posDD-B_qaDD)

	q12DD = qdDD[0]
	q22DD = qdDD[1]
	return q12DD, q22DD

def dynamic_model(q11, q12, q21, q22, q11D, q12D, q21D, q22D):
	
	M1 = np.array([[mf, 0.0], [0.0, mf]])
	
	A, B = compute_A_B(q11, q12, q21, q22)
	Ainv = np.linalg.inv(A)
	J = np.matmul(Ainv,B)
	
	[xD,yD] = dkm(q11, q12, q21, q22, q11D, q21D)
	posD = np.array([xD,yD])

	qaD = np.array([q11D,q21D])

	AD = np.array([[2*xD,2*(yD-q11D)],[2*xD,2*(yD-q21D)]])
	BD = np.array([[-2*(yD-q11D),0],[0,-2*(yD-q21D)]])	
	
	K = np.matmul(AD,posD) + np.matmul(BD,qaD)
	b = np.matmul(-Ainv,K)
	
	M = M1+mp*np.matmul(J.transpose(),J)
	
	c = mp*np.matmul(J.transpose(),b)
	
	return M, c