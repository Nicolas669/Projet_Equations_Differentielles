_author_ = "Nicolas Coruzzi"
_filename_ = "Exercice_3"
_creationdate_ = "09/05/20"

import matplotlib.pyplot as plt
import numpy as np

mu=10
Tmax=4
N=4000
temps=[i*((Tmax-0)/(N)) for i in range(N+1)]
q0=[0,3]
p0=[1,0]

def dq_H(q,p):
    q1=q[0]
    q2=q[1]
    return np.array([(mu*q1)/((q1*q1+q2*q2)**(3/2)),(mu*q2)/((q1*q1+q2*q2)**(3/2))])

def dp_H(q,p):
    p1=p[0]
    p2=p[1]
    return np.array([p1,p2])

def dq_H_sep(q):
    #a utiliser sur qn+1
    q1=q[0]
    q2=q[1]
    return np.array([(mu*q1)/((q1*q1+q2*q2)**(3/2)),(mu*q2)/((q1*q1+q2*q2)**(3/2))])

def dp_H_sep(p):
    #a utiliser sur pn
    p1=p[0]
    p2=p[1]
    return np.array([p1,p2])

def euler_hamiltonian_solve(dq_H,dp_H,p0,q0,temps):
    q=np.zeros((2,N+1))
    q[:,0] =q0
    p=np.zeros((2,N+1))
    p[:,0]=p0
    for n in range(N):
        tn=temps[n]
        tn1=temps[n+1]
        dtn=tn1-tn
        q[:,n+1]=q[:,n]+dtn*dp_H(q[:,n],p[:,n])
        p[:,n+1]=p[:,n]-dtn*dq_H(q[:,n],p[:,n])
    return [q,p]

def rk4_hamiltonian_solve(dq_H,dp_H,p0,q0,temps):
    q=np.zeros((2,N+1))
    q[:,0] =q0
    p=np.zeros((2,N+1))
    p[:,0]=p0
    for n in range(N):
        tn=temps[n]
        tn1=temps[n+1]
        dtn=tn1-tn
        [kq0,kp0]=[dp_H(q[:,n],p[:,n]),-dq_H(q[:,n],p[:,n])]
        [kq1,kp1]=[dp_H(q[:, n]+(dtn/2)*kq0,p[:,n]+(dtn/2)*kp0),-dq_H(q[:,n]+(dtn/2)*kq0,p[:,n]+(dtn/2)*kp0)]
        [kq2,kp2]=[dp_H(q[:,n]+(dtn/2)*kq1,p[:,n]+(dtn/2)*kp1),-dq_H(q[:,n]+(dtn/2)*kq1,p[:,n]+(dtn/2)*kp1)]
        [kq3,kp3]=[dp_H(q[:,n]+dtn*kq2,p[:,n]+dtn*kp2),-dq_H(q[:,n]+dtn*kq2,p[:,n]+dtn*kp2)]
        q[:,n+1]=q[:,n]+(dtn/6)*(kq0+2*kq1+2*kq2+kq3)
        p[:,n+1]=p[:,n]+(dtn/6)*(kp0+2*kp1+2*kp2+kp3)
    return [q,p]

def euler_symplect_hamiltonian_separable_solve(dq_H_sep,dp_H_sep,p0,q0,temps):
    q=np.zeros((2,N+1))
    q[:,0]=q0
    p=np.zeros((2,N+1))
    p[:,0]=p0
    for n in range(N):
        tn=temps[n]
        tn1=temps[n+1]
        dtn=tn1-tn
        q[:,n+1]=q[:,n]+dtn*dp_H_sep(p[:, n])
        p[:,n+1]=p[:,n]-dtn*dq_H_sep(q[:,n+1])
    return [q, p]

def verlet_hamiltonian_separable_solve(dq_H_sep,dp_H_sep,p0,q0,temps):
    q=np.zeros((2,N+1))
    q[:,0]=q0
    p=np.zeros((2,N+1))
    p[:,0]=p0
    for n in range(N):
        tn=temps[n]
        tn1=temps[n+1]
        dtn=tn1-tn
        q_n_1_2=q[:,n]+(dtn/2)*dp_H_sep(p[:,n])
        p[:,n+1]=p[:,n]-dtn*dq_H_sep(q_n_1_2)
        q[:,n+1]=q_n_1_2+(dtn/2)*dp_H_sep(p[:,n+1])
    return [q, p]

def energy(q,p):
    return (((1/2)*(p[0]*p[0]+p[1]*p[1]))-(mu/((q[0]*q[0]+q[1]*q[1])**(1/2))))

(Q_euler,P_euler)=euler_hamiltonian_solve(dq_H,dp_H,p0,q0,temps)
q1_euler=Q_euler[0:1,:][0]
q2_euler=Q_euler[1:2,:][0]

(Q_rk4,P_rk4)=rk4_hamiltonian_solve(dq_H,dp_H,p0,q0,temps)
q1_rk4=Q_rk4[0:1,:][0]
q2_rk4=Q_rk4[1:2,:][0]

(Q_euler_symplect,P_euler_symplect)=euler_symplect_hamiltonian_separable_solve(dq_H_sep,dp_H_sep,p0,q0,temps)
q1_euler_symplect=Q_euler_symplect[0:1,:][0]
q2_euler_symplect=Q_euler_symplect[1:2,:][0]

(Q_verlet,P_verlet)=verlet_hamiltonian_separable_solve(dq_H_sep,dp_H_sep,p0,q0,temps)
q1_verlet=Q_verlet[0:1,:][0]
q2_verlet=Q_verlet[1:2,:][0]


plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Courbe parametrique N4000 T4',fontsize=18)
plt.xlabel('temps',fontsize=14)
plt.ylabel('q',fontsize=14)
plt.plot(temps,q1_euler,label='Euler q1(t)')
plt.plot(temps,q2_euler,label='Euler q2(t)')
plt.plot(temps,q1_rk4,label='Rk4 q1(t)')
plt.plot(temps,q2_rk4,label='Rk4 q2(t)')
plt.plot(temps,q1_euler_symplect,label='Euler symp q1(t)')
plt.plot(temps,q2_euler_symplect,label='Euler symp q2(t)')
plt.plot(temps,q1_verlet,label='Verlet q1(t)')
plt.plot(temps,q2_verlet,label='Verlet q2(t)')
plt.legend()
plt.savefig('Courbe_parametrique_N4000_T4')


plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Comparaison Hamiltonien N4000 T4',fontsize=18)
plt.xlabel('q1',fontsize=14)
plt.ylabel('q2',fontsize=14)
plt.plot(q1_euler,q2_euler,label='euler')
plt.plot(q1_rk4,q2_rk4,label='rk4')
plt.plot(q1_euler_symplect,q2_euler_symplect,label='euler symp')
plt.plot(q1_verlet,q2_verlet,label='verlet')
plt.legend()
plt.savefig('Comparaison_Hamiltonien_N4000_T4')
#plt.show()

H_euler=[]
H_rk4=[]
H_euler_symplect=[]
H_verlet=[]
for i in range(len(temps)):
    H_euler+=[energy(Q_euler[:,i:i+1],P_euler[:,i:i+1])]
    H_rk4+=[energy(Q_rk4[:,i:i+1],P_rk4[:,i:i+1])]
    H_euler_symplect+=[energy(Q_euler_symplect[:,i:i+1],P_euler_symplect[:,i:i+1])]
    H_verlet+=[energy(Q_verlet[:,i:i+1],P_verlet[:,i:i+1])]

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Conservation energie Hamiltonien N4000 T4',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Energie',fontsize=14)
plt.plot(temps,H_euler,label='euler')
plt.plot(temps,H_rk4,label='rk4')
plt.plot(temps,H_euler_symplect,label='euler symp')
plt.plot(temps,H_verlet,label='verlet')
plt.legend()
plt.savefig('Conservation_energie_Hamiltonien_N4000_T4')
plt.show()