_author_ = "Nicolas Coruzzi"
_filename_ = "Exercice_2_pendule"
_creationdate_ = "08/05/20"

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import sin,pi,cos
import numpy as np

#CAS PENDULE##############################

t0=0
tf=10
theta_0=pi/4
theta_prime_0=0
g=9.81
l=1
n=5000

precision=10**(-10)
nb_iteration=10
Y_0=[theta_0,theta_prime_0]


def pendule_et_Jac(Y,t):
    theta=Y[0]
    theta_prime=Y[1]
    if t==0:
        F=np.array([theta_0,theta_prime_0])
        JacF = np.array([[0., 1.], [(-g/l)*cos(theta_0), 0.]])
        return (F, JacF)
    F = np.array([theta_prime,(-g/l)*sin(theta)])
    JacF = np.array([[0., 1.], [(-g/l)*cos(theta), 0.]])
    return (F, JacF)

def pendule(Y,t):
    theta=Y[0]
    theta_prime=Y[1]
    if t==0:
        return np.array([theta_0,theta_prime_0])
    return np.array([theta_prime,(-g/l)*sin(theta)])

def newton_raphson_vectoriel(Y0,G_method,U_n,t,FONC):
    i=0
    X=Y0 #X est un vecteur 2x1 ici
    while (abs(np.dot((G_method(X,U_n,FONC,t)[0]).T,np.linalg.inv(G_method(X,U_n,FONC,t)[1]))).all()>=precision) and (i<=nb_iteration):
        X=(X-(np.dot((G_method(X,U_n,FONC,t)[0]).T,np.linalg.inv(G_method(X,U_n,FONC,t)[1]))))
        i+=1
    return X

def Euler_implicite_pendule(Y0,t0,tf,n,FONC):
    tps=[]
    h=(tf-t0)/float(n)
    Y = np.zeros((2, n+1))
    Y[:, 0] = Y0

    def G_euler_imp(U,U_n,FONC,t):
        G=U-h*(FONC(U,t))[0]-U_n
        dG=np.identity(2, float)-h*(FONC(U,t))[1]
        return (G,dG)
    tn=0
    tps+=[tn]
    for n in range(n):
        tn = tn + h
        tps+=[tn]
        Y[:, n+1] = newton_raphson_vectoriel(Y[:,0],G_euler_imp,Y[:,n],tn,FONC)
    return (Y,tps)

def Runge_kutta_ordre4_pendule(Y0,t0,tf,n,FONC):
    tps=[]
    h=(tf-t0)/float(n)
    Y = np.zeros((2, n+1))
    Y[:, 0] = Y0

    def G_rk4(U,U_n,FONC,t):
        U_n1=U_n
        U_n2=U_n1+(h/2)*(FONC(U_n1,t))[0]
        U_n3=U_n1+(h/2)*(FONC(U_n2,t))[0]
        U_n4=U_n1+h*(FONC(U_n3,t))[0]

        G=U-U_n-(h/6)*(FONC(U_n1,t)[0]+2*FONC(U_n2,t)[0]+2*FONC(U_n3,t)[0]+FONC(U_n4,t)[0])
        dG=np.identity(2,float)-(h/6)*(FONC(U_n1,t)[1]+2*FONC(U_n2,t)[1]+2*FONC(U_n3,t)[1]+FONC(U_n4,t)[1])
        return (G, dG)
    tn=0
    tps+=[tn]
    for n in range(n):
        tn = tn + h
        tps+=[tn]
        Y[:, n+1] = newton_raphson_vectoriel(Y[:,0],G_rk4,Y[:,n],tn,FONC)
    return (Y,tps)

def Verlet_pendule(Y0,t0,tf,n):
    t=t0
    x0=Y0[0]
    y0=Y0[1]
    x=x0
    y=y0
    h=(tf-t0)/float(n)
    temps=[t0]
    X=[x0]
    Y=[y0]
    for i in range(n):
        t=t+h
        temps+=[t]
        x,y =x+h*y-((h*h*g)/(2*l))*sin(x+(h/2)*y),y-((h*g)/l)*sin(x+(h/2)*y)
        X+=[x]
        Y+=[y]

    return temps,X,Y

def energie(X,Y):
    tps=len(X)
    E=[0]*tps
    for i in range(tps):
        E[i]=(1/2)*l*l*Y[i]*Y[i]+g*l*(1-cos(X[i]))
    return E

Temps=[i*((tf-t0)/(n)) for i in range(n+1)]
pendule_odeint=odeint(pendule,Y_0,Temps)

(sol_imp,temps_imp)=Euler_implicite_pendule(Y_0,t0,tf,n,pendule_et_Jac)
(sol_rk4,temps_rk4)=Runge_kutta_ordre4_pendule(Y_0,t0,tf,n,pendule_et_Jac)
(temps_Verlet,sol_Verlet_x,sol_Verlet_y)=Verlet_pendule(Y_0,t0,tf,n)

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Euler implicite Pendule n=5000',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps,pendule_odeint[:,0:1],label='scipy angle')
plt.plot(Temps,pendule_odeint[:,1:2],label='scipy vitesse angulaire')
plt.plot(temps_imp,sol_imp[0],label='angle')
plt.plot(temps_imp,sol_imp[1],label='vitesse angulaire')
plt.plot(temps_imp,energie(sol_imp[0],sol_imp[1]), label='energie mecanique')
plt.legend(bbox_to_anchor=(0.7, 0.02), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_pendule_euler_implicite_n5000")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Runge Kutta ordre 4 Pendule n=5000',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps,pendule_odeint[:,0:1],label='scipy angle')
plt.plot(Temps,pendule_odeint[:,1:2],label='scipy vitesse angulaire')
plt.plot(temps_rk4,sol_rk4[0],label='angle')
plt.plot(temps_rk4,sol_rk4[1],label='vitesse angulaire')
plt.plot(temps_rk4,energie(sol_rk4[0],sol_rk4[1]), label='energie mecanique')
plt.legend(bbox_to_anchor=(0.7, 0.02), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_pendule_runge_kutta4_n5000")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Verlet Pendule n=5000',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps,pendule_odeint[:,0:1],label='scipy angle')
plt.plot(Temps,pendule_odeint[:,1:2],label='scipy vitesse angulaire')
plt.plot(temps_Verlet,sol_Verlet_x,label='angle')
plt.plot(temps_Verlet,sol_Verlet_y,label='vitesse angulaire')
plt.plot(temps_Verlet,energie(sol_Verlet_x,sol_Verlet_y), label='energie mecanique')
plt.legend(bbox_to_anchor=(0.7, 0.02), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_pendule_verlet_n5000")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Plan de phase Euler implicite Pendule n=5000',fontsize=18)
plt.xlabel('Angle',fontsize=14)
plt.ylabel('Vitesse angulaire',fontsize=14)
plt.plot(pendule_odeint[:,0:1],pendule_odeint[:,1:2],label='scipy')
plt.plot(sol_imp[0],sol_imp[1],label='E_imp')
plt.legend(bbox_to_anchor=(0.6, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("plan_de_phase_euler_imp_pendule_n5000")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Plan de phase Runge kutta 4 Pendule n=5000',fontsize=18)
plt.xlabel('Angle',fontsize=14)
plt.ylabel('Vitesse angulaire',fontsize=14)
plt.plot(pendule_odeint[:,0:1],pendule_odeint[:,1:2],label='scipy')
plt.plot(sol_rk4[0],sol_rk4[1],label='Rk4')
plt.legend(bbox_to_anchor=(0.6, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("plan_de_phase_rk4_pendule_n5000")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Plan de phase Verlet Pendule n=5000',fontsize=18)
plt.xlabel('Angle',fontsize=14)
plt.ylabel('Vitesse angulaire',fontsize=14)
plt.plot(pendule_odeint[:,0:1],pendule_odeint[:,1:2],label='scipy')
plt.plot(sol_Verlet_x,sol_Verlet_y,label='Verlet')
plt.legend(bbox_to_anchor=(0.6, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("plan_de_phase_verlet_pendule_n5000")


plt.show()
