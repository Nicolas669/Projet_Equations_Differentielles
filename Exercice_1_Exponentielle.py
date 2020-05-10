_author_ = "Nicolas Coruzzi"
_filename_ = "Exercice_1_Exponentielle"
_creationdate_ = "07/05/20"

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import atan
from copy import deepcopy

#CAS EXPONENTIELLE##############################

t0=0
tf=20
#Temps=[i for i in range(t0,tf+1)]
y0=1

def exponentielle(Y,t):
    if t==0:
        return y0
    return Y

def Euler_explicite_exponentielle(f,t0,tf,y0,n):
    t=t0
    y=y0
    h=(tf-t0)/float(n)
    temps=[t0]
    resultat=[y0]
    for i in range(n):
        t=t+h
        temps+=[t]
        F=f(y,t)
        y=y+h*F
        resultat+=[y]

    return temps,resultat

def Runge_kutta_ordre2_exponentielle(f,t0,tf,y0,n):
    t=t0
    y=y0
    h=(tf-t0)/float(n)
    temps=[t0]
    resultat=[y0]
    for i in range(n):
        t=t+h
        temps+=[t]
        F=f(y+(h/2)*f(y,t),t+(h/2))
        y=y+h*F
        resultat+=[y]

    return temps,resultat

def Runge_kutta_ordre4_exponentielle(f,t0,tf,y0,n):
    t=t0
    y=y0
    h=(tf-t0)/float(n)
    temps=[t0]
    resultat=[y0]
    for i in range(n):
        t=t+h
        temps+=[t]
        k1=f(y,t)
        k2=f(y+(h*k1)/2,t+(h/2))
        k3=f(y+(h*k2)/2,t+(h/2))
        k4=f(y+h*k3,t+h)
        y=y+(h/6)*(k1+2*k2+2*k3+k4)
        resultat+=[y]

    return temps,resultat

def erreur(u,v,t_tot,nb):
    #u est avec scipy
    #v est avec une methode num
    delta_t=t_tot/nb
    res=0
    for i in range (len(u)):
        res+=(u[i]-v[i])**2
    #on renvoie l'erreur et le delta_t
    return ((delta_t*res)**(1/2),delta_t)

def calcul_pente_gauche(abscisses,ordonnees):
    #s est les ordonnees
    #t est les abscisses
    s=deepcopy(abscisses)
    t=deepcopy(ordonnees)
    for i in range(len(s)):
        if s[i]==min(s) and t[i]==min(t):
            a=t[i]
            fa=s[i]
            del s[i]
            del t[i]
            for j in range(len(s)):
                if s[j]==min(s) and t[j]==min(t):
                    b=t[j]
                    fb=s[j]
                    return atan((fb-fa)/(b-a))

n1=20
n2=30
n3=100
n4=400

Temps1=[i*((tf-t0)/(n1)) for i in range(n1+1)]
exponentielle_odeint1=odeint(exponentielle,y0,Temps1)

Temps2=[i*((tf-t0)/(n2)) for i in range(n2+1)]
exponentielle_odeint2=odeint(exponentielle,y0,Temps2)

Temps3=[i*((tf-t0)/(n3)) for i in range(n3+1)]
exponentielle_odeint3=odeint(exponentielle,y0,Temps3)

Temps4=[i*((tf-t0)/(n4)) for i in range(n4+1)]
exponentielle_odeint4=odeint(exponentielle,y0,Temps4)

exponentielle_euler_explicite1=Euler_explicite_exponentielle(exponentielle,t0,tf,y0,n1)
exponentielle_euler_explicite2=Euler_explicite_exponentielle(exponentielle,t0,tf,y0,n2)
exponentielle_euler_explicite3=Euler_explicite_exponentielle(exponentielle,t0,tf,y0,n3)
exponentielle_euler_explicite4=Euler_explicite_exponentielle(exponentielle,t0,tf,y0,n4)

exponentielle_runge_ordre2_1=Runge_kutta_ordre2_exponentielle(exponentielle,t0,tf,y0,n1)
exponentielle_runge_ordre2_2=Runge_kutta_ordre2_exponentielle(exponentielle,t0,tf,y0,n2)
exponentielle_runge_ordre2_3=Runge_kutta_ordre2_exponentielle(exponentielle,t0,tf,y0,n3)
exponentielle_runge_ordre2_4=Runge_kutta_ordre2_exponentielle(exponentielle,t0,tf,y0,n4)

exponentielle_runge_ordre4_1=Runge_kutta_ordre4_exponentielle(exponentielle,t0,tf,y0,n1)
exponentielle_runge_ordre4_2=Runge_kutta_ordre4_exponentielle(exponentielle,t0,tf,y0,n2)
exponentielle_runge_ordre4_3=Runge_kutta_ordre4_exponentielle(exponentielle,t0,tf,y0,n3)
exponentielle_runge_ordre4_4=Runge_kutta_ordre4_exponentielle(exponentielle,t0,tf,y0,n4)


plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Exponentielle n=20',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps1,exponentielle_odeint1,label='scipy')
plt.plot(exponentielle_euler_explicite1[0],exponentielle_euler_explicite1[1],".-",label='Euler')
plt.plot(exponentielle_runge_ordre2_1[0],exponentielle_runge_ordre2_1[1],"+-",label='Rk2')
plt.plot(exponentielle_runge_ordre4_1[0],exponentielle_runge_ordre4_1[1],"x-",label='Rk4')
plt.legend(bbox_to_anchor=(0.01, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_exponentielle_n20")


plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Exponentielle n=30',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps2,exponentielle_odeint2,label='scipy')
plt.plot(exponentielle_euler_explicite2[0],exponentielle_euler_explicite2[1],".-",label='Euler')
plt.plot(exponentielle_runge_ordre2_2[0],exponentielle_runge_ordre2_2[1],"+-",label='Rk2')
plt.plot(exponentielle_runge_ordre4_2[0],exponentielle_runge_ordre4_2[1],"x-",label='Rk4')
plt.legend(bbox_to_anchor=(0.01, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_exponentielle_n30")


plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Exponentielle n=100',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps3,exponentielle_odeint3,label='scipy')
plt.plot(exponentielle_euler_explicite3[0],exponentielle_euler_explicite3[1],".-",label='Euler')
plt.plot(exponentielle_runge_ordre2_3[0],exponentielle_runge_ordre2_3[1],"+-",label='Rk2')
plt.plot(exponentielle_runge_ordre4_3[0],exponentielle_runge_ordre4_3[1],"x-",label='Rk4')
plt.legend(bbox_to_anchor=(0.01, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_exponentielle_n100")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Exponentielle n=400',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps4,exponentielle_odeint4,label='scipy')
plt.plot(exponentielle_euler_explicite4[0],exponentielle_euler_explicite4[1],".-",label='Euler')
plt.plot(exponentielle_runge_ordre2_4[0],exponentielle_runge_ordre2_4[1],"+-",label='Rk2')
plt.plot(exponentielle_runge_ordre4_4[0],exponentielle_runge_ordre4_4[1],"x-",label='Rk4')
plt.legend(bbox_to_anchor=(0.01, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_exponentielle_n400")

#plt.show()

#CALCUL DE L'ERREUR#######################

err_eul_exp=[]
del_eul_exp=[]
erreur_delta_euler_explicite_n1=erreur(exponentielle_odeint1,exponentielle_euler_explicite1[1],tf,n1)
err_eul_exp+=[erreur_delta_euler_explicite_n1[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n1[0]]
erreur_delta_euler_explicite_n2=erreur(exponentielle_odeint2,exponentielle_euler_explicite2[1],tf,n2)
err_eul_exp+=[erreur_delta_euler_explicite_n2[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n2[0]]
erreur_delta_euler_explicite_n3=erreur(exponentielle_odeint3,exponentielle_euler_explicite3[1],tf,n3)
err_eul_exp+=[erreur_delta_euler_explicite_n3[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n3[0]]
erreur_delta_euler_explicite_n4=erreur(exponentielle_odeint4,exponentielle_euler_explicite4[1],tf,n4)
err_eul_exp+=[erreur_delta_euler_explicite_n4[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n4[0]]

err_rk2=[]
del_rk2=[]
erreur_delta_runge_ordre2_n1=erreur(exponentielle_odeint1,exponentielle_runge_ordre2_1[1],tf,n1)
err_rk2+=[erreur_delta_runge_ordre2_n1[1]]
del_rk2+=[erreur_delta_runge_ordre2_n1[0]]
erreur_delta_runge_ordre2_n2=erreur(exponentielle_odeint2,exponentielle_runge_ordre2_2[1],tf,n2)
err_rk2+=[erreur_delta_runge_ordre2_n2[1]]
del_rk2+=[erreur_delta_runge_ordre2_n2[0]]
erreur_delta_runge_ordre2_n3=erreur(exponentielle_odeint3,exponentielle_runge_ordre2_3[1],tf,n3)
err_rk2+=[erreur_delta_runge_ordre2_n3[1]]
del_rk2+=[erreur_delta_runge_ordre2_n3[0]]
erreur_delta_runge_ordre2_n4=erreur(exponentielle_odeint4,exponentielle_runge_ordre2_4[1],tf,n4)
err_rk2+=[erreur_delta_runge_ordre2_n4[1]]
del_rk2+=[erreur_delta_runge_ordre2_n4[0]]

err_rk4=[]
del_rk4=[]
erreur_delta_runge_ordre4_n1=erreur(exponentielle_odeint1,exponentielle_runge_ordre4_1[1],tf,n1)
err_rk4+=[erreur_delta_runge_ordre4_n1[1]]
del_rk4+=[erreur_delta_runge_ordre4_n1[0]]
erreur_delta_runge_ordre4_n2=erreur(exponentielle_odeint2,exponentielle_runge_ordre4_2[1],tf,n2)
err_rk4+=[erreur_delta_runge_ordre4_n2[1]]
del_rk4+=[erreur_delta_runge_ordre4_n2[0]]
erreur_delta_runge_ordre4_n3=erreur(exponentielle_odeint3,exponentielle_runge_ordre4_3[1],tf,n3)
err_rk4+=[erreur_delta_runge_ordre4_n3[1]]
del_rk4+=[erreur_delta_runge_ordre4_n3[0]]
erreur_delta_runge_ordre4_n4=erreur(exponentielle_odeint4,exponentielle_runge_ordre4_4[1],tf,n4)
err_rk4+=[erreur_delta_runge_ordre4_n4[1]]
del_rk4+=[erreur_delta_runge_ordre4_n4[0]]

pente_euler_exp=calcul_pente_gauche(err_eul_exp,del_eul_exp)
pente_rk2=calcul_pente_gauche(err_rk2,del_rk2)
pente_rk4=calcul_pente_gauche(err_rk4,del_rk4)

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Erreur exponentielle Euler',fontsize=18)
plt.annotate(' pente: '+ str(pente_euler_exp),
             fontsize=14, xy = (0.5,0.3),
             xycoords='figure fraction', xytext = (0.5,0.3),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))
plt.xlabel('Delta',fontsize=14)
plt.ylabel('Erreur',fontsize=14)
plt.plot(del_eul_exp,err_eul_exp)
plt.xscale('log')
plt.yscale('log')
plt.savefig("erreur_exponentielle_euler")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Erreur exponentielle Rk2',fontsize=18)
plt.annotate(' pente: '+ str(pente_rk2),
             fontsize=14, xy = (0.5,0.3),
             xycoords='figure fraction', xytext = (0.5,0.3),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))
plt.xlabel('Delta',fontsize=14)
plt.ylabel('Erreur',fontsize=14)
plt.plot(del_rk2,err_rk2)
plt.xscale('log')
plt.yscale('log')
plt.savefig("erreur_exponentielle_rk2")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Erreur exponentielle Rk4',fontsize=18)
plt.annotate(' pente: '+ str(pente_rk4),
             fontsize=14, xy = (0.5,0.3),
             xycoords='figure fraction', xytext = (0.5,0.3),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))
plt.xlabel('Delta',fontsize=14)
plt.ylabel('Erreur',fontsize=14)
plt.plot(del_rk4,err_rk4)
plt.xscale('log')
plt.yscale('log')
plt.savefig("erreur_exponentielle_rk4")

plt.show()