_author_ = "Nicolas Coruzzi"
_filename_ = "Exercice_2_cas1"
_creationdate_ = "07/05/20"

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import atan
from copy import deepcopy

#CAS 1##############################

t0=0
tf=20
y0=1
#alpha doit etre strictement positif
alpha=3

def cas_1(Y,t):
    if t==0:
        return y0
    return -alpha*Y

def Euler_explicite_cas_1(f,t0,tf,y0,n):
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

def Euler_implicite_cas_1(t0,tf,y0,n):
    t=t0
    y=y0
    h=(tf-t0)/float(n)
    temps=[t0]
    resultat=[y0]
    for i in range(n):
        t=t+h
        temps+=[t]
        y=y/(1+h*alpha)
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

n3=100
n4=400

Temps3=[i*((tf-t0)/(n3)) for i in range(n3+1)]
cas_1_odeint3=odeint(cas_1,y0,Temps3)

Temps4=[i*((tf-t0)/(n4)) for i in range(n4+1)]
cas_1_odeint4=odeint(cas_1,y0,Temps4)

cas_1_euler_implicite3=Euler_implicite_cas_1(t0,tf,y0,n3)
cas_1_euler_implicite4=Euler_implicite_cas_1(t0,tf,y0,n4)

cas_1_euler_explicite3=Euler_explicite_cas_1(cas_1,t0,tf,y0,n3)
cas_1_euler_explicite4=Euler_explicite_cas_1(cas_1,t0,tf,y0,n4)


plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Exercice 2 Cas 1 n=100',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps3,cas_1_odeint3,label='scipy')
plt.plot(cas_1_euler_implicite3[0],cas_1_euler_implicite3[1],".-",label='E_imp')
plt.plot(cas_1_euler_explicite3[0],cas_1_euler_explicite3[1],"+-",label='E_exp')
plt.legend(bbox_to_anchor=(0.6, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_cas_1_n100")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Exercice 2 Cas 1 n=400',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps4,cas_1_odeint4,label='scipy')
plt.plot(cas_1_euler_implicite4[0],cas_1_euler_implicite4[1],".-",label='E_imp')
plt.plot(cas_1_euler_explicite4[0],cas_1_euler_explicite4[1],"+-",label='E_exp')
plt.legend(bbox_to_anchor=(0.6, 0.74), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_cas_1_n400")

#plt.show()

#CALCUL DE L'ERREUR#######################

err_eul_imp=[]
del_eul_imp=[]
erreur_delta_euler_implicite_n3=erreur(cas_1_odeint3,cas_1_euler_implicite3[1],tf,n3)
err_eul_imp+=[erreur_delta_euler_implicite_n3[1]]
del_eul_imp+=[erreur_delta_euler_implicite_n3[0]]
erreur_delta_euler_implicite_n4=erreur(cas_1_odeint4,cas_1_euler_implicite4[1],tf,n4)
err_eul_imp+=[erreur_delta_euler_implicite_n4[1]]
del_eul_imp+=[erreur_delta_euler_implicite_n4[0]]

err_eul_exp=[]
del_eul_exp=[]
erreur_delta_euler_explicite_n3=erreur(cas_1_odeint3,cas_1_euler_explicite3[1],tf,n3)
err_eul_exp+=[erreur_delta_euler_explicite_n3[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n3[0]]
erreur_delta_euler_explicite_n4=erreur(cas_1_odeint4,cas_1_euler_explicite4[1],tf,n4)
err_eul_exp+=[erreur_delta_euler_explicite_n4[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n4[0]]

pente_euler_exp=calcul_pente_gauche(err_eul_exp,del_eul_exp)
pente_euler_imp=calcul_pente_gauche(err_eul_imp,del_eul_imp)

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Erreur cas 1 Euler explicite',fontsize=18)
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
plt.savefig("erreur_cas_1_euler_exp")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Erreur cas 1 Euler implicite',fontsize=18)
plt.annotate(' pente: '+ str(pente_euler_imp),
             fontsize=14, xy = (0.5,0.3),
             xycoords='figure fraction', xytext = (0.5,0.3),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))
plt.xlabel('Delta',fontsize=14)
plt.ylabel('Erreur',fontsize=14)
plt.plot(del_eul_imp,err_eul_imp)
plt.xscale('log')
plt.yscale('log')
plt.savefig("erreur_cas_1_euler_imp")

plt.show()