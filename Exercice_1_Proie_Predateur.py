_author_ = "Nicolas Coruzzi"
_filename_ = "Exercice_1_Proie_Predateur"
_creationdate_ = "07/05/20"

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import atan
from copy import deepcopy

#CAS PROIE_PREDATEUR############################

t0=0
tf=20

alpha=0.25
beta=-0.01
gamma=-1
delta=0.01
x0=80
y0=30

def proie_predateur(Y,t):
    x=Y[0]
    y=Y[1]
    return [alpha*x+beta*x*y,gamma*y+delta*x*y]

def Euler_explicite_proie_predateur(f,t0,tf,x0,y0,n) :
    t=t0
    x=x0
    y=y0
    h=(tf-t0)/float(n)
    temps=[t0]
    proie=[x0]
    predateur=[y0]
    for i in range(n):
        t=t+h
        temps+=[t]
        F,G=f((x,y),t)
        x,y=x+h*F,y+h*G
        proie+=[x]
        predateur+=[y]

    return temps,proie,predateur

def Runge_kutta_ordre2_proie_predateur(f,t0,tf,x0,y0,n):
    t=t0
    x=x0
    y=y0
    h=(tf-t0)/float(n)
    temps=[t0]
    proie=[x0]
    predateur=[y0]
    for i in range(n):
        t=t+h
        temps+=[t]
        F,G=f((x+(h/2)*f((x,y),t)[0],y+(h/2)*f((x,y),t)[1]),t+(h/2))
        x,y =x+h*F,y+h*G
        proie+=[x]
        predateur+=[y]

    return temps,proie,predateur

def Runge_kutta_ordre4_proie_predateur(f,t0,tf,x0,y0,n):
    t=t0
    x=x0
    y=y0
    h=(tf-t0)/float(n)
    temps=[t0]
    proie=[x0]
    predateur=[y0]
    for i in range(n):
        t=t+h
        temps+=[t]
        k1=f((x,y),t)
        k2=f((x+(h*k1[0])/2,y+(h*k1[1])/2),t+(h/2))
        k3=f((x+(h*k2[0])/2,y+(h*k2[1])/2),t+(h/2))
        k4=f((x+(h*k3[0]),y+(h*k3[1])),t+(h/2))
        x,y =x+(h/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),y+(h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        proie+=[x]
        predateur+=[y]

    return temps,proie,predateur

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
proie_predateur_odeint1=odeint(proie_predateur,[x0,y0],Temps1)

Temps2=[i*((tf-t0)/(n2)) for i in range(n2+1)]
proie_predateur_odeint2=odeint(proie_predateur,[x0,y0],Temps2)

Temps3=[i*((tf-t0)/(n3)) for i in range(n3+1)]
proie_predateur_odeint3=odeint(proie_predateur,[x0,y0],Temps3)

Temps4=[i*((tf-t0)/(n4)) for i in range(n4+1)]
proie_predateur_odeint4=odeint(proie_predateur,[x0,y0],Temps4)


proie_predateur_euler_explicite1=Euler_explicite_proie_predateur(proie_predateur,t0,tf,x0,y0,n1)
proie_predateur_euler_explicite2=Euler_explicite_proie_predateur(proie_predateur,t0,tf,x0,y0,n2)
proie_predateur_euler_explicite3=Euler_explicite_proie_predateur(proie_predateur,t0,tf,x0,y0,n3)
proie_predateur_euler_explicite4=Euler_explicite_proie_predateur(proie_predateur,t0,tf,x0,y0,n4)

proie_predateur_Runge_kutta_ordre2_1=Runge_kutta_ordre2_proie_predateur(proie_predateur,t0,tf,x0,y0,n1)
proie_predateur_Runge_kutta_ordre2_2=Runge_kutta_ordre2_proie_predateur(proie_predateur,t0,tf,x0,y0,n2)
proie_predateur_Runge_kutta_ordre2_3=Runge_kutta_ordre2_proie_predateur(proie_predateur,t0,tf,x0,y0,n3)
proie_predateur_Runge_kutta_ordre2_4=Runge_kutta_ordre2_proie_predateur(proie_predateur,t0,tf,x0,y0,n4)

proie_predateur_Runge_kutta_ordre4_1=Runge_kutta_ordre4_proie_predateur(proie_predateur,t0,tf,x0,y0,n1)
proie_predateur_Runge_kutta_ordre4_2=Runge_kutta_ordre4_proie_predateur(proie_predateur,t0,tf,x0,y0,n2)
proie_predateur_Runge_kutta_ordre4_3=Runge_kutta_ordre4_proie_predateur(proie_predateur,t0,tf,x0,y0,n3)
proie_predateur_Runge_kutta_ordre4_4=Runge_kutta_ordre4_proie_predateur(proie_predateur,t0,tf,x0,y0,n4)

plt.figure(figsize=(14,10), dpi=80)
plt.title('Proie-Predateur n=20',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps1,proie_predateur_odeint1[:,0:1],label='proie scipy')
plt.plot(Temps1,proie_predateur_odeint1[:,1:2],label='predateur scipy')
plt.plot(proie_predateur_euler_explicite1[0],proie_predateur_euler_explicite1[1],".-",label='proie Euler')
plt.plot(proie_predateur_euler_explicite1[0],proie_predateur_euler_explicite1[2],".-",label='predateur Euler')
plt.plot(proie_predateur_Runge_kutta_ordre2_1[0],proie_predateur_Runge_kutta_ordre2_1[1],"+-",label='proie Rk2')
plt.plot(proie_predateur_Runge_kutta_ordre2_1[0],proie_predateur_Runge_kutta_ordre2_1[2],"+-",label='predateur Rk2')
plt.plot(proie_predateur_Runge_kutta_ordre4_1[0],proie_predateur_Runge_kutta_ordre4_1[1],"x-",label='proie Rk4')
plt.plot(proie_predateur_Runge_kutta_ordre4_1[0],proie_predateur_Runge_kutta_ordre4_1[2],"x-",label='predateur Rk4')
plt.legend(bbox_to_anchor=(0.01, 0.68), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_proie_predateur_n20")


plt.figure(figsize=(14,10), dpi=80)
plt.title('Proie-Predateur n=30',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps2,proie_predateur_odeint2[:,0:1],label='proie scipy')
plt.plot(Temps2,proie_predateur_odeint2[:,1:2],label='predateur scipy')
plt.plot(proie_predateur_euler_explicite2[0],proie_predateur_euler_explicite2[1],".-",label='proie Euler')
plt.plot(proie_predateur_euler_explicite2[0],proie_predateur_euler_explicite2[2],".-",label='predateur Euler')
plt.plot(proie_predateur_Runge_kutta_ordre2_2[0],proie_predateur_Runge_kutta_ordre2_2[1],"+-",label='proie Rk2')
plt.plot(proie_predateur_Runge_kutta_ordre2_2[0],proie_predateur_Runge_kutta_ordre2_2[2],"+-",label='predateur Rk2')
plt.plot(proie_predateur_Runge_kutta_ordre4_2[0],proie_predateur_Runge_kutta_ordre4_2[1],"x-",label='proie Rk4')
plt.plot(proie_predateur_Runge_kutta_ordre4_2[0],proie_predateur_Runge_kutta_ordre4_2[2],"x-",label='predateur Rk4')
plt.legend(bbox_to_anchor=(0.01, 0.68), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_proie_predateur_n30")


plt.figure(figsize=(14,10), dpi=80)
plt.title('Proie-Predateur n=100',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps3,proie_predateur_odeint3[:,0:1],label='proie scipy')
plt.plot(Temps3,proie_predateur_odeint3[:,1:2],label='predateur scipy')
plt.plot(proie_predateur_euler_explicite3[0],proie_predateur_euler_explicite3[1],".-",label='proie Euler')
plt.plot(proie_predateur_euler_explicite3[0],proie_predateur_euler_explicite3[2],".-",label='predateur Euler')
plt.plot(proie_predateur_Runge_kutta_ordre2_3[0],proie_predateur_Runge_kutta_ordre2_3[1],"+-",label='proie Rk2')
plt.plot(proie_predateur_Runge_kutta_ordre2_3[0],proie_predateur_Runge_kutta_ordre2_3[2],"+-",label='predateur Rk2')
plt.plot(proie_predateur_Runge_kutta_ordre4_3[0],proie_predateur_Runge_kutta_ordre4_3[1],"x-",label='proie Rk4')
plt.plot(proie_predateur_Runge_kutta_ordre4_3[0],proie_predateur_Runge_kutta_ordre4_3[2],"x-",label='predateur Rk4')
plt.legend(bbox_to_anchor=(0.01, 0.68), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_proie_predateur_n100")


plt.figure(figsize=(14,10), dpi=80)
plt.title('Proie-Predateur n=400',fontsize=18)
plt.xlabel('Temps',fontsize=14)
plt.ylabel('Solution',fontsize=14)
plt.plot(Temps4,proie_predateur_odeint4[:,0:1],label='proie scipy')
plt.plot(Temps4,proie_predateur_odeint4[:,1:2],label='predateur scipy')
plt.plot(proie_predateur_euler_explicite4[0],proie_predateur_euler_explicite4[1],".-",label='proie Euler')
plt.plot(proie_predateur_euler_explicite4[0],proie_predateur_euler_explicite4[2],".-",label='predateur Euler')
plt.plot(proie_predateur_Runge_kutta_ordre2_4[0],proie_predateur_Runge_kutta_ordre2_4[1],"+-",label='proie Rk2')
plt.plot(proie_predateur_Runge_kutta_ordre2_4[0],proie_predateur_Runge_kutta_ordre2_4[2],"+-",label='predateur Rk2')
plt.plot(proie_predateur_Runge_kutta_ordre4_4[0],proie_predateur_Runge_kutta_ordre4_4[1],"x-",label='proie Rk4')
plt.plot(proie_predateur_Runge_kutta_ordre4_4[0],proie_predateur_Runge_kutta_ordre4_4[2],"x-",label='predateur Rk4')
plt.legend(bbox_to_anchor=(0.01, 0.68), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("comparaison_proie_predateur_n400")

#plt.show()

#CALCUL DE L'ERREUR#######################

err_eul_exp_proie=[]
err_eul_exp_predateur=[]
del_eul_exp=[]
erreur_delta_euler_explicite_n1_proie=erreur(proie_predateur_odeint1[:,0:1],proie_predateur_euler_explicite1[1],tf,n1)
erreur_delta_euler_explicite_n1_predateur=erreur(proie_predateur_odeint1[:,1:2],proie_predateur_euler_explicite1[2],tf,n1)
err_eul_exp_proie+=[erreur_delta_euler_explicite_n1_proie[1]]
err_eul_exp_predateur+=[erreur_delta_euler_explicite_n1_predateur[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n1_proie[0]]
erreur_delta_euler_explicite_n2_proie=erreur(proie_predateur_odeint2[:,0:1],proie_predateur_euler_explicite2[1],tf,n2)
erreur_delta_euler_explicite_n2_predateur=erreur(proie_predateur_odeint2[:,1:2],proie_predateur_euler_explicite2[2],tf,n2)
err_eul_exp_proie+=[erreur_delta_euler_explicite_n2_proie[1]]
err_eul_exp_predateur+=[erreur_delta_euler_explicite_n2_predateur[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n2_proie[0]]
erreur_delta_euler_explicite_n3_proie=erreur(proie_predateur_odeint3[:,0:1],proie_predateur_euler_explicite3[1],tf,n3)
erreur_delta_euler_explicite_n3_predateur=erreur(proie_predateur_odeint3[:,1:2],proie_predateur_euler_explicite3[2],tf,n3)
err_eul_exp_proie+=[erreur_delta_euler_explicite_n3_proie[1]]
err_eul_exp_predateur+=[erreur_delta_euler_explicite_n3_predateur[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n3_proie[0]]
erreur_delta_euler_explicite_n4_proie=erreur(proie_predateur_odeint4[:,0:1],proie_predateur_euler_explicite4[1],tf,n4)
erreur_delta_euler_explicite_n4_predateur=erreur(proie_predateur_odeint4[:,1:2],proie_predateur_euler_explicite4[2],tf,n4)
err_eul_exp_proie+=[erreur_delta_euler_explicite_n4_proie[1]]
err_eul_exp_predateur+=[erreur_delta_euler_explicite_n4_predateur[1]]
del_eul_exp+=[erreur_delta_euler_explicite_n4_proie[0]]



err_rk2_proie=[]
err_rk2_predateur=[]
del_rk2=[]
erreur_delta_runge_ordre2_n1_proie=erreur(proie_predateur_odeint1[:,0:1],proie_predateur_Runge_kutta_ordre2_1[1],tf,n1)
erreur_delta_runge_ordre2_n1_predateur=erreur(proie_predateur_odeint1[:,1:2],proie_predateur_Runge_kutta_ordre2_1[2],tf,n1)
err_rk2_proie+=[erreur_delta_runge_ordre2_n1_proie[1]]
err_rk2_predateur+=[erreur_delta_runge_ordre2_n1_predateur[1]]
del_rk2+=[erreur_delta_runge_ordre2_n1_proie[0]]
erreur_delta_runge_ordre2_n2_proie=erreur(proie_predateur_odeint2[:,0:1],proie_predateur_Runge_kutta_ordre2_2[1],tf,n2)
erreur_delta_runge_ordre2_n2_predateur=erreur(proie_predateur_odeint2[:,1:2],proie_predateur_Runge_kutta_ordre2_2[2],tf,n2)
err_rk2_proie+=[erreur_delta_runge_ordre2_n2_proie[1]]
err_rk2_predateur+=[erreur_delta_runge_ordre2_n2_predateur[1]]
del_rk2+=[erreur_delta_runge_ordre2_n2_proie[0]]
erreur_delta_runge_ordre2_n3_proie=erreur(proie_predateur_odeint3[:,0:1],proie_predateur_Runge_kutta_ordre2_3[1],tf,n3)
erreur_delta_runge_ordre2_n3_predateur=erreur(proie_predateur_odeint3[:,1:2],proie_predateur_Runge_kutta_ordre2_3[2],tf,n3)
err_rk2_proie+=[erreur_delta_runge_ordre2_n3_proie[1]]
err_rk2_predateur+=[erreur_delta_runge_ordre2_n3_predateur[1]]
del_rk2+=[erreur_delta_runge_ordre2_n3_proie[0]]
erreur_delta_runge_ordre2_n4_proie=erreur(proie_predateur_odeint4[:,0:1],proie_predateur_Runge_kutta_ordre2_4[1],tf,n4)
erreur_delta_runge_ordre2_n4_predateur=erreur(proie_predateur_odeint4[:,1:2],proie_predateur_Runge_kutta_ordre2_4[2],tf,n4)
err_rk2_proie+=[erreur_delta_runge_ordre2_n4_proie[1]]
err_rk2_predateur+=[erreur_delta_runge_ordre2_n4_predateur[1]]
del_rk2+=[erreur_delta_runge_ordre2_n4_proie[0]]


err_rk4_proie=[]
err_rk4_predateur=[]
del_rk4=[]
erreur_delta_runge_ordre4_n1_proie=erreur(proie_predateur_odeint1[:,0:1],proie_predateur_Runge_kutta_ordre4_1[1],tf,n1)
erreur_delta_runge_ordre4_n1_predateur=erreur(proie_predateur_odeint1[:,1:2],proie_predateur_Runge_kutta_ordre4_1[2],tf,n1)
err_rk4_proie+=[erreur_delta_runge_ordre4_n1_proie[1]]
err_rk4_predateur+=[erreur_delta_runge_ordre4_n1_predateur[1]]
del_rk4+=[erreur_delta_runge_ordre4_n1_proie[0]]
erreur_delta_runge_ordre4_n2_proie=erreur(proie_predateur_odeint2[:,0:1],proie_predateur_Runge_kutta_ordre4_2[1],tf,n2)
erreur_delta_runge_ordre4_n2_predateur=erreur(proie_predateur_odeint2[:,1:2],proie_predateur_Runge_kutta_ordre4_2[2],tf,n2)
err_rk4_proie+=[erreur_delta_runge_ordre4_n2_proie[1]]
err_rk4_predateur+=[erreur_delta_runge_ordre4_n2_predateur[1]]
del_rk4+=[erreur_delta_runge_ordre4_n2_proie[0]]
erreur_delta_runge_ordre4_n3_proie=erreur(proie_predateur_odeint3[:,0:1],proie_predateur_Runge_kutta_ordre4_3[1],tf,n3)
erreur_delta_runge_ordre4_n3_predateur=erreur(proie_predateur_odeint3[:,1:2],proie_predateur_Runge_kutta_ordre4_3[2],tf,n3)
err_rk4_proie+=[erreur_delta_runge_ordre4_n3_proie[1]]
err_rk4_predateur+=[erreur_delta_runge_ordre4_n3_predateur[1]]
del_rk4+=[erreur_delta_runge_ordre4_n3_proie[0]]
erreur_delta_runge_ordre4_n4_proie=erreur(proie_predateur_odeint4[:,0:1],proie_predateur_Runge_kutta_ordre4_4[1],tf,n4)
erreur_delta_runge_ordre4_n4_predateur=erreur(proie_predateur_odeint4[:,1:2],proie_predateur_Runge_kutta_ordre4_4[2],tf,n4)
err_rk4_proie+=[erreur_delta_runge_ordre4_n4_proie[1]]
err_rk4_predateur+=[erreur_delta_runge_ordre4_n4_predateur[1]]
del_rk4+=[erreur_delta_runge_ordre4_n4_proie[0]]


pente_euler_exp=calcul_pente_gauche(err_eul_exp_proie,del_eul_exp)
pente_rk2=calcul_pente_gauche(err_rk2_proie,del_rk2)
pente_rk4=calcul_pente_gauche(err_rk4_proie,del_rk4)

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Erreur proie-predateur Euler',fontsize=18)
plt.annotate(' pente: '+ str(pente_euler_exp),
             fontsize=14, xy = (0.5,0.3),
             xycoords='figure fraction', xytext = (0.5,0.3),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))
plt.xlabel('Delta',fontsize=14)
plt.ylabel('Erreur',fontsize=14)
plt.plot(del_eul_exp,err_eul_exp_proie,label='proie')
plt.plot(del_eul_exp,err_eul_exp_predateur,label='predateur')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0.01, 0.68), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("erreur_proie_predateur_euler")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Erreur proie-predateur Rk2',fontsize=18)
plt.annotate(' pente: '+ str(pente_rk2),
             fontsize=14, xy = (0.5,0.3),
             xycoords='figure fraction', xytext = (0.5,0.3),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))
plt.xlabel('Delta',fontsize=14)
plt.ylabel('Erreur',fontsize=14)
plt.plot(del_rk2,err_rk2_proie,label='proie')
plt.plot(del_rk2,err_rk2_predateur,label='predateur')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0.01, 0.68), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("erreur_proie_predateur_rk2")

plt.figure(figsize=(9.5,6.5), dpi=80)
plt.title('Erreur proie-predateur Rk4',fontsize=18)
plt.annotate(' pente: '+ str(pente_rk4),
             fontsize=14, xy = (0.5,0.3),
             xycoords='figure fraction', xytext = (0.5,0.3),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))
plt.xlabel('Delta',fontsize=14)
plt.ylabel('Erreur',fontsize=14)
plt.plot(del_rk4,err_rk4_proie,label='proie')
plt.plot(del_rk4,err_rk4_predateur,label='predateur')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0.01, 0.68), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.savefig("erreur_proie_predateur_rk4")

plt.show()