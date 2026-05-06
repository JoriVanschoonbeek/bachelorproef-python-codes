import numpy as np
import matplotlib.pyplot as plt

def Fisher_kpp_snelheid_impliciet(D,r,v,L,N,T,t):
    dx = L/(N-1)        #de lengte van een discretisatiestap
    dt = T/t            #t is het aantal tijdstappen
    x = np.linspace(0,L,N)
    u = 1/(1+np.exp(x-10))

    diffusie = (dt * D) / dx**2
    advectie = (dt * v) / dx
    
    waarde_boven_nevendiagonaal = -diffusie
    waarde_onder_nevendiagonaal = -diffusie - advectie

    for k in range(t):              
        hoofddiagonaal = 1 + 2*diffusie + advectie - dt*r*(1-u)
        A = np.diag(hoofddiagonaal) + np.diag(waarde_boven_nevendiagonaal*np.ones(N-1),1) + np.diag(waarde_onder_nevendiagonaal*np.ones(N-1),-1)
        A[0,1] = A[0,1] - diffusie - advectie                      
        A[N-1,N-2] = A[N-1,N-2] - diffusie

        u_new = np.linalg.solve(A,u)
        u = u_new.copy()

        plt.clf()
        plt.plot(x, u_new)
        plt.ylim(0, 1.1)
        plt.xlim(0, L+0.01)
        plt.xlabel("x")
        plt.ylabel("u(x,t)")
        plt.title(f"Fisher-KPP impliciet, tijdstap = {k}")
        plt.draw()
        plt.pause(0.00001)

        #afbeelding opslaan
        if k == 0:
            plt.savefig("impliciet_met_advectie_t=0.png")
        if k == 83:
            plt.savefig("impliciet_met_advectie_t=83.png")
        if k == 120:
            plt.savefig("impliciet_met_advectie_t=120.png")
        if k == 249:
            plt.savefig("impliciet_met_advectie_t=249.png")
            
#print(Fisher_kpp_snelheid_impliciet(1, 1, 0.20, 20, 50, 10, 250))
