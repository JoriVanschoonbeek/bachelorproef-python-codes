import numpy as np
import matplotlib.pyplot as plt

def Fisher_KPP_vergelijking_met_snelheid_expliciet(D, r, v, L, N, T, t):
    #T is de volledige tijd
    #L is de lengte
    #N is het aantal discretisatiepunten
    #dt is de lengte van een tijdstap
    #v is de advectie

    dx = L/(N-1)        #de lengte van een discretisatiestap
    dt = T/t            #t is het aantal tijdstappen

    x = np.linspace(0,L,N)
    u = 1/(1+np.exp(x-10))

    u_new = np.zeros(N)

    stabiliteitsvoorwaarde = (2*D*dt)/dx**2 + v*dt/dx 

    term = D*dt/(dx**2)

    for k in range(t):
        for i in range(1,N-1):
            u_new[i] = dt * (D * ((u[i-1] - 2*u[i] + u[i+1])/(dx**2)) + r*u[i]*(1-u[i]) - (v/dx)*(u[i] - u[i-1])) + u[i]
        u_new[0] = dt * (D * ((2*u[1] - 2*u[0])/(dx**2)) + r*u[0]*(1-u[0]) - (v/dx)*(u[0] - u[1])) + u[0]
        u_new[N-1] = dt * (D * ((2*u[N-2] - 2*u[N-1])/(dx**2))  + r*u[N-1]*(1-u[N-1]) - (v/dx)*(u[N-1] - u[N-2])) + u[N-1]

        u = u_new.copy()
        u_new = np.zeros(N)
     
        if k % max(1, t // 100) == 0 or k == t - 1:
            plt.clf()
            plt.plot(x, u)
            plt.ylim(0, 1.1)
            plt.xlim(0, L+0.01)
            plt.xlabel("x")
            plt.ylabel("u(x,t)")
            plt.title(f"Fisher-KPP expliciet met advectieterm, tijdstap = {k}")
            plt.pause(0.001)

        #afbeeldingen opslaan voor inputwaarden die voldoen aan stabiliteitsvoorwaarde
        if stabiliteitsvoorwaarde < 1:
            if k == 0:
                plt.savefig("expliciet_advectie_voldaan_svw_t=0.png")
            if k == 750:
                plt.savefig("expliciet_advectie_voldaan_svw_t=750.png")
            if k == 1500:
                plt.savefig("expliciet_advectie_voldaan_svw_t=1500.png")
            if k == 3000:
                plt.savefig("expliciet_advectie_voldaan_svw_t=3000.png")
        #afbeeldingen opslaan voor inputwaarden die niet voldoen aan stabiliteitsvoorwaarde
        else:
            if k == 0:
                plt.savefig("expliciet_advectie_onvoldaan_svw_t=0.png")
            if k == 7:
                plt.savefig("expliciet_advectie_onvoldaan_svw_t=7.png")
            if k == 10:
                plt.savefig("expliciet_advectie_onvoldaan_svw_t=10.png")
            if k == 15:
                plt.savefig("expliciet_advectie_onvoldaan_svw_t=15.png")
#print(Fisher_KPP_vergelijking_met_snelheid_expliciet(1,1,2,20,200,10,5000))
#print(Fisher_KPP_vergelijking_met_snelheid_expliciet(1, 2, 5, 20, 100, 5, 50))

