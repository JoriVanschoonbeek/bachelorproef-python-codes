import numpy as np
import matplotlib.pyplot as plt

def Fisher_KPP_vergelijking_expliciet(D, r, L, N, T, t):
    #T is de volledige tijd
    #L is de lengte
    #N is het aantal discretisatiepunten
    #dt is de lengte van een tijdstap

    dx = L/(N-1)                     #de lengte van een discretisatiestap
    dt = T/t                         #t is het aantal tijdstappen

    x = np.linspace(0,L,N)           #opstellen van een net 
    u = 1/(1+np.exp(x-10))

    stabiliteitsvoorwaarde = min((dx**2)/(2*D), 2/r)
    u_new = np.zeros(N)

    term = D*dt/(dx**2)

    for k in range(t):
        for i in range(1,N-1):
            u_new[i] = dt * (D * ((u[i-1] - 2*u[i] + u[i+1])/(dx**2)) + r*u[i]*(1-u[i])) + u[i]
        u_new[0] = dt * (D * ((2*u[1] - 2*u[0])/(dx**2)) + r*u[0]*(1-u[0])) + u[0]
        u_new[N-1] = dt * (D * ((2*u[N-2] - 2*u[N-1])/(dx**2))  + r*u[N-1]*(1-u[N-1])) + u[N-1]

        u = u_new.copy()
        u_new = np.zeros(N)
     
        if k % max(1, t // 100) == 0 or k == t - 1:
            plt.clf()
            plt.plot(x, u)
            plt.ylim(0, 1.1)
            plt.xlim(0, L+0.01)
            plt.xlabel("x")
            plt.ylabel("u(x,t)")
            plt.title(f"Fisher-KPP expliciet, tijdstap = {k}")
            plt.draw()
            plt.pause(0.00011)


        #afbeeldingen opslaan indien inputs niet voldoen aan de stabiliteitsvoorwaarde
        if dt > stabiliteitsvoorwaarde:
            if k == 0:
                plt.savefig("expliciet_onvoldaan_svw_t=0.png")
            if k == 6:
                plt.savefig("expliciet_onvoldaan_svw_t=6.png")
            if k == 12:
                plt.savefig("expliciet_onvoldaan_svw_t=12.png")
            if k == 15:
                plt.savefig("expliciet_onvoldaan_svw_t=15.png")
        
        #afbeeldingen opslaan indien inputs niet voldoen aan de stabiliteitsvoorwaarde
        if dt < stabiliteitsvoorwaarde:
            if k == 0:
                plt.savefig("expliciet_voldaan_svw_t=0.png")
            if k == 500:
                plt.savefig("expliciet_voldaan_svw_t=500.png")
            if k == 1000:
                plt.savefig("expliciet_voldaan_svw_t=1000.png")
            if k == 3000:
                plt.savefig("expliciet_voldaan_svw_t=3000.png")

#voldoet aan stabiliteitsvoorwaarde
#print(Fisher_KPP_vergelijking_expliciet(1,1,20,200,32,7000))

#voldoet niet aan stabilitetisvoorwaarde
#print(Fisher_KPP_vergelijking_expliciet(D=1,r=1,L=20,N=100,T=10, t=50 )) # ⇒ dt = 0.2 (veel te groot!)


#randvoorwaarden nakijken!!!!!
