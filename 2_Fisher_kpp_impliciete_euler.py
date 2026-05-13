import numpy as np
import matplotlib.pyplot as plt

def Fisher_KPP_vergelijking_impliciet(D,r,L,N,T,t):
    dx = L/(N-1)        #de lengte van een discretisatiestap
    dt = T/t            #t is het aantal tijdstappen
    x = np.linspace(0,L,N)
    u = 1/(1+np.exp(x-10))

    term = (dt * D) / dx**2

    waarde_nevendiagonaal = -term

    for k in range(t):              
        hoofddiagonaal = 1 + 2*term - dt*r*(1-u)
        A = np.diag(hoofddiagonaal) + np.diag(waarde_nevendiagonaal*np.ones(N-1),1) + np.diag(waarde_nevendiagonaal*np.ones(N-1),-1)
        A[0,1] *= 2                          
        A[N-1,N-2] *= 2

        u_new = np.linalg.solve(A,u)
        u = u_new.copy()

        if k % max(1, t // 100) == 0 or k == t - 1:
            plt.clf()
            plt.plot(x, u_new)
            plt.ylim(0, 1.1)
            plt.xlim(0, L+0.01)
            plt.xlabel("x")
            plt.ylabel("u(x,t)")
            plt.title(f"Fisher-KPP impliciet, tijdstap = {k}")
            plt.pause(0.00011)

        #afbeeldingen opslaan
        if k == 0:
            plt.savefig("impliciet_t=0.png")
        if k == 85:
            plt.savefig("impliciet_t=85.png")
        if k == 170:
            plt.savefig("impliciet_t=170.png")
        if k == 260:
            plt.savefig("impliciet_t=240.png")
# print(Fisher_KPP_vergelijking_impliciet(1, 3, 20, 50, 10, 500))
# print(Fisher_KPP_vergelijking_impliciet(1, 1, 20, 200, 32, 7000))
