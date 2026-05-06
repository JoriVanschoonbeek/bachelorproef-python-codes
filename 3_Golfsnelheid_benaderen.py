import numpy as np
import matplotlib.pyplot as plt

x = [0, 0.25, 0.5, 0.75, 1, 1.25]
u = [1, 0.9, 0.7, 0.4, 0.2, 0]
def find_wave_position(u_array, x_array):
    x_positie = 0
    for i in range(len(u_array)-1):
        if u_array[i] >= 0.5 and u_array[i+1] <= 0.5:
            u_left = u_array[i]
            u_right =u_array[i+1]
            x_left = x_array[i]
            x_right =x_array[i+1]
            x_positie = x_left + ((0.5 - u_left)/(u_right - u_left))*(x_right - x_left)
    if x_positie == 0:
        return 0
    return x_positie
#print(find_wave_position(u,x))


def Golfsnelheid_benaderen(D,r,L,N,T,t):
    dx = L/(N-1)        #de lengte van een discretisatiestap
    dt = T/t            #t is het aantal tijdstappen
    x = np.linspace(0,L,N)
    u = 1/(1+np.exp(x-5))

    term = (dt * D) / dx**2
    
    waarde_nevendiagonaal = -term
    x_positie = find_wave_position(u,x)

    speedlijst = []

    cnt = -1
    for k in range(t):  
        cnt += 1            
        hoofddiagonaal = 1 + 2*term - dt*r*(1-u)
        A = np.diag(hoofddiagonaal) + np.diag(waarde_nevendiagonaal*np.ones(N-1),1) + np.diag(waarde_nevendiagonaal*np.ones(N-1),-1)
        A[0,1] *= 2                          
        A[N-1,N-2] *= 2

        u_new = np.linalg.solve(A,u)
        u = u_new.copy()
        if cnt % 10 == 0:
            x_positie_new = find_wave_position(u,x)
            speed = (x_positie_new - x_positie)/(10*dt)
            x_positie = x_positie_new
            speedlijst.append(speed)
        
    speedlijst2 = []
    for i in speedlijst:
        if i > 0:
            speedlijst2.append(i)
    return round(np.mean(speedlijst2),2), round(2*np.sqrt(3),2)
#print(Golfsnelheid_benaderen(1, 3, 20, 50, 10, 500))
