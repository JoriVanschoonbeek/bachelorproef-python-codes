import numpy as np
import matplotlib.pyplot as plt

def Fisher_kpp_2D(D, r, L, N, T, t):
    h = L / (N - 1)        
    dt = T / t

    # moet groter zijn dan 0.25
    stabiliteitsvoorwaarde = D * dt / (h**2)

    #opstellen van het grid
    x = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, x, indexing="ij")
    
    # beginvoorwaarde
    u = -np.exp((-((X - L/2)**2 + (Y - L/2)**2))/10) +1

    # Turn on interactive mode for smooth live plotting
    plt.ion()
    fig = plt.figure()

    # 3. Vectorized Time Loop
    for z in range(t):
        #toevoegen van kopie 2de waarde voor de randvoorwaarde, en voorlaatste waarde na de 2e randvoorwaarde
        u_pad = np.pad(u, pad_width=1, mode='reflect')
        
        #berekenen laplaciaan
        laplacian = (u_pad[:-2, 1:-1] + u_pad[2:, 1:-1] + 
                     u_pad[1:-1, :-2] + u_pad[1:-1, 2:] - 
                     4 * u) / h**2
        
        # Fisher-KPP update stap
        u = u + dt * (D * laplacian + r * u * (1 - u))
        
        if z % max(1, t // 50) == 0 or z == t - 1:
            plt.clf()
            fixed_levels = np.linspace(0, 1, 21)
            plt.contourf(X, Y, u, levels=fixed_levels, cmap='viridis')
            plt.colorbar(label="u(x,y)")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.title(f"Fisher-KPP 2D, tijdstap = {z}")
            plt.pause(0.01)
        if z == 0:
            plt.savefig("explicie_2D_t=0.png")
        if z == 666:
            plt.savefig("expliciet_2D_t=666.png")
        if z == 1333:
            plt.savefig("expliciet_2D_t=1333.png")
        if z == 2400:
            plt.savefig("expliciet_2D_t=2400.png")
            
    plt.ioff() 
    plt.show()

# Example usage (ensure D*dt/h^2 <= 0.25):
Fisher_kpp_2D(D=0.1, r=1.0, L=10, N=50, T=32, t=13000)
