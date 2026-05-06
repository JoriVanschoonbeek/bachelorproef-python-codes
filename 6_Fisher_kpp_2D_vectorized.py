import numpy as np
import matplotlib.pyplot as plt

def Fisher_kpp_2D(D, r, L, N, T, t):
    h = L / (N - 1)        
    dt = T / t

    # 1. Stability Check
    cfl = D * dt / (h**2)
    if cfl > 0.25:
        print(f"Warning: Explicit method may be unstable. D*dt/h^2 = {cfl:.3f} (Should be <= 0.25)")

    # 2. Grid Setup
    x = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, x, indexing="ij")
    
    # Initial condition
    u = -np.exp((-((X - L/2)**2 + (Y - L/2)**2))/10) +1

    # Turn on interactive mode for smooth live plotting
    plt.ion()
    fig = plt.figure()

    # 3. Vectorized Time Loop
    for z in range(t):
        # Pad the array to elegantly handle zero-flux (Neumann) boundary conditions.
        # 'reflect' mirrors the edge values, acting exactly like ghost points.
        u_pad = np.pad(u, pad_width=1, mode='reflect')
        
        # Calculate the 2D Laplacian using array slicing (Central Difference)
        laplacian = (u_pad[:-2, 1:-1] + u_pad[2:, 1:-1] + 
                     u_pad[1:-1, :-2] + u_pad[1:-1, 2:] - 
                     4 * u) / h**2
        
        # Fisher-KPP update step
        u = u + dt * (D * laplacian + r * u * (1 - u))
        
        # Plotting (Update every few steps to avoid plotting bottleneck)
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