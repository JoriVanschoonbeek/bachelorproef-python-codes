import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import block_diag

def matrix_2D(D, r, L, N, T ,t):
    dt = T/t
    dx = L/(N-1)

    a = -4*(D*dt)/(dx**2) + 1
    b = r * dt

    vector_a = np.full(N, a)
    vector_b = np.full(N-1, b)

    A = np.zeros((N,N)) + np.diag(vector_a) + np.diag(vector_b,-1) + np.diag(vector_b,1)
    A[0][1] *= 2
    A[N-1][N-2] *= 2
    
    M = block_diag(*([A] * N))
    vector_2b_b = N*[2*b] + (N**2-2*N)*[b]
    vector_b_2b = (N**2-2*N)*[b] + N*[2*b]
    M = M + np.diag(vector_2b_b,N) + np.diag(vector_b_2b,-N)
    return M
print(matrix_2D(D=0.1, r=1.0, L=10, N=50, T=32, t=13000))

