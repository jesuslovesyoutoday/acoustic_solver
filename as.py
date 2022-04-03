import numpy as np
from matplotlib import pyplot as plt


#TODO: проверки вводимых значений, значения по умолчанию

def init_print(string):
    
    """

    Prints given string in the box
    
    Parameters:
        string (str): A string to print

    """
    print("\n")
    print("##############################################")
    print("#                                            #")
    print(string)
    print("#                                            #")
    print("##############################################")
    print("\n")


def calculate(u, C):

    """

    Calculates displacement of the wave

    Parameters:
        u (np.array): Matrix with initial conditions
        C    (float): Courant number (C < 1)

    Returns:
        u (np.array): Matrix of displacement

    """
    # Calculating n = 1
    for m in range(1, M-1):
        u[1, m] = (u[0, m] 
        - 0.5 * C**2 * (u[0, m+1] - 2*u[0, m] + u[0, m-1]))

    # Calculating n = 2,...,N 
    for n in range(1, N-1):
        for m in range(1, M-1):
            u[n+1, m] = (-u[n-1, m] +2*u[n, m] 
            + C**2 * (u[n, m+1] - 2*u[n, m] + u[n,m-1]))
    return(u)

def graph(u, dx):
    
    """

    Draws displacement of the wave

    Parameters:
        u (np.array): Matrix of displacement
        dx   (float): Mesh step
    
    """
    x = []
    for i in range(M):
        x.append(i*dx)
    i = 0
    for frame in u:
        y = frame
        fig, axs = plt.subplots()
        axs.plot(x, y)
        plt.ylim(-0.1, 0.1)
        plt.grid()
        fig.savefig(str(i) + ".png")
        i += 1
    

init_print("#            Enter Freuency: (Hz)            #")
f = float(input())

init_print("#         Enter period of time: (s)          #")
T = float(input())

init_print("#          Enter path length: (m)            #")
L = float(input())

init_print("#            Enter time step: (s)            #")
dt = float(input())

# C <= 1 - ?
init_print("#                  Enter CFL:                #")
C = float(input())

n = 1
c = 343

# Wave length
l = c/(n*f)

# Phase velocity
v = l*f

# Mesh step
dx = (v*dt)/C

# Amount of time steps
N = int(T/dt)

# Amount of space steps
M = int(L/dx)

# Displacement function
u = np.zeros((N, M), dtype=np.float)

# Init conditions
for i in range(M):
    u[0, i] = i/100

U = calculate(u, C)
print(U)
graph(U, dx)
