import vtk
import math
import numpy as np
from matplotlib import pyplot as plt


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

def input_values():

    init_print("#            Enter Freuency: (Hz)            #")
    f = float(input())

    init_print("#         Enter period of time: (s)          #")
    T = float(input())

    init_print("#          Enter path length: (m)            #")
    L = float(input())

    init_print("#            Enter time step: (s)            #")
    dt = float(input())

    # C <= 1 
    init_print("#                  Enter CFL:                #")
    C = float(input())


def calculate(u, C, M, N):

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

def graph(u, u1, dx, dt, M, N):
    
    """
    Draws displacement of the wave

    Parameters:
        u (np.array): Matrix of displacement
        dx   (float): Mesh step
        dt   (float): Time step
        M      (int): Amount of space steps
        N      (int): Amount of time steps
    """
    x = []
    t = []
    for i in range(M):
        x.append(i*dx)
    for i in range(N):
        t.append(i*dt)
    i = 0
    for i in range(0, N, 3):
        y = u[i]
        y1 = u1[i]
        fig, axs = plt.subplots()
        axs.plot(x, y,'-', x, y1, '--')
        #axs.plot(x, y)
        plt.ylim(-np.amax(u), np.amax(u))
        axs.set_title("t = " + '%.4f'%t[i])
        plt.grid()
        fig.savefig("anim/" + str(i) + ".png")

def snapshots(u, dx, M, N):
    k = 1
    structuredGrid = vtk.vtkStructuredGrid()
    points = vtk.vtkPoints()
    #U = vtk.vtkDoubleArray()
    #U.SetName("u")
    points_1 = vtk.vtkPoints()
    du = np.amax(u)/N

    for n in range(N):
        for m in range(M):
            points.InsertNextPoint(dx * m, du * n, 0)
            points_1.InsertNextPoint(dx * m, u[n, m], 0)
            #U.InsertNextValue(u[n, m] * 1000)
        structuredGrid.SetDimensions(N, M, 1)
        structuredGrid.SetPoints(points)
        #structuredGrid.GetPointData().AddArray(U)
        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetInputDataObject(structuredGrid)
        writer.SetFileName("wave-step-" + str(n) + ".vts")
        writer.Write()

def test(u, A, L, v, dt, dx, N, M):
    
    """
    Calculates analytical solution (with no initial
    conditions)

    Parameters:
        u (np.array): Empty matrix of displacement
        A    (float): Amplitude of wave
        L    (float): Path's length
        v    (float): Phase velocity
        dt   (float): Time step
        dx   (float): Space step
    
    Returns:
        u (np.array): Filled matrix of displacement
    """
    for n in range(1, N-1):
        for m in range(1, M-1):
            x = m*dx
            t = n*dt
            u[n, m] = (A * math.sin((math.pi * x)/L) 
                    * math.cos((math.pi * v * t)/L))
    return(u)


def test_1(u, dx, dt, L, N, M):
    for n in range(0, N-1):
        for m in range(0, M-1):
            x = m*dx
            t = n*dt
            u[n, m] = x * (L - x) * math.sin(t)
    return(u)

def main():

    #input_values()

    # Default values
    f = 1000        # Frequency
    T = 10          # Time period
    L = 3430        # Path length
    dt = 0.01        # Time step
    C = 0.1           # CFL < 1
    n = 1           # Refractive index
    c = 343         # Sound speed

    l = c/(n*f)     # Wave length

    v = l*f         # Phase velocity

    dx = (v*dt)/C   # Mesh step

    N = int(T/dt)   # Amount of time steps

    M = int(L/dx)   # Amount of space steps

    u = np.zeros((N, M), dtype=np.float) # Displacement function
    u1 = np.zeros((N, M), dtype=np.float)


    for i in range(M):                   # Init conditions
        u[0, i] = math.sin(dx*i*math.pi/L)
        #u[0, i] = i/100
        #u1[0, i] = i/100
    U = calculate(u, C, M, N)

    A = np.amax(U)                       # Amplitude of the wave


    U1 = test(u1, A, L, v, dt, dx, N, M)
    U1[0, M-1] = 0
   

    #U1 = test_1(u1, dx, dt, L, N, M)
    #for n in range(N):
    #    u1[n, 0] = 0

    graph(U, U1, dx, dt, M, N)

    #snapshots(U, dx, M, N)
    
if __name__ == '__main__':
    main()
