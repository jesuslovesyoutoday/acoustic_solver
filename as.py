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
    
class Wave():
    
    def __init__(self, f, T, L, dt, C, n):
    
        """
        Initiation of parameters of the wave
        
        Parameters:
            f  (float) - frequency (Hz)
            T  (float) - time period (s)
            L  (float) - path's length (m)
            dt (float) - time step (s)
            C  (float) - courant number (<1)
            n  (float) - refractive index
        """
    
        self.f  = f
        self.T  = T
        self.L  = L
        self.dt = dt
        self.C  = C
        self.n0  = n
        
        self.c  = 343             # Sound speed
        
        self.l  = self.c/(n*f)    # Wave length

        self.v  = self.c/n        # Phase velocity

        self.dx = (self.v*dt)/C   # Mesh step

        self.N  = int(T/dt)       # Amount of time steps

        self.M  = int(L/self.dx)  # Amount of space steps
        
        # Refractive index
        self.n = np.zeros((self.N, self.M), dtype=np.float)
        
        # Displacement function
        self.u  = np.zeros((self.N, self.M), dtype=np.float)
        # Displacement function for analitical solution  
        self.u1 = np.zeros((self.N, self.M), dtype=np.float)  
                                                              
        self.A  = np.amax(self.u) # Amplitude of the wave
        
    def set_n(self):
    
        """
        Specifies the distribution of the refractive index
        
        """
        for n in range(self.N):
            for m in range(self.M):
                self.n[n, m] = 1 + m/100  
        
    def init_conditions(self):
        
        """
        Sets initial conditions of the wave
            
        """
        for i in range(self.M):
            self.u[0, i] = math.sin(self.dx*i*math.pi/self.L)
                
                                                    
    def calculate(self):
        
        """
        Calculates displacement of the wave
            
        """
            
        # Calculating n = 1
        for m in range(1, self.M-1):
            self.u[1, m] = (self.u[0, m] 
            - 0.5 * self.C**2 * (self.u[0, m+1] 
            - 2*self.u[0, m] + self.u[0, m-1]))

        # Calculating n = 2,...,N 
        for n in range(1, self.N-1):
            for m in range(1, self.M-1):
                self.u[n+1, m] = (-self.u[n-1, m] +2*self.u[n, m] 
              + self.C**2 * (self.u[n, m+1] - 2*self.u[n, m] 
              + self.u[n,m-1]))
                    
    def snapshots(self, file_name):

        self.structuredGrid = vtk.vtkStructuredGrid()
        self.points = vtk.vtkPoints()
        self.U = vtk.vtkDoubleArray()
        self.U.SetNumberOfComponents(3)
        self.U.SetName("u")
        self.grid = np.mgrid[0:self.M, 0:self.N]
            
        for m in range(self.M):
            for n in range(self.N):
                self.points.InsertNextPoint(self.grid[0][m, n], 
                                            self.grid[1][m, n], 0)
                self.U.InsertNextTuple((0, self.u[n,m], 0))
        self.structuredGrid.SetDimensions(self.M, self.N, 1)
        self.structuredGrid.SetPoints(self.points)
        self.structuredGrid.GetPointData().AddArray(self.U)
        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetInputDataObject(self.structuredGrid)
        writer.SetFileName(file_name)
        writer.Write()
            
    def test(self):
    
        """
        Calculates analytical solution (with no initial
        conditions)
            
        """
        self.A  = np.amax(self.u)
        for n in range(1, self.N-1):
            for m in range(1, self.M-1):
                x = m*self.dx
                t = n*self.dt
                self.u1[n, m] = (self.A * math.sin((math.pi * x)/self.L) 
                        * math.cos((math.pi * self.v * t)/self.L))
        self.u1[0, self.M - 1] = 0                  


def graph(u, u1, dx, dt, M, N):
    
    """
    Draws displacement of the wave

    Parameters:
        u (np.array): First matrix
        u1(np.array): Second matrix
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
    for i in range(0, N, 30):
        y = u[i]
        #y1 = u1[i, 0 : len(x)]
        y1 = u1[i]
        fig, axs = plt.subplots()
        axs.plot(x, y,'-', x, y1, '--')
        #axs.plot(x, y)
        plt.ylim(-np.amax(u), np.amax(u))
        axs.set_title("t = " + '%.4f'%t[i])
        plt.grid()
        fig.savefig("anim_new/" + str(i) + ".png")                    
                    

def main():

    #input_values()

    # Default values
    f = 2000        # Frequency
    T = 10          # Time period
    L = 3430        # Path's length
    dt = 0.001       # Time step
    C = 0.1         # CFL < 1
    """n = 1           # Refractive index

    u_1 = Wave(f, T, L, dt, C, n)

    u_1.init_conditions()
    u_1.calculate()
    u_1.snapshots("wave_1.vts")
    u_1.test()"""
    
    n = 1
    u_2 = Wave(f, T, L, dt, C, n)

    u_2.init_conditions()
    u_2.calculate()
    u_2.test()
    #u_2.snapshots("wave_n1_finer_mesh.vts")
    
    
    graph(u_2.u, u_2.u1, u_2.dx, dt, u_2.M, u_2.N)

if __name__ == '__main__':
    main()                  
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
