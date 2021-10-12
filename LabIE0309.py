#Lab IE0309

import numpy as np

class Node:
    def __init__(self): 

        self.voltage = 0
class Impedance:
    def __init__(self, value, node_a,node_b):
        self.value = value
        self.node_a = node_a
        self.node_b = node_b
    def get_current(self):
        return (self.node_a.voltage-self.node_b.voltage)/self.value
    
class V_source:

    def __init__(self, value, node_plus, node_minus):

        self.value = value
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.current = 0

    def get_current(self):
        return self.current

class I_source:
    def __init__(self, value, node_top, node_tail):
       self.value = value
       self.node_top = node_top
       self.node_tail = node_tail

    def get_current(self):
        return self.value

class Circuit:

    def __init__(self):

        self.nodes = []
        self.impedances = []
        self.v_sources = []
        self.i_sources = []

    def add_node(self):

        n = Node()
        self.nodes.append(n)

        return n
    
    def add_impedance(self, value, node_a, node_b):
        z = Impedance(value, node_a, node_b)
        self.impedances.append(z)
        
        return z

    def add_v_source(self, value, node_plus, node_minus):
        vs = V_source(value, node_plus, node_minus)
        self.v_sources.append(vs)

        return vs

    def add_i_source(self, value, node_top, node_tail):
        Is = I_source(value, node_top, node_tail)
        self.i_sources.append(Is)
        
        return Is

    def node2ind(self, node):

        return self.nodes.index(node) - 1

    def solve(self):
        
        N = len(self.nodes) #Numero de nodos
        F = len(self.v_sources) #Numero de fuentes
        I = len(self.i_sources) #Numero de fuentes de corriente
        
        #Build b
        b = np.zeros([N - 1 + F, 1], dtype= np.complex_)

        if F > 0:
            b[-F:, 0] = [source.value for source in self.v_sources]

        #if I > 0:
            for (i, Is) in enumerate(self.i_sources):
                i = self.node2ind(Is.node_top)
                j = self.node2ind(Is.node_tail)
                if i !=-1:
                    b[i,0] += Is.value
                if j != -1:
                    b[j,0] -= Is.value

        #Build A
        A = np.zeros([N -1 + F, N-1 +F], dtype = np.complex_)

        # Add terms due to impedances
        for z in self.impedances:
            i = self.node2ind(z.node_a)
            j= self.node2ind(z.node_b)
            if i == -1:
                A[j, j] += 1.0/z.value
            elif j == -1:
                A[i, i]+= 1.0/z.value
            else:
                A[i, i] += 1.0/z.value
                A[j, j] += 1.0/z.value
                A[i,j] -= 1.0/z.value
                A[j,i] -= 1.0/z.value

        #Add terms due to voltages sources
        for (k, vs) in enumerate(self.v_sources):
            i = self.node2ind(vs.node_plus)
            j= self.node2ind(vs.node_minus)
            if i != -1:
                A[i, N - 1 + k] = -1 #KCL
                A[N - 1 + k, i] = 1 #KVL
            if j != -1:
                A[j, N - 1 + k] = 1 #KCL
                A[N - 1 + k, j] = -1 #KVL
                
        # Solve
        x = np.linalg.solve(A,b)

        #Save results
        for(i,n) in enumerate(self.nodes[1:]):
            n.voltage = x[i,0]

        for(i, Is) in enumerate(x[-F:,0]):
            self.v_sources[i].current=Is

    def measure_v(self, node_plus, node_minus):
        return node_plus.voltage - node_minus.voltage

    def measure_i(self, X):
        return X.get_current()
    
def pol2rect(mag, degs):
    rads = np.deg2rad(degs)
    real = mag*np.cos(rads)
    imag = mag*np.sin(rads)
    return real +1j*imag

def rect2pol(z):
    rads = np.angle(z)
    degs = np.angle(z)
    degs = np.rad2deg(rads)
    mag = np.abs(z)
    return mag, degs

c = Circuit()

#Nodes
GND = c.add_node()
n1= c.add_node()
n2 = c.add_node()
n3 = c.add_node()

#Sources
v_1 = c.add_v_source(pol2rect(40, 90), n2, GND)
I_1 = c.add_i_source(pol2rect(3,0), n1, n3)

#Impedances
R1 = 5
R2 = 8
R3 = 10
R4 = 20
L1 = 1j*4 
L2 = 1j*15
C = -1j*2

z1 = c.add_impedance(R1, n1, n2)
z2 = c.add_impedance(R2+C, n2, n3)
z3 = c.add_impedance(R3+L1, n3, GND)
z4 = c.add_impedance(R4+L2, n1, GND)

c.solve()
Io = c.measure_i(z4)
Io_mag, Io_ang = rect2pol(Io)
print('Magnitud', round(Io_mag, 2))
print('Fase', round(Io_ang, 2))






























    
        
