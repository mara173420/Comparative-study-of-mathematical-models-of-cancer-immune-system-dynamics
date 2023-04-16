
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from scipy import io
from scipy.special import logsumexp


def rungeKutta(ODE, init, dt):
    #Inputs:
    #f: function we want to integrate
    #init: list with all the needed parameters (V, m, n, h)
    #I: applied current
    #h: step
    k1 = dt * ODE(init)
    val_k2 = init + k1/2
    k2 = dt * ODE(val_k2)
    val_k3 = init + k2/2
    k3 = dt * ODE(val_k3)
    val_k4 = init + k3
    k4 = dt * ODE(val_k4)
    
    k = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0
    
    init = init + k
    
    #Return the values of the chosen variables for each call of the function
    return init

#Variables

def drug_therapy(init):  
 T = init[0]
 L = init[1]
 C = init[2]
 I2 = init[3]
 
 a=0.13 #tumor growth rate
 b=2.3e-9 #reciprocal carrying capacity
 c=4.4e-9 #tumor cell kill rate by immune cells
 d=7.3e6 #CTL immune cell induction rate
 e=9.9e-9 #CTL proliferation rate induced by IL-2
 f=0.33 #CTL immune cell death rate
 g=1.6e7 #antigen presentation
 j=3.3e-9 #rate of consumption of IL-2 by CTL
 k=1.8e-8 #inactivation of IL-2 molecules
 l=3e6 #half-saturation constant
 Mt=0.9 #tumor cell kill rate via chemotherapy
 Ml=0.6 #CTL immune cell kill rate via chemotherapy
 p=6.4 #decay rate of chemotherapy
 Vc = 1 #chemotherapy induction rate
 Vi = 0 #immunotherapy induction rate

 dT = a*T*(1-b*T)-c*T*L-Mt*(1-np.exp(-C))*T
 dL = d + e*L*I2-f*L-Ml*(1-np.exp(-C))*L
 dC = Vc - p*C
 dI2 = Vi + ((g*T)/(T+l))-j*L*I2-k*T*I2
 
 
 return np.array([dT, dL, dC, dI2])

def combined(init):  
 T = init[0]
 L = init[1]
 C = init[2]
 I2 = init[3]
 
 a=0.13 #tumor growth rate
 b=2.3e-9 #reciprocal carrying capacity
 c=4.4e-9 #tumor cell kill rate by immune cells
 d=7.3e6 #CTL immune cell induction rate
 e=9.9e-9 #CTL proliferation rate induced by IL-2
 f=0.33 #CTL immune cell death rate
 g=1.6e7 #antigen presentation
 j=3.3e-9 #rate of consumption of IL-2 by CTL
 k=1.8e-8 #inactivation of IL-2 molecules
 l=3e6 #half-saturation constant
 Mt=0.9 #tumor cell kill rate via chemotherapy
 Ml=0.6 #CTL immune cell kill rate via chemotherapy
 p=6.4 #decay rate of chemotherapy
 Vc = 1 #chemotherapy induction rate
 Vi = 10e6 #immunotherapy induction rate

 dT = a*T*(1-b*T)-c*T*L-Mt*(1-np.exp(-C))*T
 dL = d + e*L*I2-f*L-Ml*(1-np.exp(-C))*L
 dC = Vc - p*C
 dI2 = Vi + ((g*T)/(T+l))-j*L*I2-k*T*I2
 
 
 return np.array([dT, dL, dC, dI2])

def immuno(init):  
 T = init[0]
 L = init[1]
 #C = init[2]
 I2 = init[2]
 
 a=0.13 #tumor growth rate
 b=2.3e-9 #reciprocal carrying capacity
 c=4.4e-9 #tumor cell kill rate by immune cells
 d=7.3e6 #CTL immune cell induction rate
 e=9.9e-9 #CTL proliferation rate induced by IL-2
 f=0.33 #CTL immune cell death rate
 g=1.6e7 #antigen presentation
 j=3.3e-9 #rate of consumption of IL-2 by CTL
 k=1.8e-8 #inactivation of IL-2 molecules
 l=3e6 #half-saturation constant
 Mt=0.9 #tumor cell kill rate via chemotherapy
 Ml=0.6 #CTL immune cell kill rate via chemotherapy
 p=6.4 #decay rate of chemotherapy
 Vc = 0 #chemotherapy induction rate
 Vi = 10e6 #immunotherapy induction rate

 dT = a*T*(1-b*T)-c*T*L
 dL = d +e*L*I2-f*L
 #dC = Vc - p*C
 dI2 = Vi + ((g*T)/(T+l))-j*L*I2-k*T*I2
 
 return np.array([dT, dL, dI2])

def growth(init):  
 T =  init
 
 a = 0.13
 b = 2.3e-9
 
 dT = a*T*(1-b*T)

 return dT

#Initial conditions
init = 3e7 #initial tumor

dt = 1
tvec = np.arange(0, 250, dt)

#Empty lists for the variables 

T_g = np.zeros(len(tvec))

for i, t in enumerate(tvec):
  init = rungeKutta(growth, init, dt) #Integrate

  T_g[i] = init #Store variables in empty lists


#Initial conditions

init = [3e7, 2.25e7, 0, 2.4e7] #initial values

dt = 0.1
tvec = np.arange(0, 250, dt)

#Empty lists for the variables 

T_drug = np.zeros(len(tvec))
L = np.zeros(len(tvec))
C = np.zeros(len(tvec))
I2 = np.zeros(len(tvec))


for i, t in enumerate(tvec):
  init = rungeKutta(drug_therapy, init, dt) #Integrate

  T_drug[i], L[i], C[i], I2[i] = init #Store variables in empty lists
  
init = [3e7, 2.25e7, 0, 2.4e7] #initial values

dt = 0.1
tvec = np.arange(0, 250, dt)

#Empty lists for the variables 

T_com = np.zeros(len(tvec))
L = np.zeros(len(tvec))
C = np.zeros(len(tvec))
I2 = np.zeros(len(tvec))


for i, t in enumerate(tvec):
  init = rungeKutta(combined, init, dt) #Integrate

  T_com[i], L[i], C[i], I2[i] = init #Store variables in empty lists
  
init = [3e7, 2.25e7, 2.4e7] #initial values

dt = 0.1
tvec = np.arange(0, 250, dt)

#Empty lists for the variables 

T_imm = np.zeros(len(tvec))
L = np.zeros(len(tvec))
#C = np.zeros(len(tvec))
I2 = np.zeros(len(tvec))


for i, t in enumerate(tvec):
  init = rungeKutta(immuno, init, dt) #Integrate

  T_imm[i], L[i], I2[i] = init #Store variables in empty lists  
  

 
fig = plt.figure()

plt.plot(T_imm)
plt.plot(T_drug)
plt.plot(T_com)
plt.legend(["Immunotherapy", "Drug Therapy", "Combined Therapy"])
plt.draw()
plt.title("Cancer Treatment Therapies")
plt.xlim([0, 250])
plt.xlabel("Time (t), days")
plt.ylabel("Tumor Cell Concentrations")
plt.show()




