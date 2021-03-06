import numpy as np
import matplotlib.pyplot as plt
import math

import os
here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)


# Subfunctions

# function Re = Re_t(eD)
# returns the Reynold number at fully developed turbulence
def Re_t(eD):

    phi = (-0.5/np.log10(eD/3.7))**2
    Re = 200/(eD*math.sqrt(phi))
    return Re

# function phi = blasius(Re)
# returns the friction factor according to Blasius (smooth pipe)
def blasius(Re):

    tol = 1e-5
    phi = 0  
    phi_new = 0.008 # initial value
    while(abs(phi-phi_new)>tol):
        phi = phi_new;
        phi_new = (1/(2.*np.log10(Re*math.sqrt(phi))-0.8))**2
    return phi

# function phi = cw(eD,Re)
# returns the friction factor according to Coolebrok and White
def cw(eD,Re):

    tol = 1e-5
    phi = 0    
    phi_new = 0.008 # initial value
    while(abs(phi-phi_new)>tol):
        phi = phi_new
        phi_new = (-0.5/np.log10(eD/3.7+2.51/(Re*math.sqrt(phi))))**2;
    return phi


# Plot Moody diagram

#plt.figure(figsize = (10,9))

Re1 = np.linspace(2000,200000,300)
Re2 = np.linspace(200000,1e8,300)
Re = np.concatenate((Re1,Re2))
m = len(Re)
eD = np.array([0.00001,0.00005,0.0001,0.0002,0.0004,0.0006,0.0008,0.001,0.002,0.004,0.006,0.008,0.01,0.015,0.02,0.03,0.04,0.05])
n = len(eD)
phi = np.zeros((m,n))
bla = np.zeros(m)
sat = np.zeros(n)
for i in range (n):
    for j in range (m):
        phi[j][i] = cw(eD[i],Re[j])
    sat[i] = Re_t(eD[i])

for i in range(m):
    bla[i] = blasius(Re[i])

Re_lam = np.linspace(500,2000,100)
phi_lam = 64/Re_lam
fig, ax = plt.subplots(1,1)
ax.loglog(Re,bla)

ax.plot(Re_lam,phi_lam)
ax.plot(sat,phi[-1,:])
ax.plot(Re,phi)
ax.axis([500, 1e8 ,0.008 ,0.1])
#grid on
edlab = np.chararray(n+3)
edlab = ['smooth pipe','laminar','saturation limit']

for i in range(n):
    edlab.append(str(eD[i]))
plt.legend(edlab,loc = (1,0))


####exeeeempleee ######"


Re_fonce = np.array([ 5242.75106656,  9642.91714027, 13855.84210447, 18724.11095199,
       21345.48648527])

phi_fonce_blas = np.ones(len(Re_fonce))
for i in range(len(Re_fonce)):
    phi_fonce_blas[i] = blasius(Re_fonce[i])

phi_fonce = np.array([0.03732542, 0.02961576, 0.02700069, 0.02464263, 0.023463727])

phi_clair = np.array([-0.05411661,  0.01599677,  0.0271176 ,  0.01909234,  0.02285257])
Re_clair = np.array([ 5242.75106656,  9642.91714027, 13855.84210447, 18724.11095199,
       21345.48648527])


plt.plot(Re_fonce,phi_fonce,"+k")
#plt.plot(Re_clair, phi_clair, "+k")







plt.xlabel('Re [-]')
plt.ylabel(r'$\Phi$ [-]', rotation = 0)


#plt.savefig('./images/moody_fonce.pdf')

plt.show()
















    
