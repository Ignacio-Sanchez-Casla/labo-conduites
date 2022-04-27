# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:29:24 2022

@author: Nacho
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

import os
here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)

# Plot parameters
save = True
plt.rc('axes', axisbelow=True)
def params():
    plt.figure(figsize=(7,5))
    plt.grid(True, linestyle='-.')
   
    
rho = 1e3 # [kg/m^3]
nu = 1e-6 # [??]
g = 9.81 # [m/s^2]

L = 0.914 # [m]

D_fonce = 13.6e-3 # [m]
A_fonce = np.pi * D_fonce**2 / 4 #[m^2]

D_clair = 26.2e-3 # [m]
A_clair = np.pi * D_clair**2 / 4 #[m^2]


# Donnees experimentales

debit = np.array([0.056, 0.103, 0.148, 0.2, 0.228]) / 1e3

z1 = np.array([603, 600, 594, 572, 550])
z2 = np.array([572, 521, 441, 312, 213])

z3 = np.array([410, 408, 405, 396, 554])
z4 = np.array([391, 357, 309, 236, 222])

z9 = np.array([454, 450, 445, 444, 440])
z8 = np.array([452, 452, 452, 453, 454])

z13 = np.array([444, 436, 426, 404, 374])
z14 = np.array([423, 384, 318, 220, 122])


""" PERTES DE CHARGE GÉNÉRALES """
print('1. Pertes de charge générales')
print('-----------------------------')

QL = debit**2 * L

# Conduite foncée
J_gen_fonce = (z3 - z4) * 1e-3

# Conduite claire
J_gen_clair = (z8 - z9) * 1e-3


def pdc_gen(J, couleur):
    
    if couleur == 'clair':
        print()
        print('Conduite Bleu Clair:')
        A = A_clair
    else:
        print()
        print('Conduite Bleu Fonce:')
        A = A_fonce

    vitesse = debit / (rho * A)
    Reynolds = vitesse * D_fonce / nu
    phi = 2 * g * D_fonce * J / (vitesse**2 * L)
    # Plotter ces valeurs sur le diagramme de Moody
       
    # Plot des mesures
    params()
    plt.scatter(QL, J, label = 'Mesures', color = 'black')
    plt.ylabel(r'$J_{gen}$', fontsize = 15)
    plt.xlabel(r'$Q^2 \cdot L$')
    
    # Regression lineaire
    K, intercept, r, _, _ = stats.linregress(QL, J)
    
    # Plot de la régression lineaire
    x = np.linspace(0, np.max(QL), 100)
    y = x * K + intercept
    plt.plot(x, y, label = 'Regression linéaire', color = 'black')
    plt.legend()
    
    if save:
        name = './images/pertes_generales_' + couleur + '.pdf'
        plt.savefig(name, bbox_inches='tight')
        
    plt.show()


    print('K = ', K)
    print('R^2 =', r**2)

    return K


K_fonce = pdc_gen(J_gen_fonce, 'fonce')
K_clair = pdc_gen(J_gen_clair, 'clair')


""" PERTES DE CHARGE DANS LES COUDES """
print()
print('2. Pertes de charge dans les coudes')
print('-----------------------------------')

def pdc_coude(J, coude):
    
    if coude == 'C':
        print()
        print('Coude C - Foncé')
        A = A_clair
    else:
        print()
        print('Coude J - Clair')
        A = A_fonce
    
    K_prime =  1 / (2 * g * A**2)
    QN = debit**2 * K_prime
    
    params()
    plt.scatter(QN, J, label = 'Mesures', color = 'black')
    plt.ylabel(r'$J_{loc}$', fontsize = 15)
    plt.xlabel(r'$Q^2 \cdot N$')
    
    N, intercept, r, _, _ = stats.linregress(QN, J)
    
    
    x = np.linspace(0, np.max(QN), 100)
    y = x * N + intercept
    
    plt.plot(x, y, label = 'Regression linéaire', color = 'black')
    plt.legend()
    
    if save:
        name = './images/pertes_locales_' + coude + '.pdf'
        plt.savefig(name, bbox_inches='tight')
    
    
    plt.show()
    
    print('K_prime =', K_prime)
    print('N =', N)
    print('R^2 =', r**2)


    
    
# Coude C - Fonce
J_tot_C = (z1 - z2) * 1e-3
J_gen_C = debit**2 * K_fonce * L
J_coude_C = J_tot_C - J_gen_C

pdc_coude(J_coude_C, 'C')

# Coude J - Clair
J_tot_J = (z13 - z14) * 1e-3
J_gen_J = debit**2 * K_clair * L
J_coude_J = J_tot_J - J_gen_J

pdc_coude(J_coude_J, 'J')





