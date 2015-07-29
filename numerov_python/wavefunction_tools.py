from __future__ import division
import numpy as np

from scipy import integrate 
from numpy import linalg as LA

def converged(wavefunction, e_level):

        node_counter = 0
        num_vals = len(wavefunction)

        for i in xrange(2, num_vals):
                if (wavefunction[i] * wavefunction[i-1] < 0):
                        node_counter = node_counter + 1

        if (node_counter == e_level):
#		print wavefunction[num_vals - 1]
                
                if (wavefunction[num_vals - 1] ** 2 < 0.001):
			return True, 0

        if (node_counter <= e_level):
                return False, 0

        return False, 1

def find_wavefunction(energy_level, potential, x_domain, reduced_mass):

        energy_guess = 0.1
	increment = 0.1

        found = False

	step_size = x_domain[1] - x_domain[0]

        while found == False:

        	psi = np.empty(len(potential), dtype='float64')
	        psi[0]=0
	        psi[1]=0.0001

		k_2 = 2 * reduced_mass * (energy_guess - potential)

                for i in xrange(2, len(potential)):
 
                        psi[i]=  ((2 * (1 - (5.0/12.0 * step_size * step_size * k_2[i-1])) * psi[i-1])\
				 - (1 + (step_size * step_size * k_2[i-2] /12.0)) * psi[i-2])\
				 / (1 + (step_size * step_size * k_2[i]/12.0))

#      g = ((2 * potential) - (2 * energy_guess)) * reduced_mass

#               for i in xrange(2, len(potential)):

#                       psi[i]=  (2 * psi[i-1] - psi[i-2] + 5*g[i-1]*psi[i-1]\
#                                *(step_size**2)/6 + g[i-2]*psi[i-2]*(step_size**2)/12)\
#                                /(1 - g[i]*(step_size**2)/12)

		norm = 0.0
		for i in xrange(0, len(psi)):
			if (x_domain[i] > 0.9 and x_domain[i] < 2.94):
				norm = norm + (step_size * psi[i]**2)

		psi = psi/(norm**0.5)

                found, lower_upper = converged(psi, energy_level)

                if (found == True):
                        break

                if (lower_upper == 0):
                        energy_guess = energy_guess + increment

                if (lower_upper == 1):
                        energy_guess = energy_guess - increment
                        increment = increment / 10.0

#	prob = [a**2 for a in psi]	
#	area = integrate.cumtrapz(prob, x_domain)
#	print area[len(area) - 1]

	return psi, energy_guess

def find_transition_dipole(ground_state, excited_state, dipole_moment, x_domain):
	
	au_to_debye = 2.541746 
	tdip = 0.0
	step_size = x_domain[1] - x_domain[0]
	for i in xrange(0, len(ground_state)):

		if (x_domain[i] > 0.9 and x_domain[i] < 2.94):
			tdip = tdip + (ground_state[i] * excited_state[i] * dipole_moment[i] * step_size)
	
	return tdip
