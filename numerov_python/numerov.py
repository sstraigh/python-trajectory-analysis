from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys

from numpy import linalg as LA
from scipy import integrate
from scipy import interpolate 
from scipy import constants
from wavefunction_tools import converged
from wavefunction_tools import find_wavefunction
from wavefunction_tools import find_transition_dipole

hartree_to_wavenumbers = 219474.63
kcal_to_hartree = 0.00159362
angstroms_to_bohr = 1.889725989
kcal_to_wavenumbers = 349.75
reduced_mass = 3260.9

freq_file = open("freq_file_py.dat", 'w')
tdip_file = open("trans_dip_py.dat", 'w')

if (len(sys.argv) == 1):
	print "Usage: python numerov.py BOX_WIDTH"
	sys.exit()

BOX_WIDTH = float(sys.argv[1])
x = np.linspace(0.5, 1.8, 2000)

x_points, pot = np.loadtxt("STATIS_POT", unpack=True)
tmp0, dipx, dipy, dipz, dip_mag = np.loadtxt("dipole_surface.dat", unpack=True)
tmp1, tmp2, Ox, Oy, Oz, Dx, Dy, Dz, Hx, Hy, Hz = np.loadtxt("fort.20001", unpack=True)

#loop over all frames now they are stored as variables

points_per_frame = int(x_points[0])
num_frames = int(len(Ox) / points_per_frame)

for i in xrange(0, num_frames):

	#first, fill and spline necessary arrays
	position = []
	potential = []
	reduced_dipole = []

	unit_OD_x = []
	unit_OD_y = []
	unit_OD_z = []

	for j in xrange(0, points_per_frame):

		position.append(float(x_points[j + 1 + i*(points_per_frame + 1)]))
		potential.append(float(pot[j + 1 + i*(points_per_frame + 1)]))

		#reducing the dipole...
		dipole_mag = dip_mag[j + (i*points_per_frame)]
		dipole_z = dipz[j + (i*points_per_frame)]
		dipole_y = dipy[j + (i*points_per_frame)]
		dipole_x = dipx[j + (i*points_per_frame)]

		#dont forget about obtaining the OD unit vector!
		dx = Ox[j + (i*points_per_frame)] - Dx[j + (i*points_per_frame)]
		dy = Oy[j + (i*points_per_frame)] - Dy[j + (i*points_per_frame)]
		dz = Oz[j + (i*points_per_frame)] - Dz[j + (i*points_per_frame)]

		#taking into account periodic boundary conditions...
		if (dx >   BOX_WIDTH/2.0): dx = dx - BOX_WIDTH
		if (dx <= -BOX_WIDTH/2.0): dx = dx + BOX_WIDTH

		if (dy >   BOX_WIDTH/2.0): dy = dy - BOX_WIDTH
		if (dy <= -BOX_WIDTH/2.0): dy = dy + BOX_WIDTH

		if (dz >   BOX_WIDTH/2.0): dz = dz - BOX_WIDTH
		if (dz <= -BOX_WIDTH/2.0): dz = dz + BOX_WIDTH

		norm = (dx ** 2 + dy ** 2 + dz ** 2)**0.5

		unit_OD_x.append(float(dx/norm))
		unit_OD_y.append(float(dy/norm))
		unit_OD_z.append(float(dz/norm))

		reduced_dipole.append(float(unit_OD_x[j] * dipole_x + \
					    unit_OD_y[j] * dipole_y + \
					    unit_OD_z[j] * dipole_z))
		
	tck_pos = interpolate.splrep(position, potential, s=0)
	potential_splined = interpolate.splev(x, tck_pos, der=0)

	tck_dip = interpolate.splrep(position, reduced_dipole, s=0)
	dipole_splined = interpolate.splev(x, tck_dip, der=0)
	
	potential_splined = potential_splined * kcal_to_hartree
	potential_splined = potential_splined - min(potential_splined)
	xnew  = x * angstroms_to_bohr

	#find relevent eigenfunctions and eigenvalues

#	plt.plot(xnew, dipole_splined)
#	plt.show()


	wave0, energy0 = find_wavefunction(0, potential_splined, xnew, reduced_mass)
	wave1, energy1 = find_wavefunction(1, potential_splined, xnew, reduced_mass)
	wave2, energy2 = find_wavefunction(2, potential_splined, xnew, reduced_mass)
	
	freq_output = str(tmp1[i * points_per_frame]) + "\t" + \
		      str((energy1 - energy0) * hartree_to_wavenumbers) + "\t" + \
		      str((energy2 - energy1) * hartree_to_wavenumbers) + "\n"

	freq_file.write(freq_output)

#	prob_0 = [a**2 for a in wave0]
	
#	norm = integrate.cumtrapz(prob_0[:-100], xnew[:-100])
#	print norm[len(norm) - 1]

#	wave0 = wave0/normalized[len(normalized) - 1]
#
#	prob_1 = [a*b for a,b in zip(wave1, wave1)]
#	normalized = integrate.cumtrapz(prob_1[:-100], xnew[:-100])
#	print normalized[len(normalized) - 1]
#
#	wave1 = wave0/normalized[len(normalized) - 1]
#	
#	prob_2 = [a*b for a,b in zip(wave2, wave2)]
#	normalized = integrate.cumtrapz(prob_2[:-100], xnew[:-100])
#	print normalized[len(normalized) - 1]

#	wave2 = wave0/normalized[len(normalized) - 1]

	tdip_01 = find_transition_dipole(wave0, wave1, dipole_splined, xnew)
	tdip_12 = find_transition_dipole(wave1, wave2, dipole_splined, xnew)

	reduced_OD_x = np.mean(unit_OD_x)
	reduced_OD_y = np.mean(unit_OD_y)
	reduced_OD_z = np.mean(unit_OD_z)

	tdip_output = str(tmp1[i * points_per_frame]) + "\t" \
		      + str(tdip_01) + "\t" + str(tdip_12) + "\t" \
		      + str(reduced_OD_x) + "\t" + str(reduced_OD_y) + "\t" + str(reduced_OD_z) + "\n"

	tdip_file.write(tdip_output)

#	dipole_splined = dipole_splined * 2.541746

#	mu_0 = [a*b for a,b in zip(wave0, dipole_splined)]
#	bra_ket_integrand = [a*b for a,b in zip(mu_0, wave1)]
#
#	bra_ket = integrate.cumtrapz(bra_ket_integrand, xnew)
#
#	print bra_ket[len(bra_ket) -1]
#	
#	plt.plot(xnew, wave0, label="Psi_0")
#	plt.plot(xnew, mu_0, label="u*Psi_0") 
#	plt.show()
#
#	plt.plot(xnew, bra_ket_integrand, label="Psi_1*u*Psi_0")
#	plt.show()
#
#	plt.plot(xnew, wave0, label="Psi_0")
#	plt.plot(xnew, wave1, label="Psi_1")
#	plt.show()
#
#	plt.plot(xnew, dipole_splined, label="Reduced_dipole_spline")
#	plt.show()
#	sys.exit()

freq_file.close()
tdip_file.close()
