from __future__ import division

import numpy as np
import collections
import matplotlib.pyplot as plt
import sys

from numpy import mean

from read_trajectory import grouper
from read_trajectory import read_trajectory

from order_distances import order_distances

from stat_mech_tools import determine_positional_probability
from stat_mech_tools import stat_mech_corr_func

if (len(sys.argv) < 6):
    print "Calling: python exe.py filename lower_bound_reactant, upper_bound_reactant,"
    print "	                  lower_bound_product, upper_bound_product"
    sys.exit()

filename = sys.argv[1]
lower_bound_reactant = float(sys.argv[2])
upper_bound_reactant = float(sys.argv[3])
lower_bound_product = float(sys.argv[4])
upper_bound_product = float(sys.argv[5])

trajectory, width, num_atoms, num_frames = read_trajectory(filename)

ordered_distance = order_distances(trajectory, num_atoms, width, num_frames)

stable_reactant_probs = determine_positional_probability(ordered_distance, num_frames, lower_bound_reactant, upper_bound_reactant, num_atoms)
stable_product_probs = determine_positional_probability(ordered_distance, num_frames, lower_bound_product, upper_bound_product, num_atoms)

corr_func = stat_mech_corr_func(stable_reactant_probs, stable_product_probs, num_atoms, num_frames)

time_domain = np.linspace(0.002, 26.982, num = num_frames)

plt.plot(time_domain, corr_func)

#eq_8 = np.ones(len(corr_func))
#normalized_corr_func = corr_func / max(corr_func)
#eq_8 = tuple(np.subtract(eq_8, normalized_corr_func))

#plt.plot(time_domain, eq_8)

plt.savefig("test.jpg")
