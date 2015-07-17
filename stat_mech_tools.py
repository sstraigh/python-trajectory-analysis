import numpy as np
from numpy import mean

def determine_positional_probability(ordered_distances, frame_count,
                                   lower_bound, upper_bound, num_atoms):
    '''Determine probability for a molecule to be within upper and lower
    bounds of distance cut-offs w.r.t. the center of the solvation shell
    of interest'''

    num_molecs = int(num_atoms/3)
    total_stability_probabilities = []

    for cluster_center in xrange(0, num_molecs):

        cluster_stability_probabilities = []

        for cluster_rxtnt in xrange(1, num_molecs):

            molecular_time_dependent_stability = []

            for frame in xrange(0, frame_count):

                count_stable_position = 0

                if (ordered_distances[frame][(cluster_center*(num_molecs)+cluster_rxtnt)][0] < upper_bound \
                    and ordered_distances[frame][(cluster_center*(num_molecs))+cluster_rxtnt][0] > lower_bound):
                    count_stable_position = 1

                molecular_time_dependent_stability.append(count_stable_position)

            cluster_stability_probabilities.append(molecular_time_dependent_stability)

        total_stability_probabilities.append(cluster_stability_probabilities)

    return total_stability_probabilities

def stat_mech_corr_func(reactant_probs_per_cluster, product_probs_per_cluster, num_atms, num_frames):
    '''Compute and average all possible positional probability time-correlation-functions'''
   
    num_molecs=int(num_atms/3)

    cluster_correlations = []
   
    for cluster_center in xrange(0, num_molecs):

        reactant_correlations = []

        for reactant in xrange(0, num_molecs-1):

            correlation = np.correlate(reactant_probs_per_cluster[cluster_center][reactant],
                                       product_probs_per_cluster[cluster_center][reactant], mode = 'full')

            correlation = correlation[correlation.size/2:]

            reactant_correlations.append(correlation)

        average_reactant_corr = tuple(map(mean, zip(*reactant_correlations)))
        cluster_correlations.append(average_reactant_corr)

    return (tuple(map(mean, zip(*cluster_correlations))))
