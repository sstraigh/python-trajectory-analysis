def order_distances(traj, num_atms, box_width, num_frames):
    '''Takes 1d traj[] (atom locations per frame) and
    computes the 3d ordered_distances[frame][cluster_center][cluster_member] = 
    float(distance between cluster center and member)'''

    ordered_distances_per_frame_per_molec = []

    for frame in xrange(0, num_frames):

        frame_distances = []

        for cluster_center in xrange(0, int(num_atms/3)):
            for cluster_hydrate in xrange(0, int(num_atms/3)):

                dx = (float(traj[frame][cluster_center*3][1])-\
                     float(traj[frame][cluster_hydrate*3][1]))

                if (dx > box_width/2.0): dx = dx - box_width
                if (dx <= -box_width/2.0): dx = dx + box_width

                dy = (float(traj[frame][cluster_center*3][2])-\
                      float(traj[frame][cluster_hydrate*3][2]))

                if (dy > box_width/2.0): dy = dy - box_width
                if (dy <= -box_width/2.0): dy = dy + box_width

                dz = (float(traj[frame][cluster_center*3][3])-\
                      float(traj[frame][cluster_hydrate*3][3]))

                if (dz > box_width/2.0): dz = dz - box_width
                if (dz <= -box_width/2.0): dz = dz + box_width

                distance = (dx**2  + dy**2 + dz**2)**0.5

                next_tuple = [distance, cluster_hydrate, cluster_center]

                frame_distances.append(next_tuple)

        frame_distances = sorted(frame_distances, key = lambda x : (x[2], x[0]))

        ordered_distances_per_frame_per_molec.append(frame_distances)

    return ordered_distances_per_frame_per_molec

