def grouper(n, iterable):
    args = [iter(iterable)] * n
    return zip(*args)

def read_trajectory(trajectory_file):
    frame_coordinates = []
    filename = open(trajectory_file, 'r')
    first_line = filename.readline()
    num_atoms = int(first_line.strip().split(" ")[len(first_line.strip().split(" "))-1])
    second_line = filename.readline()
    box_width = float(second_line.strip().split(" ")[0])
    filename.seek(0,0)
    
    frame_length = num_atoms + 4
    file_length = sum(1 for line in open(trajectory_file))
    frame_count = int(file_length/frame_length)
    
    counter=1
    for line in filename:       
        if (counter > 4):
            next_line = line.strip().split(" ")
            next_tuple = []
            for i in next_line:
                if (len(str(i)) > 1):
                    next_tuple.append(i)
            frame_coordinates.append(next_tuple)        
        counter = counter + 1
        
        if (counter > (4+num_atoms)):
            counter = 1

    frame_coordinates = grouper(num_atoms, frame_coordinates)
    
    return frame_coordinates, box_width, num_atoms, frame_count
