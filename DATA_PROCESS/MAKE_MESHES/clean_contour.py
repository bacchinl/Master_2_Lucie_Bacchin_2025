import numpy as np
import matplotlib.pyplot as plt
import copy
import sys

year=sys.argv[1]

def on_segment(p, q, r):
    """Check if point q lies on segment pr"""
    if (q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and
        q[1] <= max(p[1], r[1]) and q[1] >= min(p[1], r[1])):
        return True
    return False

def orientation(p, q, r):
    """Find the orientation of the ordered triplet (p, q, r).
    0 -> p, q and r are collinear
    1 -> Clockwise
    2 -> Counterclockwise
    """
    val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
    if val == 0:
        return 0
    elif val > 0:
        return 1
    else:
        return 2

def do_intersect(p1, q1, p2, q2):
    """Check if line segment 'p1q1' and 'p2q2' intersect."""
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if o1 != o2 and o3 != o4:
        return True

    if o1 == 0 and on_segment(p1, p2, q1):
        return True

    if o2 == 0 and on_segment(p1, q2, q1):
        return True

    if o3 == 0 and on_segment(p2, p1, q2):
        return True

    if o4 == 0 and on_segment(p2, q1, q2):
        return True

    return False

def remove_loops(points):
    cleaned_points = [points[0]]
    count = 0
    for i in range(1, len(points) - 1):
        has_intersection = False
        for j in range(len(cleaned_points) - 1):
            if do_intersect(cleaned_points[j], cleaned_points[j + 1], points[i], points[i + 1]):
                has_intersection = True
                print('\tIntersection found between points {} and {}'.format(j, i))
                count+=1
                break
        if not has_intersection:
            cleaned_points.append(points[i])
    cleaned_points.append(points[-1])
    return np.array(cleaned_points), count

def make_new_contour(contours=[f'./../Contours/{year}_State/contours_1_{year}_Lucille_new.txt', f'./../Contours/{year}_State/contours_2_{year}_Lucille.txt'], iterations = 2):
    """Uses a list of contour points to remove loops and create a new contour (e.g. Contour_1_Clean.txt)"""
    
    #load the plot
    plt.figure()
        
    # Load the data from the files
    for contour_file in contours:
        print('------------------------------------')
        print(f'Working on {contour_file}')
        print('------------------------------------')
    
        points = np.loadtxt(contour_file)

        # Separate x and y coordinates for plotting
        x = points[:, 0]
        y = points[:, 1]
            
        # Remove loops from the contour points
        count = 15 # Initialize count to a large dummy number
        i = 0
        temp_points = copy.deepcopy(points)
        while count > 0:
            temp_points, count = remove_loops(temp_points)
            print('Iteration: {}'.format(i+1))
            print('----------------------------')
            i+=1
        cleaned_points = temp_points

        # Separate x and y coordinates for plotting
        x_cleaned = cleaned_points[:, 0]
        y_cleaned = cleaned_points[:, 1]
        
        # Save contour
        save_name = f'./../{contour_file.split(".")[-2]}_Clean.txt'
        print(f'Saving file {save_name}')
        np.savetxt(save_name, cleaned_points, fmt='%f', delimiter=' ')

        # Plot the raw contour points
        plt.plot(x, y, 'o-', markersize=2, label='Raw Contour')
        plt.plot(x_cleaned, y_cleaned, 'o-', markersize=2,label='Cleaned Contour')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Contour Plot')
        plt.legend()
    plt.show()
    
make_new_contour()
    
