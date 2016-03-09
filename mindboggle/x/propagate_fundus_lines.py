#!/usr/bin/env python
"""
Functions for propagating fundus line endpoints to tile the surface.

Authors:
    - Oliver Hinds, 2013  (ohinds@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def propagate_fundus_lines(surf_file, fundus_lines_file, thickness_file):
    """Propagate fundus lines to tile the surface.

    Parameters
    ----------
    surf_file: file containing the surface geometry in vtk format
    fundus_lines_file: file containing scalars representing fundus lines
    thickness_file: file containing cortical thickness scalar data
    (for masking out the medial wall only)

    Returns
    -------
    scalars indicating whether each vertex is part of the closed
    fundus lines or not
    """
    from mindboggle.mio.vtks import read_vtk, read_scalars

    points, indices, lines, faces, fundus_lines, scalar_names, num_points, \
        input_vtk = read_vtk(surf_file, return_first=True, return_array=True)

    fundus_lines, _ = read_scalars(fundus_lines_file)
    fundus_line_indices = [i for i, x in enumerate(fundus_lines) if x > 0.5]

    thickness, _ = read_scalars(thickness_file,
                             return_first=True, return_array=True)

    return propagate_fundus_lines(
        points, faces, fundus_line_indices, thickness)

def propagate_fundus_lines(points, faces, fundus_line_indices, thickness):
    """Propagate fundus lines to tile the surface.

    Parameters
    ----------
    surf_file: file containing the surface geometry in vtk format
    fundus_lines_file: file containing scalars representing fundus lines
    thickness_file: file containing cortical thickness scalar data
    (for masking out the medial wall only)

    Returns
    -------
    scalars indicating whether each vertex is part of the closed
    fundus lines or not
    """
    from mindboggle.guts.mesh import find_neighbors
    import numpy as np

    num_points = len(points)
    neighbor_lists = find_neighbors(faces, num_points)

    # Find the boundary of the cc and call that a fundus line
    cc_inds = [x for x in range(num_points) if thickness[x] < 0.001]
    cc_boundary = [x for x in cc_inds if len([y for y in neighbor_lists[x]
                                              if y not in cc_inds])]

    fundus_line_indices += cc_boundary

    endpoints = _find_fundus_line_endpoints(
        fundus_line_indices, neighbor_lists)

    closed_fundus_lines = _close_fundus_lines(points, fundus_line_indices,
                                              neighbor_lists, endpoints)
    closed_fundus_line_indices = np.where(
        np.array(closed_fundus_lines) > 0)[0].tolist()
    new_endpoints = _find_fundus_line_endpoints(closed_fundus_line_indices,
                                                neighbor_lists)

    new_closed_fundus_lines = _close_fundus_lines(
        points, closed_fundus_line_indices, neighbor_lists, new_endpoints)

    return new_closed_fundus_lines, points, faces


def _find_fundus_line_endpoints(fundus_line_indices, neighbor_lists):
    import numpy as np
    endpoints = []

    for cur_ind in fundus_line_indices:
        num_fundus_line_neighbors = _count_fundus_line_neighbors(
            cur_ind, neighbor_lists, fundus_line_indices)

        if num_fundus_line_neighbors < 2:
            endpoints.append(cur_ind)

    return endpoints

def _count_fundus_line_neighbors(index, neighbor_lists, fundus_line_indices):
    return len([x for x in neighbor_lists[index] if x in fundus_line_indices])

def _close_fundus_lines(points, fundus_line_indices,
                        neighbor_lists, endpoint_indices):
    import numpy as np
    closed_fundus_lines = np.zeros(len(points))
    closed_fundus_lines[fundus_line_indices] = 1

    for endpoint in endpoint_indices:
        path = _find_shortest_path_to_any(
            endpoint, endpoint_indices, points, neighbor_lists)
        closed_fundus_lines[path] = 1.

    return closed_fundus_lines

def _find_shortest_path_to_any(origin, targets, points, neighbor_lists):
    import numpy as np
    import Queue

    open_queue = Queue.PriorityQueue()
    distance_map = {}
    predicessor_map = {}
    closed_set = set()

    cur = [0, origin]
    open_queue.put(cur)
    distance_map[origin] = cur
    predicessor_map[origin] = None
    first_point = None

    while not open_queue.empty():
        dist, vertex = open_queue.get()
        if vertex != origin and vertex in targets:
            first_point = predicessor_map[vertex]
            break

        added_any = False
        for neighbor in neighbor_lists[vertex]:
            if neighbor not in closed_set:
                neighbor_dist = np.sum(
                    (np.array(points[vertex]) -
                     np.array(points[neighbor]))**2)
                if neighbor not in distance_map:
                    neighbor_obj = [dist + neighbor_dist, neighbor]
                    distance_map[neighbor] = neighbor_obj
                    open_queue.put(neighbor_obj)
                    predicessor_map[neighbor] = vertex

                if neighbor_dist + dist < distance_map[neighbor][0]:
                    distance_map[neighbor][0] = neighbor_dist + dist
                    predicessor_map[neighbor] = vertex
                    added_any = True

        closed_set.add(vertex)

        # Hack to force resorting of the priority queue based on
        # new values
        if added_any:
            dummy = [-1, -1]
            open_queue.put(dummy)
            open_queue.get()

    # Trace back predicessors
    if first_point is None:
        raise RuntimeError("No path found for endpoint %d" % endpoint)

    path = []
    cur_point = first_point
    while cur_point != origin:
        path.append(cur_point)
        cur_point = predicessor_map[cur_point]
    path.append(origin)
    path.reverse()
    return path

def main(argv):
    from mindboggle.mio.vtks import write_vtk
    closed_fundus_lines, points, faces = propagate_fundus_lines(argv[1],
                                                                argv[2],
                                                                argv[3])
    write_vtk(argv[4], points, faces=faces, scalars=closed_fundus_lines,
              scalar_type='int')

if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv))
