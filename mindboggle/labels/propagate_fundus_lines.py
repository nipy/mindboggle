#!/usr/bin/env python
"""
Functions for propagating fundus line endpoints to tile the surface.

Authors:
    - Oliver Hinds, 2013  (ohinds@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def propagate_fundus_lines(surf_file, fundus_lines_file):
    from mindboggle.utils.io_vtk import read_vtk, read_scalars
    from mindboggle.utils.mesh import find_neighbors
    import numpy as np

    faces, _, _, points, num_points, fundus_lines, _, _ = read_vtk(
        surf_file, return_first=True, return_array=True)

    fundus_lines, _ = read_scalars(fundus_lines_file)
    fundus_line_indices = np.where(
        np.array(fundus_lines) >= 0)[0].tolist()

    neighbor_lists = find_neighbors(faces, num_points)

    endpoints = _find_fundus_line_endpoints(
        fundus_line_indices, neighbor_lists)

    closed_fundus_lines = _close_fundus_lines(points, fundus_lines,
                                              neighbor_lists, endpoints)

    return closed_fundus_lines, points, faces


def _find_fundus_line_endpoints(fundus_line_indices, neighbor_lists):
    import numpy as np
    endpoints = []

    open_fundus_line_indices = _correct_covered_faces(
        fundus_line_indices, neighbor_lists)

    while len(open_fundus_line_indices) > 0:
        cur_ind = open_fundus_line_indices.pop()
        num_fundus_line_neighbors = _count_fundus_line_neighbors(
            cur_ind, neighbor_lists, fundus_line_indices)
        if num_fundus_line_neighbors == 1:
            endpoints.append(cur_ind)

    return endpoints

def _count_fundus_line_neighbors(index, neighbor_lists, fundus_line_indices):
    return len([x for x in neighbor_lists[index] if x in fundus_line_indices])

def _correct_covered_faces(fundus_line_indices, neighbor_lists):
    """Remove vertices from fundi if it is creating a triangle or more
    at the end of the line.
    """

    corrected_indices = fundus_line_indices
    for index in fundus_line_indices:
        if index not in corrected_indices:
            continue

        for neighbor in neighbor_lists[index]:
            if neighbor not in corrected_indices:
                continue

            for neighbor_neighbor in neighbor_lists[neighbor]:
                if (neighbor_neighbor == index or
                    neighbor not in corrected_indices):
                    continue

                # Look for a neighbor that is a neighbor of a
                # neighbor, meaning they all share a face.
                if index in neighbor_lists[neighbor_neighbor]:
                    # Remove the first face index that has only two
                    # fundus index neighbors
                    for ind in [index, neighbor, neighbor_neighbor]:
                        if _count_fundus_line_neighbors(
                            ind, neighbor_lists, corrected_indices) == 2:
                            print "removing %d" % ind
                            corrected_indices.remove(ind)
                            break

    return corrected_indices

def _close_fundus_lines(points, fundus_lines,
                        neighbor_lists, endpoint_indices):
    import numpy as np
    closed_fundus_lines = np.array([1 if x >= 0 else 0 for x in fundus_lines])

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
    from mindboggle.utils.io_vtk import write_vtk
    closed_fundus_lines, points, faces = propagate_fundus_lines(argv[1], argv[2])
    write_vtk(argv[3], points, faces=faces, scalars=closed_fundus_lines)

if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv))
