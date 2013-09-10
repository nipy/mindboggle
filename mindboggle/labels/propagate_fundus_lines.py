#!/usr/bin/env python
"""
Functions for propagating fundus line endpoints to tile the surface.

Authors:
    - Oliver Hinds, 2013  (ohinds@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def propagate_fundus_lines(surf_file, fundus_lines_file):
    from mindboggle.utils.io_vtk import read_vtk, read_scalars, write_vtk
    from mindboggle.utils.mesh import find_neighbors

    faces, _, _, points, num_points, fundus_lines, _, _ = read_vtk(
        surf_file, return_first=True, return_array=True)

    fundus_lines, _ = read_scalars(fundus_lines_file)

    neighbor_lists = find_neighbors(faces, num_points)

    endpoints = _find_fundus_line_endpoints(fundus_lines, neighbor_lists)

    closed_fundus_lines = _close_fundus_lines(points, fundus_lines,
                                              neighbor_lists, endpoints)

    import ipdb; ipdb.set_trace()


def _find_fundus_line_endpoints(fundus_lines, neighbor_lists):
    import numpy as np
    endpoints = []

    open_fundus_line_indices = np.where(fundus_lines >= 0)[0].tolist()

    while len(open_fundus_line_indices) > 0:
        cur_ind = open_fundus_line_indices.pop()
        num_fundus_line_neighbors = len([x for x in neighbor_lists[cur_ind]
                                         if fundus_lines[x] >= 0])
        if num_fundus_line_neighbors == 1:
            endpoints.append(cur_ind)

    return endpoints

def _close_fundus_lines(points, fundus_lines,
                        neighbor_lists, endpoint_indices):
    import numpy as np
    import Queue

    import ipdb; ipdb.set_trace()

    closed_fundus_lines = fundus_lines

    for endpoint in endpoint_indices:
        open_queue = Queue.PriorityQueue()
        distance_map = {}
        predicessor_map = {}
        closed_set = set()

        cur = [0, endpoint]
        open_queue.put(cur)
        distance_map[endpoint] = cur
        predicessor_map[endpoint] = None
        first_point = None

        while not open_queue.empty():
            dist, vertex = open_queue.get()
            if vertex in endpoint_indices:
                first_point = predicessor_map[vertex]
                break

            added_any = False
            for neighbor in neighbor_lists(vertex):
                if neighbor not in closed_set:
                    neighbor_dist = np.dot(points[vertex], points[neighbor])
                    if (neighbor not in distance_map or
                        distance_map[neighbor][0] > neighbor_dist + dist):
                        distance_map[neighbor_dist][0] = neighbor_dist + dist
                        predicessor_map[neighbor] = vertex
                        added_any = True

            # Hack to force resorting of the priority queue based on
            # new values
            if added_any:
                dummy = [-1, -1]
                open_queue.put(dummy)
                open_queue.get()

        # Trace back predicessors
        if first_point is None:
            raise RuntimeError("No path found for endpoint %d" % endpoint)

        print "found path from %d to %d" % (endpoint, first_point)

        cur_point = first_point
        while cur_point != endpoint:
            closed_fundus_lines[cur_point] = 22
            cur_point = predicessor_map[cur_point]

    return closed_fundus_lines

def main(argv):
    propagate_fundus_lines(argv[1], argv[2])

if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv))
