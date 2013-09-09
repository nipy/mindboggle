#!/usr/bin/env python
"""
Functions for propagating fundus line endpoints to tile the surface.

Authors:
    - Oliver Hinds, 2013  (ohinds@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def propagate_fundus_lines(fundus_lines_file):
    from mindboggle.utils.io_vtk import read_vtk, read_scalars, write_vtk
    from mindboggle.utils.mesh import find_neighbors

    faces, _, _, points, num_points, _, fundus_lines, _ = read_vtk(
        fundus_lines_file, return_first=True, return_array=True)

    neighbor_lists = find_neighbors(faces, num_points)

    endpoints = _find_fundus_line_endpoints(fundus_lines, neighbor_lists)


    import ipdb; ipdb.set_trace()


def _find_fundus_line_endpoints(fundus_lines, neighbor_lists):
    endpoints = []

    open_fundus_line_indices = np.where(fundus_lines >= 0)
    while len(open_fundus_line_indices) > 0:
        cur_ind = pop(open_fundus_line_indices)
        num_fundus_line_neighbors = len([x for x in neighbor_lists[cur_ind]
                                         if fundus_lines[x] >= 0])
        if num_fundus_line_neighbors == 1:
            endpoints.append(cur_ind)

    return endpoints

def _find_path_to_nearest_endpoint(points, faces, fundus_line_indices,
                                   endpoint_indices):
    pass


def main(argv):
    propagate_fundus_lines(argv[1])

if __name__ == "__main__":
    sys.exit(main(sys.argv))
