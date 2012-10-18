
def prune(Path, Degree, TreeNbr, Terminal, Branching, Special, Vertices):
    """
    Prune an MST by deleting edges

    Note: Only links between special vertices are left.

    Parameters
    ----------
    Path : list of lists (2-tuple) of integers
        Each element of *Path* is a list of the two terminals of each pair
        of connected vertices
        Note: vertices are indexed LOCALLY
    Degree : list of integers
        Degrees of nodes in a component
    TreeNbr : list of list of integers
        Each element is a list of neighbors of a node. All LOCAL IDs.
    Terminal : list of integers
        Local IDs of nodes that are terminals.
    Branching : list of integers
        Local IDs of nodes that are branching nodes.
    Special : list of integers
        Local IDs of nodes that are special vertices connected by MST
    Vertices: list of integers
        A list of global indices to vertices

    Returns
    -------
    Path : list of lists (2-tuple) of integers
        Each element of *Path* is a list of the two terminals of each pair of
        connected vertices
        Note: vertices are indexed LOCALLY

    """

    NodeColor = {} # to visualize visited nodes of different types
    for T in Terminal:
        if len(Branching) < 1:
            break

        if T in Special:
            continue

        Trace = [ ] # store the trace visited from a terminal node to a branching node
        Visited = []  # store nodes that have been visited
        At = T  # the node of current position in tracing
        NodeColor[Vertices[T]] = 1
        while(not At in Branching and not At in Special):
            Visited.append(At)
            for Nbr in TreeNbr[At]:
                if not Nbr in Visited: # search toward the mainstream
                    Trace.append((Nbr, At))
                    Previous = At  # the node visited before At in tracing
                    At = Nbr
                    break
        # After the while loop, At stops at a branching node.
        # Delete all links from non-special terminals
        NodeColor[Vertices[At]] = 2  # just a regular branching node
        TreeNbr[At].remove(Previous)
        for Pair in Trace:
            (Src, Dst) = Pair
            if Pair in Path:
                Path.remove(Pair)
            else:  # it is possible the order of nodes is reversed in Path
                Path.remove((Dst, Src))

        if At in Branching:  # may stop at a Special-only node
            Degree[At] -= 1
            if Degree[At] < 3:
                Branching.remove(At)

    return Path, NodeColor

class Prim:  # modified from the code without license @http://hurring.com/scott/code/python/mst_prim/v0.1/mst_prim.py

    def __init__(self, A, r):
        """
        Prepare the inputs for mst_prim
        """
        import numpy as np

        self.INFINITY = np.Inf
        self.init_adjacency(A)
        self.remove_route(A, r)
        self.degree = np.zeros(len(A))
        self.tree_nbr= [[] for x in A]

    def mst_prim(self, A, w, path, degree, tree_nbr):
        """
        'A' is the adjacency matrix
        'w' is the list of all connected vertices (in order of discovery)
        'path' is a list of tuples showing (from, to)
        """
        import sys
        import numpy as np

        sys.setrecursionlimit(30000)

        # Stop when we've added all nodes to the path
        # (number of 1-D arrays within matrix that contain nonzero elements)
        if len(w) == sum([any(np.array(x)[0]) for x in A]):
            return (A, w, path, degree, tree_nbr)

        # Find minimum path coming OUT of the known vertices
        vfrom, vto, vcost = self.find_min(A, w)

        # Increase the degree for vertices vfrom and vto
        degree[vfrom] += 1
        degree[vto] += 1

        # Update tree_nbr list for vfrom and vto
        tree_nbr[vfrom].append(vto)
        tree_nbr[vto].append(vfrom)

        # Store this vertex as part of the MST path
        w.append(vto)
        path.append((vfrom, vto))

        self.remove_route(A, vto)

        return self.mst_prim(A, w, path, degree, tree_nbr)


    def init_adjacency(self, A):
        """
        Initialize adjacency list - set 0 values to INFINITY
        """
        A[A==0] = self.INFINITY

    def remove_route(self, A, v):
        """
        Once we've added a node to our path, set all routes
        to this node equal to INFINITY - to prevent loops
        """
        A[:,v] = self.INFINITY

    def find_min(self, A, w):
        """
        Find the cheapest connection we can possibly make,
        given the partially-built MST 'w' --
        'vfrom': vertex to connect from
        'vto': vertex to connect to
        'vcost': cost of connection
        """
        import numpy as np

        vcost = self.INFINITY
        vto = vfrom = -1
        for v in w:
            # Get array offset of minimum of this vertex
            i = np.argmin(A[v,:])
            if A[v,i] < vcost:
                vcost = A[v,i]
                vto = i
                vfrom = v
        return (vfrom, vto, vcost)

def mst(adjacency_matrix, indices_to_connect):
    """
    Using Prim algorithm to connect vertices within surface region

    Parameters
    ----------
    indices_region : list of integers
        indices of vertices of region in a surface mesh
    indices_to_connect : list of integers
        indices of vertices we wish to connect
    adjacency_matrix : list of lists of integers
        adjacency matrix of vertices

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> import networkx as nx
    >>> from mindboggle.utils.io_vtk import load_scalar  #, rewrite_scalars
    >>> from mindboggle.utils.mesh_operations import find_neighbors, inside_faces
    >>> import mindboggle.utils.meshlib as ml
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> sulci_file = os.path.join(data_path, 'results', 'features',
    >>>                           '_hemi_lh_subject_MMRR-21-1', 'sulci.vtk')
    >>> points, faces, sulci, n_vertices = load_scalar(sulci_file, True)
    >>> sulcus_ID = 1
    >>> sulcus_indices = [i for i,x in enumerate(sulci) if x == sulcus_ID]
    >>> sulcus_faces = inside_faces(faces, sulcus_indices)
    >>> sulcus_neighbor_lists = find_neighbors(sulcus_faces, len(points))
    >>> G=nx.Graph()
    >>> G.add_nodes_from(sulcus_indices)
    >>> for i, sulcus_neighbor_list in enumerate(sulcus_neighbor_lists):
    >>>     G.add_edges_from([[i,x] for x in sulcus_neighbor_list])
    >>> adjacency_matrix = nx.adjacency_matrix(G, nodelist=None, weight='weight')
    >>> indices_to_connect = [0, len(sulcus_indices)-1]
    >>> adjacency_matrix2, W, Path, Degree, TreeNbr = ml.mst(adjacency_matrix, indices_to_connect)
    >>> # Write results to vtk file and view with mayavi2:
    >>> MST = np.zeros(len(points))
    >>> MST[W] = 1
    >>> rewrite_scalars(sulci_file, 'test_mst.vtk', MST, MST)
    >>> os.system('mayavi2 -m Surface -d test_mst.vtk &')

    """

    if len(indices_to_connect) > 1:
        Root = indices_to_connect[0]
        M = Prim(adjacency_matrix, Root)
        adjacency_matrix, W, Path, Degree, TreeNbr = M.mst_prim(adjacency_matrix,
            [Root], [], M.degree, M.tree_nbr)

        return W, Path, Degree, TreeNbr

if __name__ == "__main__" :
    import os
    import networkx as nx
    from mindboggle.utils.io_vtk import load_scalar  #, rewrite_scalars
    from mindboggle.utils.mesh_operations import find_neighbors, inside_faces
    import mindboggle.utils.meshlib as ml
    data_path = os.environ['MINDBOGGLE_DATA']
    sulci_file = os.path.join(data_path, 'results', 'features',
                           '_hemi_lh_subject_MMRR-21-1', 'sulci.vtk')
    points, faces, sulci, n_vertices = load_scalar(sulci_file, True)
    sulcus_ID = 1
    sulcus_indices = [i for i,x in enumerate(sulci) if x == sulcus_ID]
    sulcus_faces = inside_faces(faces, sulcus_indices)
    sulcus_neighbor_lists = find_neighbors(sulcus_faces, len(points))
    G=nx.Graph()
    G.add_nodes_from(sulcus_indices)
    for i, sulcus_neighbor_list in enumerate(sulcus_neighbor_lists):
        G.add_edges_from([[i,x] for x in sulcus_neighbor_list])
    adjacency_matrix = nx.adjacency_matrix(G, nodelist=None, weight='weight')
    indices_to_connect = [0, len(sulcus_indices)-1]
    adjacency_matrix2, W, Path, Degree, TreeNbr = ml.mst(adjacency_matrix, indices_to_connect)

    # Write results to vtk file and view with mayavi2:
    MST = np.zeros(len(points))
    MST[W] = 1
    rewrite_scalars(sulci_file, 'test_mst.vtk', MST, MST)
    os.system('mayavi2 -m Surface -d test_mst.vtk &')

    Terminal, Branching = [], []
    for vtx in xrange(0,Num):
        if Degree[vtx] ==1:
            Terminal.append(vtx)
        elif Degree[vtx] > 2:
            Branching.append(vtx)

    Path, NodeColor = prune(Path, Degree, TreeNbr, Terminal, Branching, indices_to_connect, sulcus_indices)

    #endpoints = [i for i,x in enumerate(Degree) if x == 1]
