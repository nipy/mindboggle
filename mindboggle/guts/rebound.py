#!/usr/bin/env python
"""
Python module for adjusting label borders on a surface mesh.


Authors:
    - Eliezer Stavsky, 2012  (eli.stavsky@gmail.com)
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#-----------------------------------------------------------------------------
# Import Python libraries
#-----------------------------------------------------------------------------
import os
import numpy as np
from time import time
from scipy.sparse import csr_matrix, lil_matrix

from mindboggle.mio.vtks import write_vtk
import mindboggle.guts.graph as go
import mindboggle.guts.kernels as kernels

#-----------------------------------------------------------------------------
# Base label
#-----------------------------------------------------------------------------
class Bounds:
    """
    Bounds class:
    VTK I/O methods (to supplement io_vtk.py functions)
    Label propagation methods
    Label propagation:  graph-based semi-supervised learning
    Finding label boundaries

    """

    # Initialize object method
    def __init__(self):
        """
        Initialize attributes of object.
        """
        self.Points = self.Faces = self.Labels = self.Indices = 0
        self.num_points = 0

        # For label propagation
        self.seed_labels = 0
        self.min_label = 0  # ignore -1s
        self.Polylines = self.polyline_elements = 0

        # For constructing the neighbors matrix
        self.found_neighbors = 0

    ##########################################################################
    # ------------------------------------------------------------------------
    #     Label propagation methods
    # ------------------------------------------------------------------------
    ##########################################################################

    def initialize_seed_labels(self, init='lines', fraction=.05,
                               output_filename='', verbose=False):
        """
        Initialize a set of seed labels for relabeling or label propagation.

        Options include:
          - lines: vertices in polylines
          - flanks: vertices flanking polylines
          - lines_flanks: both the polylines and the flanks
          - label_boundary: label boundary vertices
          - random: a <fraction> of random vertices

        """
        self.seed_labels = np.zeros(self.num_points)

        # To initialize with vertices flanking polylines,
        # find all vertices that are part of a triangle that includes
        # at least one polyline vertex, and store a 1 in the array self.seed_labels
        if init in ['flanks','lines_flanks']:
            if verbose:
                print('Initializing seed labels with polyline-flanking '
                      'vertices...')
            self.seed_labels[self.find_polylines_flanks()] = 1

        # To initialize with polylines, find all vertices that are part of a
        # polyline, and store a 1 in the array self.seed_labels
        if init in ['lines','lines_flanks']:
            if verbose:
                print('Initializing seed labels with polyline vertices...')
            self.seed_labels[self.polyline_elements] = 1

        # To initialize with label boundaries, call find_label_boundary()
        if init == 'label_boundary':
            if verbose:
                print('Initializing seed labels with vertices of the label '
                      'boundaries')
            self.find_label_boundary(output_filename=output_filename)
            self.seed_labels[self.label_boundary] = 1

        # To initialize with a fraction of random vertices,
        # init every 1/fraction label
        if init == 'random':
            if verbose:
                print('Initializing seed labels with random vertices...')
            if fraction > 1:
                if verbose:
                    print('Please enter a fractional number less than or '
                          'equal to 1.')
                return
            randoms = np.array([np.mod(i, int(1.0/fraction))
                                for i in range(self.num_points)])
            self.seed_labels[randoms==0] = 1

        # Replace the 1s in self.seed_labels with the seed label values
        self.seed_labels[self.seed_labels==1] = self.Labels[self.seed_labels==1]

        # Provide some statistics for what was done
        self.num_seed_labels = len(self.seed_labels[self.seed_labels>0])
        self.percent_seed_labels = (self.num_seed_labels+0.0) / self.num_points * 100

        if verbose:
            print('Percentage of seed labels: {0}'.
                  format(self.percent_seed_labels))

        return self.seed_labels

    def build_label_matrix(self):
        """
        Construct a vertices x labels array of vertex label assignment values.

        Parameters
        ----------
        array of n labels; -1 corresponds to no label

        Returns
        -------
        n x C : array
            row corresponds to vertex, column corresponds to label
            1 indicates that it is assigned the column's label
           -1 indicates that it does not have that label
            0 indicates that there is no label assigned

        """

        # Get the unique (sorted) set of labels
        self.unique_labels = np.sort(np.asarray(list(set(self.seed_labels))))

        # Number of labels and vertices
        self.unique_labels = [x for x in self.unique_labels if x >= self.min_label]
        C = len(self.unique_labels)
        n = len(self.Labels)

        # Construct n x C matrix
        self.label_matrix = np.zeros((n, C))

        # Populate the label assignment matrix with -1s and 1s for seed labels
        # For each row...
        for i in range(n):
            if self.seed_labels[i] >= self.min_label:
                self.label_matrix[i, :] = -1
                unique_label_index = np.where(self.unique_labels == self.seed_labels[i])[0]
                self.label_matrix[i, unique_label_index] = 1

        self.num_labels = C

        return self.label_matrix

    def graph_based_learning(self, method='propagate_labels', realign=False,
                             kernel=kernels.rbf_kernel,
                             sigma=10, max_iters=200, tol=.001, vis=False,
                             verbose=False):
        """
        Main function to perform graph-based learning, such as label propagation.

        Other methods include "jacobi iteration" and "label spreading"

        Parameters
        ----------
        - faces and points of a vtk surface mesh
        method: string (choice of algorithm)
        realign: boolean (use label propagation for realigning boundaries?)
        kernel: function (used in constructing affinity matrix)
        sigma: float (gaussian kernel parameter)
        vis: boolean (show progress of algorithm?)
        max_iters: int (number of times to repeat the algorithm)
        tol: float (threshold to assess convergence of the algorithm)

        Returns
        -------
        self.learned_matrix: np array (n x C array of probabilities
                                       that vertex has a given label)
        """

        # Step 1. Construct affinity matrix - compute edge weights
        if self.Points.shape and self.Indices.shape and self.Faces.shape:
            self.affinity_matrix = go.weight_graph(self.Points, self.Indices,
                self.Faces, kernel=kernel, sigma=sigma, add_to_graph=False,
                verbose=False)
        else:
            if verbose:
                raise IOError("  Missing data!")

        # Step 2. Transform column of labels into n x C matrix, one column per label
        if not realign:
            self.build_label_matrix()
        else:
            self.build_label_segment_matrix()

        # Step 3. Propagate Labels!
        if method == "propagate_labels":
            if verbose:
                print('Perform weighted average algorithm (max_iters={0})'.
                    format(max_iters))
            # Construct self.learned_matrix matrix within method
            self.propagate_labels(realign, max_iters, tol, vis=vis)
        else:
            if verbose:
                print('That algorithm is not available.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #-------------------------------------------------------------------------
    #     Label propagation:  post-processing
    #-------------------------------------------------------------------------
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def assign_max_prob_label(self, verbose=False):
        """
        Assign hard labels to vertices by finding the highest probability label
        of self.learned_matrix. Output results to a VTK file.

        This method takes self.learned_matrix and determines the most
        probable labeling for each vertex.  It outputs a separate array,
        max_prob_label, which contains one label for each vertex.

        Returns
        -------
        self.max_prob_labels: np array (size num_points,
                             containing the most probable label for each vertex)

        """

        # Go row by row (one vertex at a time),
        # and find the column with the maximum value
        try:
            max_col = np.argmax(self.learned_matrix, axis=1)
        except:
            if verbose:
                print('First call graph_based_learning().')
            return

        # Define an array called max_prob_label to store the final values
        self.max_prob_labels = np.zeros(self.num_points)

        # Use the array of unique, sorted labels to convert this matrix
        # back to the original labeling; max_col[i] is the temporary label number
        for i in range(self.num_points):
            self.max_prob_labels[i] = self.unique_labels[max_col[i]]

        return self.max_prob_labels

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #-------------------------------------------------------------------------
    #     Label propagation:  graph-based semi-supervised learning
    #-------------------------------------------------------------------------
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def propagate_labels(self, realign, max_iters, tol, vis=True,
                         verbose=False):
        """
        Run iterative weighted average algorithm to propagate labels to unlabeled vertices.

        Parameters
        ----------
        realign:    boolean (propagation is for realignment?)
        max_iters:  int (number of iterations)
        tol:        float (threshold for terminating algorithm)
        vis:        boolean (incremental VTK files to visualize
                             progress of the algorithm?)
        Returns
        -------
        self.learned_matrix: np array
             (n x C matrix of probabilities that vertex belongs to a given label)

        Features: Hard label clamps, probabilistic solution.
        See: Zhu and Ghahramani, 2002.

        """

        """ The first approach to be considered in the semi-supervised learning case
        is to propagate labels on the graph.
        A simple algorithm of this sort has been propoosed by Zhu and Ghahramani (2002),
        and starts (like all such algorithms) with a set of n vertices,
        l of which are labeled, and u unlabeled.
        The algorithm takes as its input the affinity matrix W (self.affinity_matrix).
        From the affinity matrix, one may construct the diagonal degree matrix,
        which is a measure of the total weight (or number of edges) attached to a vertex."""

        self.DDM = go.diagonal_degree_matrix(self.affinity_matrix, inverse=True)

        """ Next, we must initialize a vector to represent the results of the label
        propagation algorithm. It will contain l labels and u 0's.
        This has already been done by the function initialize_seed_labels,
        and is called self.seed_labels.
        We will just check to make sure this has been accomplished."""

        if isinstance(self.seed_labels,int):
            if verbose:
                print('Please initialize the labels by calling '
                      'self.initialize_seed_labels()')
            return

        """ Now, we can actually proceed to perform the iterative algorithm.
        At each timestep, the labels will be updated to reflect the weighted
        average of adjacent vertices. An important caveat of this algorithm
        is that the labeled vertices remain fixed, or clamped.
        They should not be changed, and will need to be reset.
        We accomplish the reset by recalling that self.seed_labels
        stores the indexes of seed labels,
        and self.Labels contains the actual labels.
        The algorithm repeats itself until either convergence or max_iters
        (which will prevent excessive computation time).
        We must also take care to solve the multi-label problem.
        To do so, we employ a one-vs-all framework, where each label is
        considered independently, and set against the rest of the labels.
        More specifically, self.label_matrix is an n x C matrix, where each
        row represents a vertex and each column represents label membership.
        We can go column by column and process the algorithm
        iteratively. So, we'd start at the first column and see which vertices
        get labeled. Then we'd move to the next column and label more vertices.
        Because it is possible (likely) that some vertices will not receive any label,
        and also to account for probabilistic labeling, we will assign a probability
        of a vertex receiving a label. Then we can report these probabilities.
        So, to begin, let us first construct this probabilistic label assignment:
        This matrix will store a 1 for 100% probability, 0 for 0%, and
        fractional values for the rest.
        We will rename self.label_matrix for this purpose."""

        if not realign:
            self.learned_matrix = self.label_matrix
        else:
            # For realignment, we need to use a different labeling matrix,
            # one which has each segment assigned a different labels
            self.learned_matrix = self.label_segment_matrix

        """ We will later change the -1s to 0s.
        As vertices get labeled, we assign a confidence measure to the labeling
        and store the value in this matrix.
        Now, let us go column by column, and run the weighted averaging algorithm.
        For each column, you're studying one label. Therefore, when updating
        self.learned_matrix, you'll be working with one column at a time too.
        If a label gets vertex, keep the fractional value, do not simply round
        to 1 to assign membership."""

        i = 0 # record of label number
        for column in self.learned_matrix.T:

            t0 = time()
            if verbose:
                print('Number of initial members for label {0}: {1}'.format(
                    i, np.nonzero(column==1)[0].size))

            # Set up indices and values to be clamped during propagation
            if not realign:
                restore_indices = self.seed_labels >= self.min_label
                restore_values = column[restore_indices]
            else:
                restore_indices = np.hstack((self.label_boundary,self.polyline_elements))
                restore_values = column[restore_indices]

            Y_hat_now = csr_matrix(column).transpose()
            converged = False
            counter = 0
            while not converged and counter < max_iters:
                """ The option will exist to visualize the proceedings of the algorithm.
                The results of a number of the iterations will be sent to vtk
                files which can then be visualized.
                For the visualization, we will construct two types of vtk files.
                The first will be the actual (manual) labels, as found in
                self.Labels, with the label of interest highlighted (=1),
                and the others blanked out (=-1).
                The other vtk files will be the result of the algorithm."""
                if vis and not realign:
                    """First, we'll find out which label we're working with,
                    by calling label_mapping. We'll then send that label to the
                    method highlight() which will do the actual work of
                    creating the vtk."""
                    label = self.unique_labels[i]

                    if not counter: # No need to do this more than once :-)
                        self.highlight(label)

                    """ Next, we'll construct vtk files, assuming that the
                    iteration step is one we care about.
                    For our purposes, let's see the early iterations in high
                    density, once low density in the middle,
                    and the last iteration before convergence or max_iters.
                    So, the numbers of interest will be when counter is between
                    0 and 10, and at 100.
                    We'll also see max_iters/2, and max_iters (or convergence).
                    """
                    # Actually, just see if counter is a multiple of a number of interest.
                    # set_of_interest = np.array([0,101,202,max_iters-2,max_iters-1])

                    """ Let's define a file for output.
                    We have the vertices and faces, and we have the labels
                    which are found in column.todense().flatten()."""

                    filename = str(label)+'_'+str(counter)+'.vtk'

                    if not np.mod(counter,1000):
                        LABELS = np.zeros(self.num_points)
                        LABELS[:] = Y_hat_now.todense().T.flatten()
                        write_vtk(filename, self.Points, self.Vertices,
                                  [], self.Faces, [LABELS], scalar_type='int')

                # column matrix
                Y_hat_next = (self.DDM * self.affinity_matrix * Y_hat_now).todense()
                # reset
                Y_hat_next[restore_indices, 0] = restore_values
                # check convergence
                converged = (np.sum(np.abs(Y_hat_now.todense() - Y_hat_next)) < tol)
                # if verbose:
                # print('Iteration number {0}, convergence = {1}'.
                # format(counter,np.sum(np.abs(column.todense() - tmp)))
                Y_hat_now = csr_matrix(Y_hat_next)
                counter += 1

            # Print out the number of iterations, so that we get a sense for future runs.
            # It is also an indication of whether the algorithm converged.

            if verbose:
                if counter == max_iters:
                    print('Done in {0:.2f} seconds (the algorithm did not '
                          'converge)'.format(time()-t0))
                else:
                    print('Done in {0:.2f} seconds ({1} iterations)'.
                        format(time()-t0, counter))

            self.learned_matrix[:,i] = Y_hat_now.todense().flatten()

            #if verbose:
            #print('There were {0} initial seed vertices for this label'.
            #format(self.count_assigned_members(i))
            #print('The file actually had {0} vertices for this label'.
            #format(self.count_real_members(self.label_mapping[i]))
            #print('Using only those vertices which crossed the threshold,
            # there are now: '.format(self.count_real_members(i))
            #pylab.plot(self.learned_matrix[:,i])
            #pylab.show()

            i += 1

        """ Before reporting the probabilistic assignment, we change all -1's,
        which indicates 0 probability that the vertex has that label.
        Note: because the labels were initially numbered -1 and 1, there will
        be 'probabilities' below 0.  So, to obtain a sort of probability
        distribution which preserves order, we will add 1 to each number,
        and then divide by 2. Thus -1 --> 0, 1 --> 1 and everything else
        keeps its order."""

        self.learned_matrix += 1
        self.learned_matrix /= 2

        """ self.learned_matrix is now complete."""
        return self.learned_matrix

    ##########################################################################
    # ------------------------------------------------------------------------
    #     Finding label boundaries
    # ------------------------------------------------------------------------
    ##########################################################################

    def find_label_boundary(self, output_filename, realigned_labels=False,
                            verbose=False):
        """
        Find the vertices for all boundaries between different labels in a mesh.

        Label boundaries are the set of all vertices
        whose neighbors do not share the same label.

        NOTE: This method is redundant and inelegant.

        Parameters
        ----------
        realigned_labels: boolean (use realigned labels, not manual labels?)
        output_filename: string (vtk file containing highlighted label boundaries)

        Returns
        -------
        self.label_boundary: numpy array (of indices of vertices which comprise the label boundary)
        self.label_boundary_file: string (VTK file name containing highlighted label boundary)
        --OR--
        self.Rlabel_boundary: numpy array (of indices of vertices which comprise the realigned label boundary)
        self.Rlabel_boundary_file: string (VTK file name containing highlighted realigned label boundary)

        """
        if not realigned_labels:
            self.label_boundary = np.zeros(self.num_points)
        else:
            self.Rlabel_boundary = np.zeros(self.num_points)

        # Loop through triangles
        for triangle in self.Faces:
            v0, v1, v2 = triangle[0], triangle[1], triangle[2]
            # If they are not all the same...
            if not realigned_labels:
                same_labels = [self.Labels[v1]==self.Labels[v2],
                               self.Labels[v0]==self.Labels[v2],
                               self.Labels[v0]==self.Labels[v1]]
            else:
                same_labels = [self.RLabels[v1]==self.RLabels[v2],
                               self.RLabels[v0]==self.RLabels[v2],
                               self.RLabels[v0]==self.RLabels[v1]]

            # Did it this way for option to modify it later. Perhaps reconsider 'not all' statement.
            if not all(same_labels):
                # Then label those vertices as part of the boundary.
                if not realigned_labels:
                    self.label_boundary[triangle] = 1
                else:
                    self.Rlabel_boundary[triangle] = 1

        # We can now output a file to show the boundary.
        if not realigned_labels:
            write_vtk(output_filename, self.Points, self.Vertices,
                               self.Faces, [self.label_boundary], scalar_type='int')
            self.label_boundary_file = output_filename
        else:
            write_vtk(output_filename, self.Points, self.Vertices,
                      [], self.Faces, [self.Rlabel_boundary], scalar_type='int')
            self.Rlabel_boundary_file = output_filename

        # Reformat to include only indices of those vertices on the boundaries.
        if not realigned_labels:
            self.label_boundary = np.nonzero(self.label_boundary==1)[0]
        else:
            self.Rlabel_boundary = np.nonzero(self.Rlabel_boundary==1)[0]

        if verbose:
            print('The label boundary array is: {0}'.
                  format(self.label_boundary))

        if not realigned_labels:
            return self.label_boundary, self.label_boundary_file
        else:
            return self.Rlabel_boundary, self.Rlabel_boundary_file

    def find_label_boundary_per_label(self):
        """
        Find label boundary for each label.

        This is a helper method for find_label_boundary_segments.

        Parameters
        ----------

        Returns
        -------
        self.label_boundary_per_label: dict (key: int label, value: list of vertices)

        """

        try:
            self.label_boundary
        except AttributeError:
            self.find_label_boundary() # get the boundaries between all labels

        self.label_boundary_per_label = {}
        setA = set(self.label_boundary)

        for Class in self.set_manual_labels:
            setB = set(np.nonzero(self.Labels==Class)[0])
            setC = setA.intersection(setB)

            self.label_boundary_per_label[Class] = list(setC)

        return self.label_boundary_per_label

    def find_label_boundary_segments(self, verbose=False):
        """
        Break up the label boundaries into segments (corresponding to same-label pairs).

        This method will output a dictionary to store label boundary segments (and subsegments).
        The key is a tuple of the currently assigned label and the adjacent label.
        The value is the set of vertices which comprise the segment.

        Note:  Two different keys will correspond to two cosegments:
        The label boundary consists of a set of vertices on the surface mesh
        two thick, following from the (>=2 label neighborhood) definition,
        so we call these "label cosegments".

        Returns
        -------
        self.label_boundary_segments: dict (key: tuple of labels, value: set of vertices)
        self.segment_file: string (pickled file containing the dictionary, for future, and faster, use)
        self.highlighted_segment_file: string (VTK file with boundary segments highlighted according to label)

        """
        self.find_label_boundary_per_label()

        # Initialize dictionary for later ease of use;
        # there's probably a simpler way to do the concatenation
        self.label_boundary_segments = {}
        for a in self.set_manual_labels:
            for b in self.set_manual_labels:
                self.label_boundary_segments[(a,b)]=[]

        # Populate the dictionary with vertices
        for Class in self.set_manual_labels:
            if verbose:
                print(Class)
            for vertex in self.label_boundary_per_label[Class]:
                neighbors = self.neighbors(vertex)
                A = set(self.Labels[neighbors])
                B = set(list([self.Labels[vertex]]))
                neighbor_labels = set.difference(A,B)
                for c in neighbor_labels:
                    self.label_boundary_segments[(Class,c)] += [vertex]

        # Trim results - delete any dict entry which has no vertices
        for key, value in list(self.label_boundary_segments.items()):
            if value == []:
                self.label_boundary_segments.pop(key)

        # Print results
        if verbose:
            for key in list(self.label_boundary_segments): #.keys():
                print('For labels: {0} {1}'.
                      format(key, self.label_boundary_segments[key]))

        # Output the results to a VTK file
        self.highlighted_segment_file = 'highlighted_segments.vtk'
        color = 1
        colored_segments = np.zeros(self.Labels.shape)
        for value in list(self.label_boundary_segments.values()):
            colored_segments[value] = color
            color += 1
        write_vtk(self.highlighted_segment_file, self.Points,
                  self.Vertices, [], self.Faces, [colored_segments],
                  scalar_type='int')

        return self.label_boundary_segments, self.highlighted_segment_file

    def find_intersections(self, segment, endpoint, verbose=False):
        """
        Find the terminal vertices which are at the intersection
        of the polylines and label boundary.

        NOTE: FIX -- no endpoints

        Parameters
        ----------
        segment: list (vertices comprising the label boundary segment)
        endpoint: list (two endpoints of the label boundary segment)

        Returns
        -------
        intersection: list (two outermost vertices of the label boundary segment which intersect the polylines)
                        If there are not two intersection, return -1s in place

        """

        if verbose:
            print('Finding intersection of segment with polylines...')
        intersection = [0,0]

        for i in range(2):
            pointer = endpoint[i]
            used_vertices = [pointer]
            neighbors = []
            while pointer not in self.polyline_elements:
                tmp0 = np.intersect1d(self.neighbors(pointer),segment)
                tmp1 = np.setdiff1d(tmp0,used_vertices)
                neighbors = neighbors + list(tmp1)

                if not neighbors:
                    pointer = -1
                    if verbose:
                        print('No neighbors left to explore.')
                    break

                pointer = neighbors.pop()
                used_vertices.append(pointer)

            intersection[i] = pointer

        if verbose:
            print('The intersection points are: {0}'.format(intersection))

        if np.product(intersection) < 0:
            if verbose:
                print(segment)
            labels = np.zeros(self.Labels.shape)
            labels[segment] = 100
            labels[endpoint] = 200
            write_vtk('debug_intersections.vtk',self.Points,
                      self.Vertices, [], self.Faces, [labels],
                      scalar_type='int')

        return intersection

    def build_label_segment_matrix(self, verbose=False):
        """
        Construct a vertices x label segments array of vertex label segment assignment values.

        This method is analogous to build_label_matrix but contains values
        associating vertices with label pairs rather than with labels.

        Take the dictionary self.label_boundary_segments and convert
        it into the label matrix self.label_segment_matrix for label propagation.
        Separate segments should be assigned separate classes.
        The graph-based learning algorithms all work on some n x C label matrix.
        They start with a "sparsely populated" matrix, and fill it with probabilities.

        NOTE:  Useful if using label propagation for realignment.

        Parameters
        ----------
        Array of n label segments.
        self.label_boundary_segments: dict
            key:  2-tuple containing the assigned label and the adjacent label
            value:  list of vertices with the assigned label

        Returns
        -------
        n x C array. Row corresponds to vertex, column corresponds to label segment.
            1:   assigned the label segment corresponding to the column
           -1:   does not have that label
            0:   no label assigned

        self.label_segment_matrix: n x C array (initial labels used in propagation)
        self.realignment_mapping: dict (key-value pairs of new-old labels)

        """

        # Step 0. Find number of segments
        self.num_segments = len(self.label_boundary_segments)

        # Step 1. Construct n x self.num_segments matrix of labels,
        # and produce label mapping dictionary
        self.label_segment_matrix = np.zeros((self.num_points,self.num_segments))

        self.realignment_mapping = {}
        label = 0
        for key, value in list(self.label_boundary_segments.items()):
            self.realignment_mapping[label] = key
            self.label_segment_matrix[value,:] = -1
            self.label_segment_matrix[value,label] = 1
            label += 1

        if verbose:
            print('Mapping is: {0}'.format(self.realignment_mapping))

        self.determine_appropriate_segments()

        return self.label_segment_matrix, self.realignment_mapping

    def same_boundary(self, vertex1, vertex2):
        """
        Determines whether two label boundary vertices belong to the same label boundary.

        Parameters
        ----------
        vertex1, vertex2: int (index of label boundary vertex)

        Returns
        -------
        same: boolean (belong to the same boundary?)

        """
        same = False

        for value in list(self.label_boundary_segments.values()):
            if vertex1 in value and vertex2 in value:
                same = True

        return same

    def determine_appropriate_segments(self, proportion = 1, dist_threshold = 8,
                                       lb_fundus_threshold = 16,
                                       num_good_vertices = 5, eps=1E-7,
                                       spread_tol = 6, verbose=False):
        """
        Determine which label boundary segments should propagate their labels.

        NOTE:  Face validity but not directed to do what we need at present.

        Parameters
        ----------
        proportion: float (threshold of acceptable ratio of differences between polylines-to-boundary distances)
        dist_threshold: float (threshold of absolute distance between polylines and boundary, above which propagation is prohibited)
        num_good_vertices: int (threshold above which a label boundary segment will be preserved)
        eps: float (numerical stability - avoid division by zero)
        pickled_filename: str (pkl file storing the distance matrix. Saves time.)

        Returns
        -------
        self.label_segment_matrix: np array (n x num_segments matrix of labels, with zeros in unusable columns)

        Explanation ::

        In this method, we will prune the realignment matrix which was just constructed from the label boundary segments.
        The pruning will follow the principles that we want the polylines to be not too far from a given label boundary segment,
        and that it should be significantly closer to one than another.

        We should also preserve the shape of the label boundary.
        And the polylines should run somewhat parallel to the label boundary.
        """

        # Step 0. Construct num_polylines_vertices x num_label_boundary_vertices np array of distances:
        if verbose:
            print('Beginning determine_appropriate_segments()...')
        t0 = time()

        if verbose:
            print(self.polyline_elements.shape)
            print(self.label_boundary.shape)

        distance_matrix = np.asarray([np.linalg.norm(self.Points[x1] -
                                                     self.Points[x2])
            for x1 in self.polyline_elements
            for x2 in self.label_boundary]).\
            reshape((self.polyline_elements.size, -1))

        if verbose:
            print('Distance Matrix has been constructed in {0}. '
                  'Bounds is {1}.'.format(time() - t0, distance_matrix.shape))

        # Step 1. For each fundus vertex, find the closest and
        # second closest label boundary vertices,

        sorted_distances = np.argsort(distance_matrix)
        if verbose:
            print('Got sorted distances. Bounds is {0}'.
                  format(sorted_distances.shape))

        closest_label_boundary = sorted_distances[:,0]
        if verbose:
            print('Got closest label boundary. Bounds is {0}. '
                  'First few values are {1}'.format(
                  closest_label_boundary.shape, closest_label_boundary[:10]))

        dir = os.getcwd()
        self.highlight_vtk_vertices(self.label_boundary[closest_label_boundary],
                                    dir + '/close_vertices.vtk')

        closest_polylines = np.argsort(distance_matrix, 0)[0,:]

        closest_distances = np.amin(distance_matrix, 1)
        if verbose:
            print('Got closest_distances. Bounds is {0}. '
                  'First few values are {1}'.format(closest_distances.shape,
                                                    closest_distances[:10]))

        second_closest_distances = np.asarray([distance_matrix[i,sorted_distances[i,1]]
                                               for i in range(self.polyline_elements.size)])
        if verbose:
            print('Got second closest distances. Bounds is {0}. '
                  'First few values are {1}'.format(
                    second_closest_distances.shape,
                    second_closest_distances[:10]))

        # Let's try using a dictionary to express the mapping relationship.
        # We will have one which maps polylines vertices to nearest label boundary vertices.
        # And we'll have one which maps lb vertices to nearest polylines vertices.

        polylines_lb = dict((self.polyline_elements[i],
                             (self.label_boundary[closest_label_boundary[i]],
                              distance_matrix[i,closest_label_boundary[i]]))
                                              for i in range(self.polyline_elements.size))
        lb_polylines = dict((self.label_boundary[i],
                             (self.polyline_elements[closest_polylines[i]],
                              distance_matrix[closest_polylines[i],i])) for i in range(self.label_boundary.size))

        if verbose:
            print('The polylines to label boundary mapping is: {0}'.format(polylines_lb))
            print('The label boundary to polylines mapping is: {0}'.format(lb_polylines))

        # Step 2. Determine which obey proper proportions and distances, using parameters
        within_distance = (closest_distances < dist_threshold)
        if verbose:
            print('Got within distance. Num satisfy is {0}. First few are {1}'.
                format(within_distance.nonzero()[0].size,
                       within_distance[:10]))

        self.highlight_vtk_vertices(self.label_boundary[closest_label_boundary[within_distance==1]],
                                dir + '/close_distance.vtk')

        within_proportion = np.bitwise_or((closest_distances / second_closest_distances > proportion),
                                          (second_closest_distances / (closest_distances+eps) > proportion))
        if verbose:
            print('Got within proportion. Num satisfy is {0}. First few are {1}'.
                format(within_proportion.nonzero()[0].size,
                       within_proportion[:10]))

        self.highlight_vtk_vertices([self.label_boundary[closest_label_boundary[within_proportion==1]]],
                                dir + '/good_proportion.vtk')

        # The following array stores the indices of the label boundary vertices which satisfy the above properties.
        satisfy_distances = self.label_boundary[closest_label_boundary[np.nonzero(np.bitwise_and(within_distance,
                                                                                                 within_proportion))]]
        if verbose:
            print('Got satisfy distances. Bounds is {0}. They are {1}'.
                  format(satisfy_distances.shape, satisfy_distances))

        self.highlight_vtk_vertices(satisfy_distances, dir + '/satisfy_distance.vtk')

        if verbose:
            print('Currently, {0} vertices satisfy the distance requirement'.
                  format(satisfy_distances.size))

        # Ok, now here comes the critical step.
        # We have the array satisfy_distances. It stores the indices of the elite vertices, those which satisfy the first two properties.
        # Now, we want to prune the array, and then augment the pruned version.
        # For pruning, we will eliminate any lb vertex whose closest fundus vertex has a large spread among the lb vertices.
        # For augmenting, we will add any vertex which maps to a fundus which maps to a qualified lb vertex on the same label boundary.

        # Pruning...
        for lbvertex in satisfy_distances:
            fundus_vertex = lb_polylines[lbvertex][0]
            fundus_index = np.nonzero(self.polyline_elements == fundus_vertex)[0]
            top_five_indices = sorted_distances[fundus_index,:5].flatten()
            top_five_lbvertices = self.label_boundary[top_five_indices]
            spread_matrix = np.zeros((top_five_lbvertices.size,top_five_lbvertices.size))
            for i in range(top_five_lbvertices.size):
                for j in range(top_five_lbvertices.size):
                    v1 = top_five_lbvertices[i]
                    v2 = top_five_lbvertices[j]
                    spread_matrix[i,j] = np.linalg.norm(self.Points[v1] - self.Points[v2])

            spread = np.max(spread_matrix)
            if spread > spread_tol:
                satisfy_distances = np.delete(satisfy_distances,np.nonzero(satisfy_distances == lbvertex))
                if verbose:
                    print('deleted vertex: {0}'.format(lbvertex))

                #### AH! I'm changing that over which I'm iterating! Fix.

        self.highlight_vtk_vertices(satisfy_distances, dir + '/satisfy_distance_pruned.vtk')

        if verbose:
            print('After pruning, {0} vertices satisfy the distance requirement'.
                  format(satisfy_distances.size))

        # Augmenting...
        for lbvertex in self.label_boundary:
            fundus_vertex, distance = lb_polylines[lbvertex]
            if distance < lb_fundus_threshold:
                mapped_lbvertex = polylines_lb[fundus_vertex][0]
                if mapped_lbvertex in satisfy_distances and self.same_boundary(mapped_lbvertex,lbvertex):
                    satisfy_distances = np.append(satisfy_distances,lbvertex)
                    if verbose:
                        print('added vertex: {0}'.format(lbvertex))

        self.highlight_vtk_vertices(satisfy_distances, dir + '/satisfy_distance_pruned_augmented.vtk')
        if verbose:
            print('After augmenting, {0} vertices satisfy the distance '
                  'requirement'.format(satisfy_distances.size))

        # Now we will see how many vertices from each label boundary segment satisfy the properties.
        # If a segment only contains a few vertices, then we won't bother propagating labels from it.

        reverse_mapping = dict((v,k) for k, v in list(self.realignment_mapping.items()))

        # Let's include some information as to which label boundaries will propagate their labels...
        vertices_to_highlight = np.zeros(self.Labels.shape)

        for key, value in list(self.label_boundary_segments.items()):
            # num_intersections = np.intersect1d(satisfy_distances, value).size + np.intersect1d(satisfy_distances, self.label_boundary_segments[key[::-1]]).size
            num_intersections = np.intersect1d(satisfy_distances, value).size
            if verbose:
                print('Number of intersections is: {0}'.
                      format(num_intersections))
            if (num_intersections < num_good_vertices):
                self.label_segment_matrix[:,reverse_mapping[key]] = 0
            else:
                vertices_to_highlight[value] = 1
                if verbose:
                    print('______________Preserving Label '
                          'Boundary Segment_____________')

        write_vtk(dir + '/propagating_regions.vtk',self.Points,
                  self.Vertices, [], self.Faces, [vertices_to_highlight],
                  scalar_type='int')

        return self.label_segment_matrix

    def analyze_label_polylines_intersections(self, verbose=False):
        """
        Find polyline vertices which intersect label boundary, and determine whether they belong to the same fundus curve.

        NOTE:  Requires sequential vertices for fundus or label boundary
               WORKING?

        Runs functions: self.find_label_boundary_segments() - if necessary

        Returns
        -------
        self.label_boundary_segments: dict (updated dict, with code at end signifying what type of label propagation to use)
            0: two polylines intersections from same fundus curve - do fill operation, in essence.
            1: all other cases - do label propagation, see other methods for details.
            This list may change as the algorithm for propagation gets more sophisticated.

        """

        try:
            self.label_boundary_segments
        except AttributeError:
            self.find_label_boundary_segments(completed=completed)

        for key, segment in list(self.label_boundary_segments.items()):
            if segment:
                endpoint = self.find_endpoints(segment)
                intersection = self.find_intersections(segment, endpoint)

                self.same_fundus = False
                if -1 not in intersection and intersection[0] != intersection[1]:
                    # It found 2 different intersections. We must therefore now assess whether these two intersections are
                    # part of the same fundus.

                    pointer = intersection[0] # This is the first intersection point.
                    if verbose:
                        print('First pointer is: {0}'.format(pointer))
                    row_avoid = [] # This will be an array of rows in self.Polylines to avoid
                    vertex_avoid = [pointer] # This will be an array of vertices to avoid

                    rows = list(np.nonzero([pointer in row for row in self.Polylines])[0]) # This is a list of rows to explore
                    if verbose:
                        print('And the list of rows to explore is: {0}'.
                              format(rows))

                    while rows:
                        path_to_follow = rows.pop(0)

                        if verbose:
                            print('Following path: {0}'.format(path_to_follow))

                        row_avoid.append(path_to_follow)

                        tmp = list(self.Polylines[path_to_follow])
                        pointer = set.difference(set(tmp), set(vertex_avoid)).pop()
                        vertex_avoid.append(pointer)
                        if verbose:
                            print('pointer is now: {0}'.format(pointer))

                        if pointer == intersection[1]:
                            # Bingo!
                            if verbose:
                                print('Bingo! Both intersections are part '
                                      'of the same fundus!')
                            self.same_fundus = True
                            break

                        rows = rows + list(set.difference(set(np.nonzero([pointer in row for row in self.Polylines])[0]),set(row_avoid)))
                        if verbose:
                            print('Rows is now: {0}'.format(rows))

                self.label_boundary_segments[key] = np.append(segment,int(self.same_fundus))

        return self.label_boundary_segments

    def neighbors(self, vertex, verbose=False):
        """
        Take a vertex as input and return an array of the vertex's neighbors,
        as defined by self.Faces.
        """
        # First check to see if the neighbors matrix was constructed.
        if not self.found_neighbors:

            if verbose:
                print('Constructing neighborhood function.')

            self.Neighbors = lil_matrix((self.num_points, self.num_points))

            for row in self.Faces:
                self.Neighbors[row[0], row[1]] = self.Neighbors[row[1], row[0]] = 1
                self.Neighbors[row[0], row[2]] = self.Neighbors[row[2], row[0]] = 1
                self.Neighbors[row[1], row[2]] = self.Neighbors[row[2], row[1]] = 1

            self.Neighbors = self.Neighbors.tocsr()
            self.found_neighbors = 1

        return np.nonzero(self.Neighbors[vertex])[1]

    ##########################################################################
    # ------------------------------------------------------------------------
    #     Master methods
    # ------------------------------------------------------------------------
    ##########################################################################

    def realign_label_boundary(self, surface_file, polylines_file,
                               label_boundary_filename, output_file_regions,
                               output_file_boundaries, max_iters,
                               verbose=False):
        """
        Complete method to realign the label boundaries.
        Calls all necessary subroutines.

        Parameters
        ----------
        surface_file: string (vtk file containing the vertices,
                              faces and manual labels of brain surface)
        polylines_file: string (vtk file containing polylines)
        label_boundary_filename: string (vtk file with initial label boundaries)
        output_file_regions: string (vtk file to contain new labels)
        max_iters: int (maximum number of iterations)

        Returns
        -------
        output_file_regions: string (vtk file containing new labels)
        output_file_boundaries: string (vtk file containing highlighted
                                        label boundaries)

        """
        t0 = time()
        self.load_vtk_surface(surface_file)
        self.load_vtk_polylines(polylines_file)
        if verbose:
            print('Imported Data in: {0}'.format(time() - t0))

        self.initialize_seed_labels(init='label_boundary',
                                    output_filename = label_boundary_filename)
        self.find_label_boundary_segments()
        self.graph_based_learning(realign=True, max_iters=max_iters)
        self.assign_realigned_labels(filename = output_file_regions)
        self.find_label_boundary(realigned_labels=True,
                                   output_filename = output_file_boundaries)

        return output_file_regions, output_file_boundaries

    ##########################################################################
    # ------------------------------------------------------------------------
    #     Post-processing methods
    # ------------------------------------------------------------------------
    ##########################################################################

    def assign_realigned_labels(self, filename, threshold=0.5, verbose=False):
        """
        Reassign labels to vertices after realignment.

        Take the self.learned_matrix weight matrix
        and determine which vertices should be relabeled.
        Output a separate array, self.RLabels,
        that contains one label for each vertex.

        NOTE: assignment tricky -- based on label segments NOT labels

        Parameters
        ----------
        filename: string (output vtk file to visualize realigned label boundaries)

        Returns
        -------
        self.RLabels: array (of size num_points which contains the most
                             probable realigned label for each vertex)
        self.RLabels_file: string (vtk file containing most probable labels)

        """

        # First check that the label propagation algorithm has been called
        try:
            a = self.learned_matrix[0]
        except:
            if verbose:
                print('First call graph_based_learning().')
            return

        self.RLabels = self.Labels.copy()

        # Collect consensus labels...
        # self.get_consensus_labels(0,0)

        # Go column by column, find entries (indices to vertices)
        # which meet a criterion.  Store these indices as values in a
        # dictionary, with the index to the column (label segment) as the key.
        counter = 0
        vertices_to_change = {}
        for column in self.learned_matrix.T:
            vertices_to_change[counter] = list(np.nonzero(column > threshold)[0])
            #if set.intersection(set(vertices_to_change[counter]), set(self.consensus_labels)):
            #	if verbose:
                  # print('You are trying to cross consensus labels!'
            #	vertices_to_change[counter] = []
            # self.RLabels[vertices_to_change] = self.realignment_mapping[counter][1]
            counter += 1

        if verbose:
            print('There are {0} regions to be relabeled.'.format(
                len(vertices_to_change)))

        # Run check_for_polylines
        vertices_to_change = self.check_for_polylines(vertices_to_change)

        if verbose:
            print('After further checks, {0} regions are going to be relabeled.'.
                format(len(vertices_to_change)))

        # Resolve label ambiguities
        # vertices_to_change = self.resolve_label_ambiguity(vertices_to_change)

        # if verbose:
        #     print('After label resolving ambiguities, \
        # {} regions are to be relabeled:'.format(len(vertices_to_change))

        # Resolve label front ambiguities
        # vertices_to_change = self.resolve_label_front(vertices_to_change)

        # if verbose:
        #     print('After resolving label front ambiguities, \
        # {} regions are to be relabeled:'.format(len(vertices_to_change))

        for key, value in list(vertices_to_change.items()):
            if verbose:
                print('For key {0}, the following vertices will be changed: {1}'.
                      format(self.realignment_mapping[key],value))

        # For vertices that have passed all checks and are to be relabeled,
        # select the second (relabel) entry in the corresponding dictionary tuple
        for key, value in list(vertices_to_change.items()):
            self.RLabels[value] = self.realignment_mapping[key][1]

        # Write VTK file with the new labels
        self.RLabels_file = filename
        write_vtk(self.RLabels_file, self.Points, self.Vertices,
                  [], self.Faces, [self.RLabels], scalar_type='int')

        return self.RLabels, self.RLabels_file

    def check_for_polylines(self, dict_of_vertices, threshold=15,
                            verbose=False):
        """
        Check which groups of vertices contain a sufficient number that
        border polylines.

        Problem:
        A group of vertices is considered for label reassignment
        but they might not be situated between a label boundary
        and a polyline representing a label-delimiting feature
        (such as a fundus in the brain cortex).

        Solution:
        Here we check to see if the vertices include a sufficient
        number that flank polylines.  If the proposed label
        reassignment would not change many polyline vertices,
        then this is probably not a good group of vertices to change.
        (Remove them from the input dictionary.)

        Parameters
        ----------
        dict_of_vertices: dict
            key is the label index
            value is the list of vertices to be reassigned to the label

        threshold: int
            minimum number of vertices that must be part of polylines boundary

        Returns
        -------
        dict_of_vertices: dict (subset of the input dict)
            key is the label index
            value is the list of vertices to be reassigned to the label

        """
        self.find_polylines_flanks()

        if verbose:
            print(self.polylines_flanks_indices)

        for key, value in list(dict_of_vertices.items()):
            if len(np.intersect1d(value,self.polylines_flanks_indices)) < threshold:
                dict_of_vertices.pop(key)

        return dict_of_vertices

    def find_polylines_flanks(self):
        """
        Find vertices that flank the polylines.

        Returns
        -------
        self.polylines_flanks_indices: np array (of indices of vertices which border but are not part of the polylines)

        """

        self.polylines_flanks_indices = np.zeros(self.Labels.shape)

        for triangles in self.Faces:
            vertex0, vertex1, vertex2 = triangles[0], triangles[1], triangles[2]
            num_points_in_polylines = (vertex0 in self.polyline_elements) + \
                                        (vertex1 in self.polyline_elements) + (vertex2 in self.polyline_elements)
            if num_points_in_polylines > 0:
                self.polylines_flanks_indices[triangles] = 1
        self.polylines_flanks_indices[self.polyline_elements] = 0

        self.polylines_flanks_indices = np.nonzero(self.polylines_flanks_indices)[0]

        return self.polylines_flanks_indices

    def resolve_label_ambiguity(self, dict_of_vertices, verbose=False):
        """
        Resolve competing labels for two groups of vertices that have some
        vertices in common.

        Problem:
        Two groups of vertices are considered for label reassignment,
        each to a different label, but they share some vertices.
        To which label should they be reassigned?

        Solution:
        Here we simply select the larger group of vertices,
        and remove the other group from the input dictionary.
        Presumably this group runs parallel to the polyline
        representing a label-delimiting feature (such as a fundus).

        Parameters
        ----------
        dict_of_vertices: dict
            key is the label index
            value is the list of vertices to be reassigned to the label

        threshold: int
            minimum number of vertices that must be part of polylines boundary

        Returns
        -------
        dict_of_vertices: dict (subset of the input dict)
            key is the label index
            value is the list of vertices to be reassigned to the label

        """

        # First we need an efficient procedure to identify the overlaps.
        # Create a num_keys x num_keys matrix indicating whether there is overlap,
        # and if so, which array is larger.

        # As a proxy for whether a fundus runs along a label boundary,
        # we could simply see which label boundary segment relabels more vertices
        if verbose:
            print('We made it here! So far so good.')

        num_keys = len(dict_of_vertices)
        overlap = np.zeros((num_keys, num_keys))

        for key1, value1 in list(dict_of_vertices.items()):
            for key2, value2, in list(dict_of_vertices.items()):
                if key1 != key2:
                    overlap_problem = np.intersect1d(value1,value2).any()
                    if overlap_problem:
                        if verbose:
                            print('The keys are: {0} {1}'.format(key1, key2))
                        # 1 indicates value1 is larger, 2 if value2 is larger.
                        overlap[key1,key2] = (len(value2) > len(value1)) + 1
                        if verbose:
                            print('The value of overlap is: {0}'.
                                  format(overlap[key1,key2]))

                        # For NOW, let us disregard the matrix overlap and simply resolve the issue right here.
                        # In the future, this matrix may be desirable.

                        if overlap[key1, key2] == 1:
                            try:
                                dict_of_vertices.pop(key2)
                            except:
                                pass
                        else:
                            try:
                                dict_of_vertices.pop(key1)
                            except:
                                pass

        return dict_of_vertices

    def resolve_label_front(self, dict_of_vertices, verbose=False):
        """
        Resolve which group of vertices on either side of a label boundary
        is to be relabeled.

        The label boundary consists of a set of vertices on the surface mesh
        two thick, following from the (>=2 label neighborhood) definition,
        so we call these "label cosegments".

        In the case where a label boundary corresponds to a label-delimiting
        feature (such as a fundus) to one side of the label boundary,
        we wish to identify which cosegment should be used to relabel vertices
        between the label boundary and the feature.

        NOTE: face validity but does what we want it to do?

        Problem:
        Two groups of vertices are considered for label reassignment,
        but they lie to either side of a label boundary.
        If both are relabeled, this can result in a "striping" artifact.
        Which of the two groups of vertices should be relabeled?

        Solution:
        Here we simply select the larger group of vertices,
        and remove the other group from the input dictionary.
        Presumably this group runs parallel to the polyline
        representing a label-delimiting feature (such as a fundus).

        Parameters
        ----------
        dict_of_vertices: dict
            key is the label index
            value is the list of vertices considered for reassignment to label

        threshold: int
            minimum number of vertices that must be part of polylines boundary

        Returns
        -------
        dict_of_vertices: dict (subset of the input dict)
            key:  label index
            value:  list of vertices to be reassigned to the label

        """

        for key1, value1 in list(dict_of_vertices.items()):
            for key2, value2 in list(dict_of_vertices.items()):
                if key1 != key2:
                    # If thay are co-segments...
                    if len(np.intersect1d(self.realignment_mapping[key1], self.realignment_mapping[key2])) == 2:
                        if verbose:
                            print('Found co-segments.')
                        # Find which array contains more polylines border vertices...
                        border1 = len(np.intersect1d(value1,self.polylines_flanks_indices))
                        border2 = len(np.intersect1d(value2,self.polylines_flanks_indices))
                        if border1 > border2:
                            try:
                                dict_of_vertices.pop(key2)
                            except:
                                pass
                        else:
                            try:
                                dict_of_vertices.pop(key1)
                            except:
                                pass

        return dict_of_vertices


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)()