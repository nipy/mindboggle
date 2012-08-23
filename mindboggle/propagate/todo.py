
def jacobi_iteration(self, alpha, max_iters, tol, eps):
    """
    Perform label propagation inspired from Jacobi iterative algorithm
    to propagate labels to unlabeled vertices.
    Features: Soft label clamps (alpha), probabilistic solution.
    See: Chapelle, ch. 11 (algorithm 11.2).

    alpha: float (clamping factor)
    diagonal: int (value for diagonal entries in weight matrix, 0 or 1)
    eps: float (for numerical stability in some algorithms)

    """

    """The next approach to be considered in the semi-supervised learning case
    is to propagate labels on the graph, using a modified version of the above algorithm.
    The main differences are soft clamping, forcing the diagonals to be equal to 0,
    and the introduction of an error term (eps) for numerical stability.

    We start with a set of n vertices, l of which are labeled, and u unlabeled.
    The algorithm takes as its input the affinity matrix W (self.affinity_matrix).
    From the affinity matrix, we construct the diagonal degree matrix,
    which is a measure of the total weight (or number of edges) which are attached to a vertex."""

    self.DDM = go.compute_diagonal_degree_matrix(self.affinity_matrix, inverse=True)

    """Next, we must initialize a vector to represent the results
    of the label propagation algorithm. It will contain l labels and u 0's.
    This has already been done by the function initialize_seed_labels(),
    and is called self.seed_labels.
    We will just check to make sure this has been accomplished."""

    if isinstance(self.seed_labels,int):
        print('Please initialize the labels by calling self.initialize_seed_labels()')
        return

    """ Now, we can actually proceed to perform the iterative algorithm.
    At each timestep, the labels will be updated to reflect the weighted average
    of adjacent vertices. An important caveat of this algorithm,
    is that the labeled vertices do not remain fixed, or clamped.
    The algorithm repeates itself until either convergence or max_iters
    (which will prevent excessive computation time).
    We must again take care to solve the multi-label problem.
    So, to begin, let us first construct this probabilistic label assignment:
    This matrix will store a 1 for 100% probability, 0 for 0%, and fractional values for the rest.
    We will rename self.label_matrix for this purpose.
    We will later change the -1s to 0s."""

    self.probabilistic_assignment = self.label_matrix

    """ Before proceeding, let us check that the parameters are valid"""

    if not (alpha < 1 and alpha > 0 and eps > 0 and isinstance(max_iters, int) and max_iters > 0 and tol > 0):
        print('You have failed to properly input parameters. Alpha must be strictly between 0 and 1, eps \ \
                and tol must be greater than 0 and max_iters must be an integer.')
        return

    """ As vertices get labeled, we assign a confidence measure to the labeling and store the value in this matrix.
    Now, let us go column by column, and run the weighted averaging algorithm.
    For each column, you're studying one label. Therefore, when updating self.probabilistic_assignment,
    you'll be working with one column at a time too.
    If a label gets vertex, keep the fractional value, do not simply round to 1 to assign membership."""

    i = 0 # record of label number

    for column in self.probabilistic_assignment.T:
        self.labeled_indices = column[self.seed_labels > 0]
        column = csr_matrix(column).transpose()
        converged = False
        counter = 0
        while not converged and counter < max_iters:
            tmp = self.DDM * self.affinity_matrix * column # column matrix
            tmp = tmp.tolil() # store results of iteration
            tmp[self.labeled_indices,0] = self.labeled_indices # reset
            converged = (np.abs(column - tmp).sum() < tol) # check convergence
            print('convergence=', np.abs(column-tmp).sum())
            column = csr_matrix(tmp)
            counter += 1

        # Print out the number of iterations, so that we get a sense for future runs.
        # It is also an indication of whether the algorithm converged.

        if counter == max_iters:
            print('The algorithm did not converge.')
        else:
            print('The algorithm converged in {0} iterations.'.format(str(counter)))
        i += 1

    """ Before reporting the probabilistic assignment, we change all -1's, indicating
    0 probability that the vertex has that label."""

    self.probabilistic_assignment[self.probabilistic_assignment==-1] = 0

    return self.probabilistic_assignment

