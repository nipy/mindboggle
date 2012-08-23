
##############################################################################
#
#   Alternative to weighted average approach to label propagation
#
##############################################################################
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

##############################################################################
#
#   Use ANTS ImageMath rather than FreeSurfer to fill a volume:
#
##############################################################################

# pipeline:

    """
    #-------------------------------------------------------------------------
    # Put surface vertices in a volume
    #-------------------------------------------------------------------------
    surf2vol = node(name='Surface_to_volume',
                    interface = fn(function = surface_to_volume,
                                   input_names = ['surface_file',
                                                  'volume_file',
                                                  'use_freesurfer'],
                                   output_names = ['output_file']))
    surf2vol.inputs.use_freesurfer = use_freesurfer
    atlasflow.add_nodes([surf2vol])
    atlasflow.connect([(vote, surf2vol, [('maxlabel_file','surface_file')])])
    if use_freesurfer:
        mbflow.connect([(vol, atlasflow,
                         [('volume_file','Surface_to_volume.volume_file')])])
    else:
        mbflow.connect([(vol, atlasflow,
                         [('volume_file','Surface_to_volume.volume_file')])])
    #-------------------------------------------------------------------------
    # Fill volume mask with surface vertex labels
    #-------------------------------------------------------------------------
    fill_maxlabels = node(name='Fill_volume_maxlabels',
                          interface = fn(function = fill_volume,
                                         input_names = ['command',
                                                        'input_file',
                                                        'mask_file'],
                                         output_names = ['output_file']))
    fill_maxlabels.inputs.command = imagemath
    atlasflow.add_nodes([fill_maxlabels])
    atlasflow.connect([(surf2vol, fill_maxlabels,
                        [('output_file', 'input_file')])])
    if use_freesurfer:
        mbflow.connect([(convertvol, atlasflow,
                         [('out_file','Fill_volume_maxlabels.mask_file')])])
    else:
        mbflow.connect([(vol, atlasflow,
                         [('volume_file','Fill_volume_maxlabels.mask_file')])])
    mbflow.connect([(atlasflow, datasink,
                     [('Fill_volume_maxlabels.output_file',
                       'labels.@maxvolume')])])
    """

# volume_functions.py:

def surface_to_volume(surface_file, volume_file, use_freesurfer):
    """
    Save the vertices of a FreeSurfer surface mesh as an image volume.

    FIX:  labels are incorrect!
    """

    from os import path, getcwd, error
    import numpy as np
    import nibabel as nb
    import vtk

    scalar_name = "Max_(majority_labels)"
    if use_freesurfer:
        trans = 128  # translation to middle of FreeSurfer conformed space

    # Load image volume
    vol = nb.load(volume_file)
    vol_shape = vol.shape
    xfm = vol.get_affine()

    # Load surface
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(surface_file)
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    d = data.GetPointData()
    labels = d.GetArray(scalar_name)

    # Create a new volume (permuted and flipped)
    V = np.zeros(vol_shape)
    npoints = data.GetNumberOfPoints()
    for i in range(npoints):
        point = data.GetPoint(i)
        label = labels.GetValue(i)
        if use_freesurfer:
            V[-point[0]+trans, -point[2]+trans, point[1]+trans] = label
        else:
            V[point[0], point[1], point[2]] = label

    # Save the image with the same affine transform
    output_file = path.join(getcwd(), surface_file.strip('.vtk')+'.nii.gz')
    img = nb.Nifti1Image(V, xfm)
    img.to_filename(output_file)

    return output_file

    """
    # Alternative (NEEDS A FIX):
    # Create a new volume (permuted and flipped)
    from apply_utils import apply_affine
    xfm2 = np.array([[-1,0,0,128],
                    [0,0,-1,128],
                    [0,1,0,128],
                    [0,0,0,1]],dtype=float)
    xyz = apply_affine(xyz[:,0], xyz[:,1], xyz[:,2], xfm2)

    V = np.zeros(vol_shape)
    for vertex in xyz:
        V[vertex[0], vertex[1], vertex[2]] = 1
    """

def fill_volume(command, input_file, mask_file):
    """
    Fill (e.g., gray matter) volume with surface labels using ANTS
    (ImageMath's PropagateLabelsThroughMask)

    Brian avants:
    The initial box labels are propagated through the gray matter with
    gm-probability dependent speed. It uses the fast marching algorithm.
    You can control how tightly the propagation follows the gray matter
    label by adjusting the speed image -- e.g. a binary speed image
    will constrain the propagated label only to the gm.

    """
    from os import path, getcwd, error
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')

    print("Fill gray matter volume with surface labels using ANTS...")

    output_file = path.join(getcwd(), input_file.strip('.nii.gz')+'.fill.nii.gz')

    args = ['3',
            output_file,
            'PropagateLabelsThroughMask',
            mask_file,
            input_file]

    cli = CommandLine(command=command)
    cli.inputs.args = ' '.join(args)
    logger.info(cli.cmdline)
    cli.run()

    return output_file
