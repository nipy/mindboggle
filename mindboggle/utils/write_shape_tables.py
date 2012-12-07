"""
Save a shape table for the fundus label boundary vertices
in each subject hemisphere.
"""
import os
from mindboggle.utils.io_vtk import load_scalars, write_vertex_shapes_table, \
    write_mean_shapes_table
from mindboggle.utils.mesh_operations import find_neighbors, detect_boundaries

do_fundus = False

# data_path = os.environ['MINDBOGGLE_DATA']
labeled_subjects_path = '/hd2/Lab/Brains/Mindboggle101/subjects'
data_path = '/Users/arno/Desktop/output_measures/results'
column_names = ['labels', 'area', 'depth', 'mean_curvature',
                'gauss_curvature', 'max_curvature', 'min_curvature',
                'freesurfer_thickness', 'freesurfer_convexity']
hemis = ['lh','rh']
list_file = os.path.join(os.environ['MINDBOGGLE'], 'info', 'atlases101.txt')
exclude_values = [-1]

# For each subject in the subjects file
fid = open(list_file, 'r')
subjects = fid.readlines()
subjects = [''.join(x.split()) for x in subjects]
for subject in subjects: #['OASIS-TRT-20-11']: #subjects[90::]:
    print(subject)
    for hemi in hemis:

            # Load labels, label boundaries (corr. to fundi), sulci, and fundi
            label_file = os.path.join(labeled_subjects_path, subject,
                                      'label', hemi+'.labels.DKT25.manual.vtk')
            points, faces, labels, n_vertices = load_scalars(label_file, True)
            neighbor_lists = find_neighbors(faces, n_vertices)
            sulcus_file = os.path.join(data_path, 'features',
                                      '_hemi_'+hemi+'_subject_' + subject, 'sulci.vtk')
            points, faces, sulci, n_vertices = load_scalars(sulcus_file, True)
            sulcus_indices = [i for i,x in enumerate(sulci) if x > -1]
            boundary_indices, label_pairs, foo = detect_boundaries(sulcus_indices, labels,
                                                neighbor_lists)
            if do_fundus:
                fundus_file = os.path.join(data_path, 'features',
                              '_hemi_'+hemi+'_subject_' + subject, 'fundi.vtk')
                points, faces, fundi, n_vertices = load_scalars(fundus_file, True)
                fundus_indices = [i for i,x in enumerate(fundi) if x > -1]

            # Shape measure files
            area_file = os.path.join(data_path, 'measures',
                '_hemi_'+hemi+'_subject_'+subject, hemi+'.pial.area.vtk')
            depth_file = os.path.join(data_path, 'measures',
                '_hemi_'+hemi+'_subject_'+subject, hemi+'.pial.depth.vtk')
            mean_curvature_file = os.path.join(data_path, 'measures',
                '_hemi_'+hemi+'_subject_'+subject, hemi+'.pial.curv.avg.vtk')
            gauss_curvature_file = os.path.join(data_path, 'measures',
                '_hemi_'+hemi+'_subject_'+subject, hemi+'.pial.curv.gauss.vtk')
            max_curvature_file = os.path.join(data_path, 'measures',
                '_hemi_'+hemi+'_subject_'+subject, hemi+'.pial.curv.max.vtk')
            min_curvature_file = os.path.join(data_path, 'measures',
                '_hemi_'+hemi+'_subject_'+subject, hemi+'.pial.curv.min.vtk')
            thickness_file = os.path.join(data_path, 'measures',
                '_hemi_'+hemi+'_subject_'+subject, 'thickness.vtk')
            convexity_file = os.path.join(data_path, 'measures',
                '_hemi_'+hemi+'_subject_'+subject, 'sulc.vtk')

            # Write label shape table
            print('  Write vertex label shape table for subject {0}'.format(subject))
            vtx_label_file = subject + '_' + hemi + '_vertex_label_shape_table.txt'
            write_vertex_shapes_table(vtx_label_file, column_names,
                range(len(points)), area_file, depth_file, mean_curvature_file,
                gauss_curvature_file, max_curvature_file, min_curvature_file,
                thickness_file, convexity_file, segments=labels, nonsegments=[-1])
            print('  Write average label shape table for subject {0}'.format(subject))
            avg_label_file = subject + '_' + hemi + '_avg_label_shape_table.txt'
            write_mean_shapes_table(avg_label_file, column_names, labels,
                                    exclude_values, area_file, depth_file,
                                    mean_curvature_file, gauss_curvature_file,
                                    max_curvature_file, min_curvature_file,
                                    thickness_file='', convexity_file='')

            # Write sulcus shape table
            print('  Write vertex sulcus shape table for subject {0}'.format(subject))
            vtx_sulcus_table = subject + '_' + hemi + '_vertex_sulcus_shape_table.txt'
            write_vertex_shapes_table(vtx_sulcus_table, column_names,
                sulcus_indices, area_file, depth_file, mean_curvature_file,
                gauss_curvature_file, max_curvature_file, min_curvature_file,
                thickness_file, convexity_file, segments=sulci, nonsegments=[-1])
            print('  Write average sulcus shape table for subject {0}'.format(subject))
            avg_sulcus_file = subject + '_' + hemi + '_avg_sulcus_shape_table.txt'
            write_mean_shapes_table(avg_sulcus_file, column_names, labels,
                                    exclude_values, area_file, depth_file,
                                    mean_curvature_file, gauss_curvature_file,
                                    max_curvature_file, min_curvature_file,
                                    thickness_file='', convexity_file='')

            # Write fundus label boundary shape table
            print('  Write vertex fundus label boundary shape table for subject {0}'.format(subject))
            vtx_label_boundary_table = subject + '_' + hemi + '_vertex_fundus_label_boundary_shape_table.txt'
            write_vertex_shapes_table(vtx_label_boundary_table, column_names,
                boundary_indices, area_file, depth_file, mean_curvature_file,
                gauss_curvature_file, max_curvature_file, min_curvature_file,
                thickness_file, convexity_file, segments=sulci, nonsegments=[-1])
            avg_label_boundary_file = subject + '_' + hemi + '_avg_fundus_label_boundary_shape_table.txt'
            print('  Write average fundus label boundary shape table for subject {0}'.format(subject))
            write_mean_shapes_table(avg_label_boundary_file, column_names, labels,
                                    exclude_values, area_file, depth_file,
                                    mean_curvature_file, gauss_curvature_file,
                                    max_curvature_file, min_curvature_file,
                                    thickness_file='', convexity_file='')

            # Write fundus shape table
            if do_fundus:
                print('  Write vertex fundus shape table for subject {0}'.format(subject))
                vtx_fundus_table = subject + '_' + hemi + '_vertex_fundus_shape_table.txt'
                write_vertex_shapes_table(vtx_fundus_table, column_names,
                    fundus_indices, area_file, depth_file, mean_curvature_file,
                    gauss_curvature_file, max_curvature_file, min_curvature_file,
                    thickness_file, convexity_file, segments=fundi, nonsegment_IDs=[-1])
                print('  Write average fundus shape table for subject {0}'.format(subject))
                avg_fundus_table = subject + '_' + hemi + '_avg_fundus_shape_table.txt'
                write_mean_shapes_table(avg_fundus_table, column_names, labels,
                                        exclude_values, area_file, depth_file,
                                        mean_curvature_file, gauss_curvature_file,
                                        max_curvature_file, min_curvature_file,
                                        thickness_file='', convexity_file='')
