#ifndef __verdict_mangle_h
#define __verdict_mangle_h

/*
This header file mangles all symbols exported from the verdict library.
It is included in all files while building the verdict library.

The following command was used to obtain the symbol list:

nm libverdict.a |grep " [TRD] "

This is the way to recreate the whole list:

nm libvtkverdict.dylib  |grep " [TRD] " | sed -e 's+^[^_]*_++g' | \
   awk '{ print "#define "$1" vtk_"$1 }' | \
   grep -v __fini | grep -v __init | sort > verdict_mangle.h.in

After this, you must manually edit the file and replace any names
mangled by the C++ compiler with their generic counterparts since
different C++ compilers have different mangling schemes.

Note that _fini and _init should be excluded because they are not functions
implemented by the library but are rather created by the linker and
used when the shared library is loaded/unloaded from an executable.

*/

/* C++ function and class names */
#define v_determinant vtk_v_determinant
#define v_interpolate vtk_v_interpolate
#define v_vector_rotate vtk_v_vector_rotate
#define VerdictVector vtk_VerdictVector
#define GaussIntegration vtk_GaussIntegration

/* C function names */
#define v_edge_length vtk_v_edge_length
#define v_edge_quality vtk_v_edge_quality
#define v_hex_condition vtk_v_hex_condition
#define v_hex_diagonal vtk_v_hex_diagonal
#define v_hex_dimension vtk_v_hex_dimension
#define v_hex_distortion vtk_v_hex_distortion
#define v_hex_edge_ratio vtk_v_hex_edge_ratio
#define v_hex_jacobian vtk_v_hex_jacobian
#define v_hex_max_aspect_frobenius vtk_v_hex_max_aspect_frobenius
#define v_hex_max_edge_ratio vtk_v_hex_max_edge_ratio
#define v_hex_med_aspect_frobenius vtk_v_hex_med_aspect_frobenius
#define v_hex_oddy vtk_v_hex_oddy
#define v_hex_quality vtk_v_hex_quality
#define v_hex_relative_size_squared vtk_v_hex_relative_size_squared
#define v_hex_scaled_jacobian vtk_v_hex_scaled_jacobian
#define v_hex_shape vtk_v_hex_shape
#define v_hex_shape_and_size vtk_v_hex_shape_and_size
#define v_hex_shear vtk_v_hex_shear
#define v_hex_shear_and_size vtk_v_hex_shear_and_size
#define v_hex_skew vtk_v_hex_skew
#define v_hex_stretch vtk_v_hex_stretch
#define v_hex_taper vtk_v_hex_taper
#define v_hex_volume vtk_v_hex_volume
#define v_knife_quality vtk_v_knife_quality
#define v_knife_volume vtk_v_knife_volume
#define v_pyramid_quality vtk_v_pyramid_quality
#define v_pyramid_volume vtk_v_pyramid_volume
#define v_quad_area vtk_v_quad_area
#define v_quad_aspect_ratio vtk_v_quad_aspect_ratio
#define v_quad_condition vtk_v_quad_condition
#define v_quad_distortion vtk_v_quad_distortion
#define v_quad_edge_ratio vtk_v_quad_edge_ratio
#define v_quad_jacobian vtk_v_quad_jacobian
#define v_quad_max_aspect_frobenius vtk_v_quad_max_aspect_frobenius
#define v_quad_max_edge_ratio vtk_v_quad_max_edge_ratio
#define v_quad_maximum_angle vtk_v_quad_maximum_angle
#define v_quad_med_aspect_frobenius vtk_v_quad_med_aspect_frobenius
#define v_quad_minimum_angle vtk_v_quad_minimum_angle
#define v_quad_oddy vtk_v_quad_oddy
#define v_quad_quality vtk_v_quad_quality
#define v_quad_radius_ratio vtk_v_quad_radius_ratio
#define v_quad_relative_size_squared vtk_v_quad_relative_size_squared
#define v_quad_scaled_jacobian vtk_v_quad_scaled_jacobian
#define v_quad_shape vtk_v_quad_shape
#define v_quad_shape_and_size vtk_v_quad_shape_and_size
#define v_quad_shear vtk_v_quad_shear
#define v_quad_shear_and_size vtk_v_quad_shear_and_size
#define v_quad_skew vtk_v_quad_skew
#define v_quad_stretch vtk_v_quad_stretch
#define v_quad_taper vtk_v_quad_taper
#define v_quad_warpage vtk_v_quad_warpage
#define v_set_hex_size vtk_v_set_hex_size
#define v_set_quad_size vtk_v_set_quad_size
#define v_set_tet_size vtk_v_set_tet_size
#define v_set_tri_normal_func vtk_v_set_tri_normal_func
#define v_set_tri_size vtk_v_set_tri_size
#define v_tet_aspect_beta vtk_v_tet_aspect_beta
#define v_tet_aspect_frobenius vtk_v_tet_aspect_frobenius
#define v_tet_aspect_gamma vtk_v_tet_aspect_gamma
#define v_tet_aspect_ratio vtk_v_tet_aspect_ratio
#define v_tet_collapse_ratio vtk_v_tet_collapse_ratio
#define v_tet_condition vtk_v_tet_condition
#define v_tet_distortion vtk_v_tet_distortion
#define v_tet_edge_ratio vtk_v_tet_edge_ratio
#define v_tet_jacobian vtk_v_tet_jacobian
#define v_tet_minimum_angle vtk_v_tet_minimum_angle
#define v_tet_quality vtk_v_tet_quality
#define v_tet_radius_ratio vtk_v_tet_radius_ratio
#define v_tet_relative_size_squared vtk_v_tet_relative_size_squared
#define v_tet_scaled_jacobian vtk_v_tet_scaled_jacobian
#define v_tet_shape vtk_v_tet_shape
#define v_tet_shape_and_size vtk_v_tet_shape_and_size
#define v_tet_volume vtk_v_tet_volume
#define v_tri_area vtk_v_tri_area
#define v_tri_aspect_frobenius vtk_v_tri_aspect_frobenius
#define v_tri_aspect_ratio vtk_v_tri_aspect_ratio
#define v_tri_condition vtk_v_tri_condition
#define v_tri_distortion vtk_v_tri_distortion
#define v_tri_edge_ratio vtk_v_tri_edge_ratio
#define v_tri_maximum_angle vtk_v_tri_maximum_angle
#define v_tri_minimum_angle vtk_v_tri_minimum_angle
#define v_tri_quality vtk_v_tri_quality
#define v_tri_radius_ratio vtk_v_tri_radius_ratio
#define v_tri_relative_size_squared vtk_v_tri_relative_size_squared
#define v_tri_scaled_jacobian vtk_v_tri_scaled_jacobian
#define v_tri_shape vtk_v_tri_shape
#define v_tri_shape_and_size vtk_v_tri_shape_and_size
#define v_wedge_quality vtk_v_wedge_quality
#define v_wedge_volume vtk_v_wedge_volume

#endif /* __verdict_mangle_h */
