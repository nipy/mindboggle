#ifndef vtk_png_mangle_h
#define vtk_png_mangle_h

/*

This header file mangles all symbols exported from the png library.
It is included in all files while building the png library.  Due to
namespace pollution, no png headers should be included in .h files in
VTK.

The following command was used to obtain the symbol list:

nm libvtkpng.a |grep " [TR] "

*/

#define png_access_version_number vtk_png_access_version_number
#define png_build_gamma_table vtk_png_build_gamma_table
#define png_build_grayscale_palette vtk_png_build_grayscale_palette
#define png_calculate_crc vtk_png_calculate_crc
#define png_check_chunk_name vtk_png_check_chunk_name
#define png_check_keyword vtk_png_check_keyword
#define png_check_sig vtk_png_check_sig
#define png_chunk_error vtk_png_chunk_error
#define png_chunk_warning vtk_png_chunk_warning
#define png_combine_row vtk_png_combine_row
#define png_convert_from_struct_tm vtk_png_convert_from_struct_tm
#define png_convert_from_time_t vtk_png_convert_from_time_t
#define png_convert_to_rfc1123 vtk_png_convert_to_rfc1123
#define png_crc_error vtk_png_crc_error
#define png_crc_finish vtk_png_crc_finish
#define png_crc_read vtk_png_crc_read
#define png_create_info_struct vtk_png_create_info_struct
#define png_create_read_struct vtk_png_create_read_struct
#define png_create_struct vtk_png_create_struct
#define png_create_write_struct vtk_png_create_write_struct
#define png_data_freer vtk_png_data_freer
#define png_decompress_chunk vtk_png_decompress_chunk
#define png_destroy_info_struct vtk_png_destroy_info_struct
#define png_destroy_read_struct vtk_png_destroy_read_struct
#define png_destroy_struct vtk_png_destroy_struct
#define png_destroy_write_struct vtk_png_destroy_write_struct
#define png_do_background vtk_png_do_background
#define png_do_bgr vtk_png_do_bgr
#define png_do_chop vtk_png_do_chop
#define png_do_dither vtk_png_do_dither
#define png_do_expand vtk_png_do_expand
#define png_do_expand_palette vtk_png_do_expand_palette
#define png_do_gamma vtk_png_do_gamma
#define png_do_gray_to_rgb vtk_png_do_gray_to_rgb
#define png_do_invert vtk_png_do_invert
#define png_do_pack vtk_png_do_pack
#define png_do_packswap vtk_png_do_packswap
#define png_do_read_filler vtk_png_do_read_filler
#define png_do_read_interlace vtk_png_do_read_interlace
#define png_do_read_invert_alpha vtk_png_do_read_invert_alpha
#define png_do_read_swap_alpha vtk_png_do_read_swap_alpha
#define png_do_read_transformations vtk_png_do_read_transformations
#define png_do_rgb_to_gray vtk_png_do_rgb_to_gray
#define png_do_shift vtk_png_do_shift
#define png_do_strip_filler vtk_png_do_strip_filler
#define png_do_swap vtk_png_do_swap
#define png_do_unpack vtk_png_do_unpack
#define png_do_unshift vtk_png_do_unshift
#define png_do_write_interlace vtk_png_do_write_interlace
#define png_do_write_invert_alpha vtk_png_do_write_invert_alpha
#define png_do_write_swap_alpha vtk_png_do_write_swap_alpha
#define png_do_write_transformations vtk_png_do_write_transformations
#define png_error vtk_png_error
#define png_flush vtk_png_flush
#define png_free vtk_png_free
#define png_free_data vtk_png_free_data
#define png_get_IHDR vtk_png_get_IHDR
#define png_get_PLTE vtk_png_get_PLTE
#define png_get_bKGD vtk_png_get_bKGD
#define png_get_bit_depth vtk_png_get_bit_depth
#define png_get_cHRM vtk_png_get_cHRM
#define png_get_cHRM_fixed vtk_png_get_cHRM_fixed
#define png_get_channels vtk_png_get_channels
#define png_get_color_type vtk_png_get_color_type
#define png_get_compression_buffer_size vtk_png_get_compression_buffer_size
#define png_get_compression_type vtk_png_get_compression_type
#define png_get_copyright vtk_png_get_copyright
#define png_get_error_ptr vtk_png_get_error_ptr
#define png_get_filter_type vtk_png_get_filter_type
#define png_get_gAMA vtk_png_get_gAMA
#define png_get_gAMA_fixed vtk_png_get_gAMA_fixed
#define png_get_hIST vtk_png_get_hIST
#define png_get_header_ver vtk_png_get_header_ver
#define png_get_header_version vtk_png_get_header_version
#define png_get_iCCP vtk_png_get_iCCP
#define png_get_image_height vtk_png_get_image_height
#define png_get_image_width vtk_png_get_image_width
#define png_get_int_32 vtk_png_get_int_32
#define png_get_interlace_type vtk_png_get_interlace_type
#define png_get_io_ptr vtk_png_get_io_ptr
#define png_get_libpng_ver vtk_png_get_libpng_ver
#define png_get_oFFs vtk_png_get_oFFs
#define png_get_pCAL vtk_png_get_pCAL
#define png_get_pHYs vtk_png_get_pHYs
#define png_get_pixel_aspect_ratio vtk_png_get_pixel_aspect_ratio
#define png_get_pixels_per_meter vtk_png_get_pixels_per_meter
#define png_get_progressive_ptr vtk_png_get_progressive_ptr
#define png_get_rgb_to_gray_status vtk_png_get_rgb_to_gray_status
#define png_get_rowbytes vtk_png_get_rowbytes
#define png_get_rows vtk_png_get_rows
#define png_get_sBIT vtk_png_get_sBIT
#define png_get_sCAL vtk_png_get_sCAL
#define png_get_sPLT vtk_png_get_sPLT
#define png_get_sRGB vtk_png_get_sRGB
#define png_get_signature vtk_png_get_signature
#define png_get_tIME vtk_png_get_tIME
#define png_get_tRNS vtk_png_get_tRNS
#define png_get_text vtk_png_get_text
#define png_get_uint_16 vtk_png_get_uint_16
#define png_get_uint_32 vtk_png_get_uint_32
#define png_get_unknown_chunks vtk_png_get_unknown_chunks
#define png_get_user_chunk_ptr vtk_png_get_user_chunk_ptr
#define png_get_user_transform_ptr vtk_png_get_user_transform_ptr
#define png_get_valid vtk_png_get_valid
#define png_get_x_offset_microns vtk_png_get_x_offset_microns
#define png_get_x_offset_pixels vtk_png_get_x_offset_pixels
#define png_get_x_pixels_per_meter vtk_png_get_x_pixels_per_meter
#define png_get_y_offset_microns vtk_png_get_y_offset_microns
#define png_get_y_offset_pixels vtk_png_get_y_offset_pixels
#define png_get_y_pixels_per_meter vtk_png_get_y_pixels_per_meter
#define png_handle_IEND vtk_png_handle_IEND
#define png_handle_IHDR vtk_png_handle_IHDR
#define png_handle_PLTE vtk_png_handle_PLTE
#define png_handle_as_unknown vtk_png_handle_as_unknown
#define png_handle_bKGD vtk_png_handle_bKGD
#define png_handle_cHRM vtk_png_handle_cHRM
#define png_handle_gAMA vtk_png_handle_gAMA
#define png_handle_hIST vtk_png_handle_hIST
#define png_handle_iCCP vtk_png_handle_iCCP
#define png_handle_oFFs vtk_png_handle_oFFs
#define png_handle_pCAL vtk_png_handle_pCAL
#define png_handle_pHYs vtk_png_handle_pHYs
#define png_handle_sBIT vtk_png_handle_sBIT
#define png_handle_sCAL vtk_png_handle_sCAL
#define png_handle_sPLT vtk_png_handle_sPLT
#define png_handle_sRGB vtk_png_handle_sRGB
#define png_handle_tEXt vtk_png_handle_tEXt
#define png_handle_tIME vtk_png_handle_tIME
#define png_handle_tRNS vtk_png_handle_tRNS
#define png_handle_unknown vtk_png_handle_unknown
#define png_handle_zTXt vtk_png_handle_zTXt
#define png_info_destroy vtk_png_info_destroy
#define png_info_init vtk_png_info_init
#define png_info_init_3 vtk_png_info_init_3
#define png_init_io vtk_png_init_io
#define png_init_read_transformations vtk_png_init_read_transformations
#define png_libpng_ver vtk_png_libpng_ver
#define png_malloc vtk_png_malloc
#define png_memcpy_check vtk_png_memcpy_check
#define png_memset_check vtk_png_memset_check
#define png_mmx_support vtk_png_mmx_support
#define png_pass_dsp_mask vtk_png_pass_dsp_mask
#define png_pass_mask vtk_png_pass_mask
#define png_pass_yinc vtk_png_pass_yinc
#define png_pass_ystart vtk_png_pass_ystart
#define png_pass_inc vtk_png_pass_inc
#define png_pass_start vtk_png_pass_start
#define png_permit_empty_plte vtk_png_permit_empty_plte
#define png_process_IDAT_data vtk_png_process_IDAT_data
#define png_process_data vtk_png_process_data
#define png_process_some_data vtk_png_process_some_data
#define png_progressive_combine_row vtk_png_progressive_combine_row
#define png_push_crc_finish vtk_png_push_crc_finish
#define png_push_crc_skip vtk_png_push_crc_skip
#define png_push_fill_buffer vtk_png_push_fill_buffer
#define png_push_handle_tEXt vtk_png_push_handle_tEXt
#define png_push_handle_unknown vtk_png_push_handle_unknown
#define png_push_handle_zTXt vtk_png_push_handle_zTXt
#define png_push_have_end vtk_png_push_have_end
#define png_push_have_info vtk_png_push_have_info
#define png_push_have_row vtk_png_push_have_row
#define png_push_process_row vtk_png_push_process_row
#define png_push_read_IDAT vtk_png_push_read_IDAT
#define png_push_read_chunk vtk_png_push_read_chunk
#define png_push_read_sig vtk_png_push_read_sig
#define png_push_read_tEXt vtk_png_push_read_tEXt
#define png_push_read_zTXt vtk_png_push_read_zTXt
#define png_push_restore_buffer vtk_png_push_restore_buffer
#define png_push_save_buffer vtk_png_push_save_buffer
#define png_read_data vtk_png_read_data
#define png_read_destroy vtk_png_read_destroy
#define png_read_end vtk_png_read_end
#define png_read_filter_row vtk_png_read_filter_row
#define png_read_finish_row vtk_png_read_finish_row
#define png_read_image vtk_png_read_image
#define png_read_info vtk_png_read_info
#define png_read_init vtk_png_read_init
#define png_read_init_2 vtk_png_read_init_2
#define png_read_init_3 vtk_png_read_init_3
#define png_read_png vtk_png_read_png
#define png_read_push_finish_row vtk_png_read_push_finish_row
#define png_read_row vtk_png_read_row
#define png_read_rows vtk_png_read_rows
#define png_read_start_row vtk_png_read_start_row
#define png_read_transform_info vtk_png_read_transform_info
#define png_read_update_info vtk_png_read_update_info
#define png_reset_crc vtk_png_reset_crc
#define png_reset_zstream vtk_png_reset_zstream
#define png_save_int_32 vtk_png_save_int_32
#define png_save_uint_16 vtk_png_save_uint_16
#define png_save_uint_32 vtk_png_save_uint_32
#define png_set_IHDR vtk_png_set_IHDR
#define png_set_PLTE vtk_png_set_PLTE
#define png_set_bKGD vtk_png_set_bKGD
#define png_set_background vtk_png_set_background
#define png_set_bgr vtk_png_set_bgr
#define png_set_cHRM vtk_png_set_cHRM
#define png_set_cHRM_fixed vtk_png_set_cHRM_fixed
#define png_set_compression_buffer_size vtk_png_set_compression_buffer_size
#define png_set_compression_level vtk_png_set_compression_level
#define png_set_compression_mem_level vtk_png_set_compression_mem_level
#define png_set_compression_method vtk_png_set_compression_method
#define png_set_compression_strategy vtk_png_set_compression_strategy
#define png_set_compression_window_bits vtk_png_set_compression_window_bits
#define png_set_crc_action vtk_png_set_crc_action
#define png_set_dither vtk_png_set_dither
#define png_set_error_fn vtk_png_set_error_fn
#define png_set_expand vtk_png_set_expand
#define png_set_filler vtk_png_set_filler
#define png_set_filter vtk_png_set_filter
#define png_set_filter_heuristics vtk_png_set_filter_heuristics
#define png_set_flush vtk_png_set_flush
#define png_set_gAMA vtk_png_set_gAMA
#define png_set_gAMA_fixed vtk_png_set_gAMA_fixed
#define png_set_gamma vtk_png_set_gamma
#define png_set_gray_1_2_4_to_8 vtk_png_set_gray_1_2_4_to_8
#define png_set_gray_to_rgb vtk_png_set_gray_to_rgb
#define png_set_hIST vtk_png_set_hIST
#define png_set_iCCP vtk_png_set_iCCP
#define png_set_interlace_handling vtk_png_set_interlace_handling
#define png_set_invalid vtk_png_set_invalid
#define png_set_invert_alpha vtk_png_set_invert_alpha
#define png_set_invert_mono vtk_png_set_invert_mono
#define png_set_keep_unknown_chunks vtk_png_set_keep_unknown_chunks
#define png_set_oFFs vtk_png_set_oFFs
#define png_set_pCAL vtk_png_set_pCAL
#define png_set_pHYs vtk_png_set_pHYs
#define png_set_packing vtk_png_set_packing
#define png_set_packswap vtk_png_set_packswap
#define png_set_palette_to_rgb vtk_png_set_palette_to_rgb
#define png_set_progressive_read_fn vtk_png_set_progressive_read_fn
#define png_set_read_fn vtk_png_set_read_fn
#define png_set_read_status_fn vtk_png_set_read_status_fn
#define png_set_read_user_chunk_fn vtk_png_set_read_user_chunk_fn
#define png_set_read_user_transform_fn vtk_png_set_read_user_transform_fn
#define png_set_rgb_to_gray vtk_png_set_rgb_to_gray
#define png_set_rgb_to_gray_fixed vtk_png_set_rgb_to_gray_fixed
#define png_set_rows vtk_png_set_rows
#define png_set_sBIT vtk_png_set_sBIT
#define png_set_sCAL vtk_png_set_sCAL
#define png_set_sPLT vtk_png_set_sPLT
#define png_set_sRGB vtk_png_set_sRGB
#define png_set_sRGB_gAMA_and_cHRM vtk_png_set_sRGB_gAMA_and_cHRM
#define png_set_shift vtk_png_set_shift
#define png_set_sig_bytes vtk_png_set_sig_bytes
#define png_set_strip_16 vtk_png_set_strip_16
#define png_set_strip_alpha vtk_png_set_strip_alpha
#define png_set_swap vtk_png_set_swap
#define png_set_swap_alpha vtk_png_set_swap_alpha
#define png_set_tIME vtk_png_set_tIME
#define png_set_tRNS vtk_png_set_tRNS
#define png_set_tRNS_to_alpha vtk_png_set_tRNS_to_alpha
#define png_set_text vtk_png_set_text
#define png_set_unknown_chunk_location vtk_png_set_unknown_chunk_location
#define png_set_unknown_chunks vtk_png_set_unknown_chunks
#define png_set_user_transform_info vtk_png_set_user_transform_info
#define png_set_write_fn vtk_png_set_write_fn
#define png_set_write_status_fn vtk_png_set_write_status_fn
#define png_set_write_user_transform_fn vtk_png_set_write_user_transform_fn
#define png_sig_cmp vtk_png_sig_cmp
#define png_start_read_image vtk_png_start_read_image
#define png_warning vtk_png_warning
#define png_write_IDAT vtk_png_write_IDAT
#define png_write_IEND vtk_png_write_IEND
#define png_write_IHDR vtk_png_write_IHDR
#define png_write_PLTE vtk_png_write_PLTE
#define png_write_bKGD vtk_png_write_bKGD
#define png_write_cHRM vtk_png_write_cHRM
#define png_write_cHRM_fixed vtk_png_write_cHRM_fixed
#define png_write_chunk vtk_png_write_chunk
#define png_write_chunk_data vtk_png_write_chunk_data
#define png_write_chunk_end vtk_png_write_chunk_end
#define png_write_chunk_start vtk_png_write_chunk_start
#define png_write_data vtk_png_write_data
#define png_write_destroy vtk_png_write_destroy
#define png_write_end vtk_png_write_end
#define png_write_filtered_row vtk_png_write_filtered_row
#define png_write_find_filter vtk_png_write_find_filter
#define png_write_finish_row vtk_png_write_finish_row
#define png_write_flush vtk_png_write_flush
#define png_write_gAMA vtk_png_write_gAMA
#define png_write_gAMA_fixed vtk_png_write_gAMA_fixed
#define png_write_hIST vtk_png_write_hIST
#define png_write_iCCP vtk_png_write_iCCP
#define png_write_image vtk_png_write_image
#define png_write_info vtk_png_write_info
#define png_write_info_before_PLTE vtk_png_write_info_before_PLTE
#define png_write_init vtk_png_write_init
#define png_write_init_2 vtk_png_write_init_2
#define png_write_init_3 vtk_png_write_init_3
#define png_write_oFFs vtk_png_write_oFFs
#define png_write_pCAL vtk_png_write_pCAL
#define png_write_pHYs vtk_png_write_pHYs
#define png_write_png vtk_png_write_png
#define png_write_row vtk_png_write_row
#define png_write_rows vtk_png_write_rows
#define png_write_sBIT vtk_png_write_sBIT
#define png_write_sCAL vtk_png_write_sCAL
#define png_write_sPLT vtk_png_write_sPLT
#define png_write_sRGB vtk_png_write_sRGB
#define png_write_sig vtk_png_write_sig
#define png_write_start_row vtk_png_write_start_row
#define png_write_tEXt vtk_png_write_tEXt
#define png_write_tIME vtk_png_write_tIME
#define png_write_tRNS vtk_png_write_tRNS
#define png_write_zTXt vtk_png_write_zTXt
#define png_zalloc vtk_png_zalloc
#define png_zfree vtk_png_zfree

#define png_zTXt vtk_png_zTXt
#define png_tRNS vtk_png_tRNS
#define png_tIME vtk_png_tIME
#define png_tEXt vtk_png_tEXt
#define png_sRGB vtk_png_sRGB
#define png_sPLT vtk_png_sPLT
#define png_sBIT vtk_png_sBIT
#define png_pHYs vtk_png_pHYs
#define png_sCAL vtk_png_sCAL
#define png_pCAL vtk_png_pCAL
#define png_oFFs vtk_png_oFFs
#define png_iTXt vtk_png_iTXt
#define png_iCCP vtk_png_iCCP
#define png_hIST vtk_png_hIST
#define png_gAMA vtk_png_gAMA
#define png_cHRM vtk_png_cHRM
#define png_bKGD vtk_png_bKGD
#define png_PLTE vtk_png_PLTE
#define png_IEND vtk_png_IEND
#define png_IDAT vtk_png_IDAT
#define png_IHDR vtk_png_IHDR
#define png_sig  vtk_png_sig
#define png_sig_bytes vtk_png_sig_bytes

#endif
