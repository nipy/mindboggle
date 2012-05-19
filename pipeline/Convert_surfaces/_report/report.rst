Node: Convert_surfaces (utility)
================================

 Hierarchy : pipeline.Convert_surfaces
 Exec ID : Convert_surfaces

Original Inputs
---------------

* function_str : S'def convert_to_vtk(input_surface_files):\n    """Measure\n\n    measure_()\n    """\n    import subprocess as sp\n    surface_files = []\n    for input_surface_file in input_surface_files:\n        surface_file = input_surface_file + \'.vtk\'\n        surface_files.append(surface_file)\n        cmd = [\'mris_convert\', input_surface_file, surface_file]\n        proc = sp.Popen(\' \'.join(cmd))\n        o, e = proc.communicate()\n        if proc.returncode > 0 :\n            raise Exception(\'\\n\'.join([cmd + \' failed\', o, e]))\n    return surface_files\n'
.
* ignore_exception : False
* input_surface_files : <undefined>

