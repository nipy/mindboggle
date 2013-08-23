# Converting Olivier's fundi into VTK format. 
# Two inputs: tex-format fundi, VTK file of original mesh
# Output: one VTK file where fundi is labeled as non -1 

# Gifti to tex conversion done by using AimsFileConvert in BrainVISA
# for i in *; do AimsFileConvert -i $i/$i\_Lwhite_sulcalines.gii -o $i/L.tex -e 1 ; done
# for i in *; do AimsFileConvert -i $i/$i\_Rwhite_sulcalines.gii -o $i/R.tex -e 1 ; done

def load_tex(File):
    """Load the 4th line of a tex file
    """
    with open(File, 'r') as f:
        for i in range(3): # skip the first 3 lines
            f.readline() 
        last_line = f.readline()

    old_fundi = map(int, last_line.split())    
    new_fundi = [x if x >0 else -1 for x in old_fundi]

    return new_fundi

def load_mesh(File):
    """Load original mesh
    """
    from mindboggle.utils.io_vtk import read_vtk
    faces, u2, u3, points, u5, u6, u7, u8 = read_vtk(File) 

    return faces, points

def merge(fundi_LUT, faces, points, output_vtk):
    from mindboggle.utils.io_vtk import write_vtk
    write_vtk(output_vtk, points, faces=faces, scalars=[fundi_LUT[2:]], scalar_names=['fundi'])

def short2long():
    """convert short folder names in Olivier's to long name in MB101
    """
    pairs = []# each tuple is (oliver name, MB101 name)
    conversion_dict = {"m":"MMRR-21-", "nt": "NKI-TRT-20-", "o":"OASIS-TRT-20-"}
    group_size = {"o":([1,2]+range(4,20+1)), "m":range(1,19+1), "nt":range(1,20+1)} # o3 is missing
    for group_short_name, num_subj in group_size.iteritems():
        pairs += [(group_short_name+(str(subj_id)), conversion_dict[group_short_name]+str(subj_id)) for subj_id in num_subj]

    return pairs

def loop_thru(pairs):
    import os.path
    import os
    for (short_name, long_name) in pairs:
        for hemisphere in ["L", "R"]:
            tex_file = os.path.join("/data/data/Mindboggle_MRI/MB101/results/Olivier_Fundi", short_name, hemisphere+".tex")
            print tex_file
            fundi_LUT = load_tex(tex_file)

            hemisphere = hemisphere.lower()
            long_path = os.path.join("/data/data/Mindboggle_MRI/MB101/results/features", "_".join(["_hemi", hemisphere+"h", "subject", long_name]))
            mesh_file = os.path.join(long_path, "folds.vtk")
            print mesh_file
            faces, points = load_mesh(mesh_file)

            output_path = os.path.join("/data/data/Mindboggle_MRI/MB101/results/Olivier_Fundi", "_".join(["_hemi", hemisphere+"h", "subject", long_name]))
            
            try:
                os.mkdir(output_path)
            except OSError:
                pass
            output_vtk = os.path.join(output_path, hemisphere+"h.pial.fundi.vtk")
            print output_vtk
#            output_vtk = "/dev/null"
            merge(fundi_LUT, faces, points, output_vtk)
            

if __name__ == "__main__":
    import sys
#    fundi_LUT = load_tex(sys.argv[1])
#    mesh = load_mesh()
    pairs = short2long()
    loop_thru(pairs)


