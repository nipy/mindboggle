Files in this folder:

1. label21.py: Traverse a path. For each subject directory in the path, load all 21 atlases-to-patient label sets (21 FreeSurfer .annot files) and vote out one label for every vertex on the hemisphere. Result saved in VTK format. 
Since the voting takes time, the result is saved as a pickle file. 
Input: atlases-to-patient manual label files. 
Output: assignment VTK files in name ?h.assign.{pial,inflated}.vtk in each subject's surf folder. 

Update 2012-02-13: Now output labels are aggregated. 
2= 2,10,23,26
3= 3, 27
18=18,19,20

2. segment.py: Assign TWO NEAREST labels to every fundus vertex. Thus, this segments fundi. Segmentation results are saved in VTK. 
Input: fundi files, neighboring files (?h.vrtx.nbr), assignment VTK files
Output: segmented fundi file with suffix seg.vtk or seg.2nd.vtk. These VTK files have additional SCALARS blocks for segmentation information. 

Update 2012-02-13: Label pairs are aggregated as follows. 
        precentral: [28,24]*,[3,24]*,[18,24]*
        postcentral: [22,29],[22,31]
        intraparietal: [29,31],[29,8]
        lateral occipital sulcus: [11,8]*, [11,29]*
        anterior occipital sulcus: [11,15]*,[11,9]
        circular sulcus: [35,30],[35,34],[35,12],[35,2],[35,24],[35,22],[35,31]
        cingulate sulcus: [2,14],[2,28],[2,17],[25,17]
        calcarine fissure: [13,25],[13,2]
        lateral H-shaped orbital sulcus: [12,18],[12,3]
        occipitotemporal sulcus: [7,9],[7,11]
        collateral sulcus: [7,6],[7,16],[7,13]
        interhemispheric fissure, dorsal margin: [17,28],[17,24],[17,22],[25,29],[5,29],[5,11]

3. shape_table.py: Generate shape table. Fundus vertexes of the same AGGREGTED label pair are grouped together.
Input: segmented fundi in VTK. 
Output: two kinds of shape tables, the average one and the individual one. In TSV. 

4. segment_Yrjo.py: Assign up to 4 nearest labels to Yrjo's pits and return distances from every pit to its 4 labels. 

Pipeline to use scripts:
1. label21.py
2. segment.py 
3. shape_table.py

=======Useful label pairs before aggregation 


Lateral surface:
frontomarginal sulcus: [28,12] 
superior frontal: [28,3] = [[28,3], [28,27]]
inferior frontal: [3,18] = [[3,18] x, [3,19] x, [3,20]x, [27,18] x, [27,19] x, [27,20]]
precentral: [28,24] x,[3,24] = [[3,24],[27,24] x], [18,24] = [[18,24], [19,24] x, [20,24] x]
central sulcus: [24,22]
postcentral: [22,29], [22,31]
intraparietal: [29,31] x,[29,8] 
primary intermediate sulcus/1st segment of the posterior superior temporal sulcus: [31,8]
sylvian fissure: [31,30] x
lateral occipital sulcus: [11,8], [11,29] x
anterior occipital sulcus: [11,15] x, [11,9] x
superior temporal sulcus: [30,15]
inferior temporal sulcus: [15,9]

PeriSylvian area (folds within the Sylvian fissure):
circular sulcus: [35,30], [35,34], [35,12], [35,2]=[[35,2] x, [35, 10] x, [35,23] x, [35,26] x], [35,24], [35,22], [35,31] x
1st transverse temporal sulcus: [30,34] x
Heschl’s sulcus: [30,34] x

Medial surface:
cingulate sulcus: [2,14]=[[2,14] x, [10,14] x, [23,14] x, [26,14] x], [2,28] =[[2,28], [10,28] x, [23,28] x, [26,28]],  [2,17] = [[2,17] x, [10,17] x, [23,17], [26,17] x], [25,17]
paracentral sulcus: [17,28] x
parietooccipital fissure: [5,25]
calcarine fissure: [13,25] x,[13,2] = [[2,13] x, [10,13] x, [23,13] x, [26,13] x]
superior rostral sulcus: [28,14] x
callosal sulcus: [2,corpus callosum]

Ventral surface:
lateral H-shaped orbital sulcus: [12,18] = [[12,18] x, [12,19], [12,20] x], [12,3] = [[12,3] x, [12,27] x]
olfactory sulcus: [12,14] x
occipitotemporal sulcus: [7,9], [7,11] x
collateral sulcus: [7,6] x, [7,16] x, [7,13]

Regions bounded by Sulcal Margins:
interhemispheric fissure, dorsal margin: [17,28] x, [17,24] x, [17,22] x, [25,29] x, [5,29] x, [5,11] x
calcarine sulcus, dorsal margin: [5,21]
calcarine sulcus, ventral margin: [21,13]
