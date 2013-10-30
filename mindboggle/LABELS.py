#!/usr/bin/env python
"""
Brain label numbers, names, and colormap related to the DKT labeling protocol.

For more information about the Desikan-Killiany-Tourville
cortical labeling protocol see http://mindboggle.info/data
and the article:

http://www.frontiersin.org/Brain_Imaging_Methods/10.3389/fnins.2012.00171/full
"101 labeled brain images and a consistent human cortical labeling protocol"
Arno Klein, Jason Tourville. Frontiers in Brain Imaging Methods. 6:171.
DOI: 10.3389/fnins.2012.00171

==============================================================================
Combined/eliminated regions:
==============================================================================
(1) Temporal (33) and frontal (32) poles, and bankstss (see note #1)
    regions eliminated, corresponding cortex absorbed by adjacent regions.
(2) Caudal (2), isthmus (10), posterior (23), and rostral anterior (26)
    cingulate combined to form single cingulate region (2).
(3) Caudal (3) and rostral (27) middle frontal regions combined
    to form single middle frontal region (3).
(4) Opercular (18), orbital (19), and triangular (20) inferior frontal regions
    combined to form a single inferior frontal region (18).
This is a perfectly reasonable aggregation of these regions and is the one
reflected in the sulcus/region pairings below. An alternative breakdown
would be to lump 19 with lateral orbitofrontal cortex (12) and use the
anterior horizontal ramus of the sylvian fissure as the boundary between
18 and 12. Anatomically, both aggregations are defensible but one or the other
may suit your needs better.

Regarding my earlier note about the lack of a full, consistent sulcal
anterior boundary for the inferior frontal gyrus:
This will be the case for several regions, i.e., in practice, many boundaries
are not formed by sulci but instead require "jumps" across gyri
(paths along regions of different direction curvature). This can be variable,
(e.g., the precentral sulcus is consistently formed by 2 or more disconnected
components) or implicit in the definition of the boundary (e.g., the anterior
boundary between orbital inferior frontal gyrus (19) and rostral middle
frontal gyrus (27) requires a "jump" over the lateral orbital gyrus.
Below, I note with a '*' those boundaries given principally by a sulcal fundus
but which frequently require "jumps" across gyri. I handle separately
definitions that explicitly rely on non-fundus boundaries, i.e., those that
rely on the margins of sulcal banks.

(5) Parahippocampal + entorhinal cortex + and lingual gyrus?

==============================================================================
Regions bounded by sulcal fundi:
==============================================================================
Lateral surface:
------------------------------------------------------------------------------
frontomarginal sulcus: [12,28]
superior frontal: [3,28],[27,28]
inferior frontal: [3,18],[3,19],[3,20], [18,27],[19,27],[20,27]
precentral: [24,28]*, [[3,24],[24,27]]*, [[18,24],[19,24],[20,24]]*
central sulcus: [22,24]
postcentral: [22,29],[22,31], not:[22,24]
intraparietal: [29,31], [8,29]
primary intermediate sulcus /
    1st segment of the posterior superior temporal sulcus: [8,31]*
sylvian fissure: [30,31]*, not:[18,30] (see note #2)
lateral occipital sulcus: [8,11]*,[11,29]*
anterior occipital sulcus: [11,15]*,[9,11]
superior temporal sulcus: [15,30]
inferior temporal sulcus: [9,15]*
------------------------------------------------------------------------------
PeriSylvian area (folds within the Sylvian fissure):
------------------------------------------------------------------------------
circular sulcus: [12,35],[30,35],[34,35], [2,35],[10,35],[23,35],[26,35],
                 [22,35], [24,35], [31,35]
1st transverse temporal sulcus: [30,34]
Heschl's sulcus: [30,34]


    # DKT31 to DKT25: [[10,23,26,27,19,20], [2,2,2,3,18,18]]


------------------------------------------------------------------------------
Medial surface:
------------------------------------------------------------------------------
cingulate sulcus: [2,14],[10,14],[14,23],[14,26] (see note #3),
                  [2,28],[10,28],[23,28],[26,28],
                  [2,17],[10,17],[17,23],[17,26], [17,25]
paracentral sulcus: [17,28]*
parietooccipital fissure: [5,25]
calcarine fissure: [13,25], [2,13],[10,13],[13,23],[13,26] not:[5,13] (note #4)
superior rostral sulcus: [14,28]
callosal sulcus: [2,4],[4,10],[4,23],[4,26]
------------------------------------------------------------------------------
Ventral surface:
------------------------------------------------------------------------------
lateral H-shaped orbital sulcus: [3,12],[12,27], [12,18],[12,19],[12,20]
olfactory sulcus: [12,14]
occipitotemporal sulcus: [7,9],[7,11]
collateral sulcus: [6,7], [7,13], [7,16]

==============================================================================
What boundaries will NEVER be derived by fundi, but instead by curvature, etc.
==============================================================================
Regions bounded by sulcal margins:
------------------------------------------------------------------------------
interhemispheric fissure, dorsal margin:
    [17,28],[17,24],[17,22],[25,29],[5,29],[5,11]
calcarine sulcus, dorsal margin: [5,21]
calcarine sulcus, ventral margin: [21,13]
------------------------------------------------------------------------------
Regions with additional non-sulcal boundaries with subcortical regions:
------------------------------------------------------------------------------
[16,6,9,30,12,14]

==============================================================================
Notes:
==============================================================================
[1] This was eliminated b/c it spanned the superior temporal sulcus fundus
    and because the anterior boundary was ambiguous.
[2] The insula lies b/w these regions and is separated from them by the
    circular sulcus which is marked by an easily distinguished fold deep
    within the Sylvian fissure.
[3] This is the case in some, but not all, hemispheres. It occurs when the
    superior rostral sulcus fails to intersect with the cingulate sulcus.
[4] The pericalcarine region lies between these 2 regions. As defined in
    "Regions bounded by sulcal margins", the pericalcarine cortex (21)
    dorsal (with 5) and ventral (with 13) boundaries are formed by the
    lateral margins of the dorsal and ventral banks of the calcarine sulcus
    rather than a sulcal fundus; because this region spans the sulcal fundus,
    we cannot simply incorporate portions of the region into the adjacent
    regions based on the fundus.

Authors:
    - Jason Tourville, 2011-2012  (jtour@bu.edu)
    - Arno Klein, 2011-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=============================================================================
# DKT protocol
#=============================================================================
class DKTprotocol:
    """
    Variables related to the Desikan-Killiany-Tourville labeling protocol.

    Returns
    -------
    sulcus_names : list of strings
        sulcus names
    unique_sulcus_label_pairs : list of unique pairs of integers
        unique label pairs corresponding to label boundaries / sulcus / fundus
    sulcus_label_pair_lists : list of two lists of multiple lists of integer pairs
        list containing left and right lists, each with multiple lists of
        integer pairs corresponding to label boundaries / sulcus / fundus
    left_sulcus_label_pair_lists : list of two lists of lists of integer pairs
    right_sulcus_label_pair_lists : list of two lists of lists of integer pairs
    label_names : list of strings
        label names
    left_label_names : list of strings
    right_label_names : list of strings
    label_numbers : list of integers
        label numbers
    left_label_numbers : list of integers
    right_label_numbers : list of integers
    cerebrum_names : list of strings
        label names
    left_cerebrum_names : list of strings
    right_cerebrum_names : list of strings
    cerebrum_numbers : list of integers
        label numbers
    left_cerebrum_numbers : list of integers
    right_cerebrum_numbers : list of integers
    cerebrum_cortex_names : list of strings
        label names for cortical regions in the DKT31 protocol
    left_cerebrum_cortex_names : list of strings
    right_cerebrum_cortex_names : list of strings
    cerebrum_cortex_numbers : list of integers
        label numbers for cortical regions in the DKT31 protocol
    left_cerebrum_cortex_numbers : list of integers
    right_cerebrum_cortex_numbers : list of integers
    cerebrum_noncortex_names : list of strings
        label names for noncortical regions in the CMA protocol
    left_cerebrum_noncortex_names : list of strings
    right_cerebrum_noncortex_names : list of strings
    medial_cerebrum_noncortex_names : list of strings
    cerebrum_noncortex_numbers : list of integers
        label numbers for noncortical regions in the CMA protocol
    left_cerebrum_noncortex_numbers : list of integers
    right_cerebrum_noncortex_numbers : list of integers
    medial_cerebrum_noncortex_numbers : list of integers
    cerebrum_cortex_names_DKT25 : list of strings
        label names for cortical regions in the DKT25 protocol
    left_cerebrum_cortex_names_DKT25 : list of strings
    right_cerebrum_cortex_names_DKT25 : list of strings
    cerebrum_cortex_numbers_DKT25 : list of integers
        label numbers for cortical regions in the DKT25 protocol
    left_cerebrum_cortex_numbers_DKT25 : list of integers
    right_cerebrum_cortex_numbers_DKT25 : list of integers
    cerebellum_names : list of strings
        label names
    left_cerebellum_names : list of strings
    right_cerebellum_names : list of strings
    cerebellum_numbers : list of integers
        label numbers
    left_cerebellum_numbers : list of integers
    right_cerebellum_numbers : list of integers
    cerebellum_cortex_names : list of strings
        label names for cortical regions in the DKT31 protocol
    left_cerebellum_cortex_names : list of strings
    right_cerebellum_cortex_names : list of strings
    cerebellum_cortex_numbers : list of integers
        label numbers for cortical regions in the DKT31 protocol
    left_cerebellum_cortex_numbers : list of integers
    right_cerebellum_cortex_numbers : list of integers
    cerebellum_noncortex_names : list of strings
        label names for noncortical regions in the CMA protocol
    left_cerebellum_noncortex_names : list of strings
    right_cerebellum_noncortex_names : list of strings
    medial_cerebellum_noncortex_names : list of strings
    cerebellum_noncortex_numbers : list of integers
        label numbers for noncortical regions in the CMA protocol
    left_cerebellum_noncortex_numbers : list of integers
    right_cerebellum_noncortex_numbers : list of integers
    medial_cerebellum_noncortex_numbers : list of integers

    Examples
    --------
    >>> from mindboggle.LABELS import DKTprotocol
    >>> dkt = DKTprotocol()
    >>> dkt.left_cerebrum_names

    """
    #-------------------------------------------------------------------------
    # DKT cerebral cortical labeling protocol -- 31 labels:
    #-------------------------------------------------------------------------
    left_cerebrum_cortex_numbers_names = [
       #[3,       "left cerebral cortex"],
        [1002,    "left caudal anterior cingulate"],
        [1003,    "left caudal middle frontal"],
        [1005,    "left cuneus"],
        [1006,    "left entorhinal"],
        [1007,    "left fusiform"],
        [1008,    "left inferior parietal"],
        [1009,    "left inferior temporal"],
        [1010,    "left isthmus cingulate"],
        [1011,    "left lateral occipital"],
        [1012,    "left lateral orbitofrontal"],
        [1013,    "left lingual"],
        [1014,    "left medial orbitofrontal"],
        [1015,    "left middle temporal"],
        [1016,    "left parahippocampal"],
        [1017,    "left paracentral"],
        [1018,    "left pars opercularis"],
        [1019,    "left pars orbitalis"],
        [1020,    "left pars triangularis"],
        [1021,    "left pericalcarine"],
        [1022,    "left postcentral"],
        [1023,    "left posterior cingulate"],
        [1024,    "left precentral"],
        [1025,    "left precuneus"],
        [1026,    "left rostral anterior cingulate"],
        [1027,    "left rostral middle frontal"],
        [1028,    "left superior frontal"],
        [1029,    "left superior parietal"],
        [1030,    "left superior temporal"],
        [1031,    "left supramarginal"],
        [1034,    "left transverse temporal"],
        [1035,    "left insula"]]

    right_cerebrum_cortex_numbers_names = [
       #[42,      "right cerebral cortex"],
        [2002,    "right caudal anterior cingulate"],
        [2003,    "right caudal middle frontal"],
        [2005,    "right cuneus"],
        [2006,    "right entorhinal"],
        [2007,    "right fusiform"],
        [2008,    "right inferior parietal"],
        [2009,    "right inferior temporal"],
        [2010,    "right isthmus cingulate"],
        [2011,    "right lateral occipital"],
        [2012,    "right lateral orbitofrontal"],
        [2013,    "right lingual"],
        [2014,    "right medial orbitofrontal"],
        [2015,    "right middle temporal"],
        [2016,    "right parahippocampal"],
        [2017,    "right paracentral"],
        [2018,    "right pars opercularis"],
        [2019,    "right pars orbitalis"],
        [2020,    "right pars triangularis"],
        [2021,    "right pericalcarine"],
        [2022,    "right postcentral"],
        [2023,    "right posterior cingulate"],
        [2024,    "right precentral"],
        [2025,    "right precuneus"],
        [2026,    "right rostral anterior cingulate"],
        [2027,    "right rostral middle frontal"],
        [2028,    "right superior frontal"],
        [2029,    "right superior parietal"],
        [2030,    "right superior temporal"],
        [2031,    "right supramarginal"],
        [2034,    "right transverse temporal"],
        [2035,    "right insula"]]

    left_cerebrum_cortex_numbers = \
        [x[0] for x in left_cerebrum_cortex_numbers_names]
    right_cerebrum_cortex_numbers = \
        [x[0] for x in right_cerebrum_cortex_numbers_names]
    left_cerebrum_cortex_names = \
        [x[1] for x in left_cerebrum_cortex_numbers_names]
    right_cerebrum_cortex_names = \
        [x[1] for x in right_cerebrum_cortex_numbers_names]
    cerebrum_cortex_numbers = \
        left_cerebrum_cortex_numbers + right_cerebrum_cortex_numbers
    cerebrum_cortex_names = \
        left_cerebrum_cortex_names + right_cerebrum_cortex_names

    #-------------------------------------------------------------------------
    # DKT cerebral cortical labeling protocol -- 25 labels:
    #-------------------------------------------------------------------------
    # Region numbers:
    # DKT31 to DKT25: [[10,23,26,27,19,20], [2,2,2,3,18,18]]
    left_cerebrum_cortex_numbers_DKT25 = left_cerebrum_cortex_numbers[:]
    for n in [1010, 1019, 1020, 1023, 1026, 1027]:
        left_cerebrum_cortex_numbers_DKT25.remove(n)
    right_cerebrum_cortex_numbers_DKT25 = right_cerebrum_cortex_numbers[:]
    for n in [2010, 2019, 2020, 2023, 2026, 2027]:
        right_cerebrum_cortex_numbers_DKT25.remove(n)
    cerebrum_cortex_numbers_DKT25 = left_cerebrum_cortex_numbers_DKT25 + \
                                    right_cerebrum_cortex_numbers_DKT25

    # Consolidate region labels:
    left_cerebrum_cortex_names_DKT25 = []
    for ilabel, label_number in enumerate(left_cerebrum_cortex_numbers):
        if label_number == 1002:
            left_cerebrum_cortex_names_DKT25.append('left cingulate')
        elif label_number == 1003:
            left_cerebrum_cortex_names_DKT25.append('left middle frontal')
        elif label_number == 1018:
            left_cerebrum_cortex_names_DKT25.append('left inferior frontal')
        else:
            left_cerebrum_cortex_names_DKT25.\
                append(left_cerebrum_cortex_names[ilabel])
    right_cerebrum_cortex_names_DKT25 = []
    for ilabel, label_number in enumerate(right_cerebrum_cortex_numbers):
        if label_number == 2002:
            right_cerebrum_cortex_names_DKT25.append('right cingulate')
        elif label_number == 2003:
            right_cerebrum_cortex_names_DKT25.append('right middle frontal')
        elif label_number == 2018:
            right_cerebrum_cortex_names_DKT25.append('right inferior frontal')
        else:
            right_cerebrum_cortex_names_DKT25.\
                append(right_cerebrum_cortex_names[ilabel])
    cerebrum_cortex_names_DKT25 = \
        left_cerebrum_cortex_names_DKT25 + right_cerebrum_cortex_names_DKT25

    #-------------------------------------------------------------------------
    # Cerebral noncortex label numbers and names
    # These labels were converted from Neuromorphometrics BrainCOLOR subcortex
    # labels to be consistent with FreeSurferColorLUT.txt labels.
    #
    # Two labels did not have counterparts in FreeSurfer:
    #         [75, "left basal forebrain"],
    #         [76, "right basal forebrain"]
    # and were reassigned to unused numbers in FreeSurferColorLUT.txt:
    #         [91, "left basal forebrain"],
    #         [92, "right basal forebrain"]
    #-------------------------------------------------------------------------
    left_cerebrum_noncortex_numbers_names = [
        [2, "left cerebral white matter"],
        [4, "left lateral ventricle"],
        [5, "left inferior lateral ventricle"],
        [10, "left thalamus proper"],
        [11, "left caudate"],
        [12, "left putamen"],
        [13, "left pallidum"],
        [17, "left hippocampus"],
        [18, "left amygdala"],
        [25, "left lesion"],
        [26, "left accumbens area"],
        [28, "left ventral DC"],
        [30, "left vessel"],
        [91, "left basal forebrain"]]
    right_cerebrum_noncortex_numbers_names = [
        [41, "right cerebral white matter"],
        [43, "right lateral ventricle"],
        [44, "right inferior lateral ventricle"],
        [49, "right thalamus proper"],
        [50, "right caudate"],
        [51, "right putamen"],
        [52, "right pallidum"],
        [53, "right hippocampus"],
        [54, "right amygdala"],
        [57, "right lesion"],
        [58, "right accumbens area"],
        [60, "right ventral DC"],
        [62, "right vessel"],
        [92, "right basal forebrain"]]
    medial_cerebrum_noncortex_numbers_names = [
        [16, "Brain stem"],
        [24, "CSF"],
        [14, "3rd ventricle"],
        [15, "4th ventricle"],
        [72, "5th ventricle"],
        [85, "optic chiasm"]]

    left_cerebrum_noncortex_numbers = \
        [x[0] for x in left_cerebrum_noncortex_numbers_names]
    right_cerebrum_noncortex_numbers = \
        [x[0] for x in right_cerebrum_noncortex_numbers_names]
    medial_cerebrum_noncortex_numbers = \
        [x[0] for x in medial_cerebrum_noncortex_numbers_names]
    left_cerebrum_noncortex_names = \
        [x[1] for x in left_cerebrum_noncortex_numbers_names]
    right_cerebrum_noncortex_names = \
        [x[1] for x in right_cerebrum_noncortex_numbers_names]
    medial_cerebrum_noncortex_names = \
        [x[1] for x in medial_cerebrum_noncortex_numbers_names]
    cerebrum_noncortex_numbers = \
        left_cerebrum_noncortex_numbers + right_cerebrum_noncortex_numbers + \
        medial_cerebrum_noncortex_numbers
    cerebrum_noncortex_names = \
        left_cerebrum_noncortex_names + right_cerebrum_noncortex_names + \
        medial_cerebrum_noncortex_names

    left_cerebrum_numbers = \
        left_cerebrum_cortex_numbers + left_cerebrum_noncortex_numbers
    right_cerebrum_numbers = \
        right_cerebrum_cortex_numbers + right_cerebrum_noncortex_numbers
    cerebrum_numbers = \
        left_cerebrum_numbers + right_cerebrum_numbers + \
        medial_cerebrum_noncortex_numbers
    left_cerebrum_names = \
        left_cerebrum_cortex_names + left_cerebrum_noncortex_names
    right_cerebrum_names = \
        right_cerebrum_cortex_names + right_cerebrum_noncortex_names
    cerebrum_names = \
        left_cerebrum_names + right_cerebrum_names + \
        medial_cerebrum_noncortex_names

    #-------------------------------------------------------------------------
    # Cerebellar label numbers and names
    # These labels were converted from Neuromorphometrics BrainCOLOR subcortex
    # labels to be consistent with FreeSurferColorLUT.txt labels.
    #
    # Three labels did not have counterparts in FreeSurfer:
    #         [71, "cerebellar vermal lobules I-V"],
    #         [72, "cerebellar vermal lobules VI-VII"],
    #         [73, "cerebellar vermal lobules VIII-X"],
    # and were reassigned to unused numbers in FreeSurferColorLUT.txt:
    #         [630, "cerebellar vermal lobules I-V"],
    #         [631, "cerebellar vermal lobules VI-VII"],
    #         [632, "cerebellar vermal lobules VIII-X"],
    #-------------------------------------------------------------------------
    left_cerebellum_cortex_numbers_names = [
        [6, "left cerebellum exterior"]]
    right_cerebellum_cortex_numbers_names = [
        [45, "right cerebellum exterior"]]
    left_cerebellum_noncortex_numbers_names = [
        [7, "left cerebellum white matter"]]
    right_cerebellum_noncortex_numbers_names = [
        [46, "right cerebellum white matter"]]
    medial_cerebellum_noncortex_numbers_names = [
        [630, "cerebellar vermal lobules I-V"],
        [631, "cerebellar vermal lobules VI-VII"],
        [632, "cerebellar vermal lobules VIII-X"]]

    left_cerebellum_cortex_numbers = [x[0]
                        for x in left_cerebellum_cortex_numbers_names]
    right_cerebellum_cortex_numbers = [x[0]
                        for x in right_cerebellum_cortex_numbers_names]
    left_cerebellum_cortex_names = [x[1]
                        for x in left_cerebellum_cortex_numbers_names]
    right_cerebellum_cortex_names = [x[1]
                        for x in right_cerebellum_cortex_numbers_names]
    left_cerebellum_noncortex_numbers = [x[0]
                        for x in left_cerebellum_noncortex_numbers_names]
    right_cerebellum_noncortex_numbers = [x[0]
                        for x in right_cerebellum_noncortex_numbers_names]
    medial_cerebellum_noncortex_numbers = [x[0]
                        for x in medial_cerebellum_noncortex_numbers_names]
    left_cerebellum_noncortex_names = [x[1]
                        for x in left_cerebellum_noncortex_numbers_names]
    right_cerebellum_noncortex_names = [x[1]
                        for x in right_cerebellum_noncortex_numbers_names]
    medial_cerebellum_noncortex_names = [x[1]
                        for x in medial_cerebellum_noncortex_numbers_names]

    left_cerebellum_numbers = left_cerebellum_cortex_numbers + \
                                    left_cerebellum_noncortex_numbers
    right_cerebellum_numbers = right_cerebellum_cortex_numbers + \
                                     right_cerebellum_noncortex_numbers
    left_cerebellum_names = left_cerebellum_cortex_names + \
                                  left_cerebellum_noncortex_names
    right_cerebellum_names = right_cerebellum_cortex_names + \
                                   right_cerebellum_noncortex_names

    cerebellum_cortex_numbers = left_cerebellum_cortex_numbers + \
                               right_cerebellum_cortex_numbers
    cerebellum_noncortex_numbers = left_cerebellum_noncortex_numbers + \
                                  right_cerebellum_noncortex_numbers + \
                                  medial_cerebellum_noncortex_numbers
    cerebellum_cortex_names = left_cerebellum_cortex_names + \
                             right_cerebellum_cortex_names
    cerebellum_noncortex_names = left_cerebellum_noncortex_names + \
                                right_cerebellum_noncortex_names + \
                                medial_cerebellum_noncortex_names
    cerebellum_numbers = cerebellum_cortex_numbers + \
                              cerebellum_noncortex_numbers
    cerebellum_names = cerebellum_cortex_names + \
                            cerebellum_noncortex_names

    label_numbers = cerebrum_numbers + cerebellum_numbers
    label_names = cerebrum_names + cerebellum_names

    #-------------------------------------------------------------------------
    # Sulcus names from the DKT labeling protocol:
    #-------------------------------------------------------------------------
    sulcus_names = [
        "frontomarginal sulcus",
        "superior frontal sulcus",
        "inferior frontal sulcus",
        "precentral sulcus",
        "central sulcus",
        "postcentral sulcus",
        "intraparietal sulcus",
        "primary intermediate sulcus/1st segment of post. sup. temporal sulcus",
        "sylvian fissure",
        "lateral occipital sulcus",
        "anterior occipital sulcus",
        "superior temporal sulcus",
        "inferior temporal sulcus",
        "circular sulcus",
        "1st transverse temporal sulcus and Heschl's sulcus",
        "cingulate sulcus",
        "paracentral sulcus",
        "parietooccipital fissure",
        "calcarine fissure",
        "superior rostral sulcus",
        "callosal sulcus",
        "lateral H-shaped orbital sulcus",
        "olfactory sulcus",
        "occipitotemporal sulcus",
        "collateral sulcus"]

    #-------------------------------------------------------------------------
    # Lists of label pairs that define sulcus boundaries (or fundi)
    # according to the DKT labeling protocol.
    # 1000 [lh] or 2000 [rh] are added to these numbers below.
    #-------------------------------------------------------------------------
    pair_lists = [
        [[12,28]],
        [[3,28], [27,28]],
        [[3,18],[3,19],[3,20], [18,27],[19,27],[20,27]],
        [[24,28], [3,24],[24,27], [18,24],[19,24],[20,24]],
        [[22,24]],
        [[22,29], [22,31]],
        [[29,31], [8,29]],
        [[8,31]],
        [[30,31]],
        [[8,11], [11,29]],
        [[11,15], [9,11]],
        [[15,30]],
        [[9,15]],
        [[12,35], [30,35], [34,35], [2,35],[10,35],[23,35],[26,35],
            [22,35], [24,35], [31,35]],
        [[30,34]],
        [[2,14],[10,14],[14,23],[14,26], [2,28],[10,28],[23,28],[26,28],
            [2,17],[10,17],[17,23],[17,26], [17,25]],
        [[17,28]],
        [[5,25]],
        [[13,25], [2,13],[10,13],[13,23],[13,26]],
        [[14,28]],
        [[2,4], [4,10], [4,23], [4,26]],
        [[3,12],[12,27], [12,18],[12,19],[12,20]],
        [[12,14]],
        [[7,9], [7,11]],
        [[6,7], [7,16], [7,13]]]

    relabel = True

    left_sulcus_label_pair_lists = []
    right_sulcus_label_pair_lists = []
    unique_sulcus_label_pairs = []  # unique sorted label pairs
    for pair_list in pair_lists:
        left_pairs = []
        right_pairs = []
        for pair in pair_list:
            left_pair = [1000 + pair[0], 1000 + pair[1]]
            right_pair = [2000 + pair[0], 2000 + pair[1]]
            left_pairs.append(left_pair)
            right_pairs.append(right_pair)
            if relabel:
                if left_pair not in unique_sulcus_label_pairs:
                    unique_sulcus_label_pairs.append(left_pair)
                if right_pair not in unique_sulcus_label_pairs:
                    unique_sulcus_label_pairs.append(right_pair)
            else:
                if pair not in unique_sulcus_label_pairs:
                    unique_sulcus_label_pairs.append(pair)
        left_sulcus_label_pair_lists.append(left_pairs)
        right_sulcus_label_pair_lists.append(right_pairs)
    if relabel:
        sulcus_label_pair_lists = left_sulcus_label_pair_lists + \
                                  right_sulcus_label_pair_lists
    else:
        sulcus_label_pair_lists = pair_lists

    #=========================================================================
    # Mindboggle colormap
    #=========================================================================
    number_name_rgb = [
        [0,   "Unknown",    [0, 0, 0]],
        [1,   "Left-Cerebral-Exterior",    [70, 130, 180]],
        [2,   "Left-Cerebral-White-Matter",    [245, 245, 245]],
        [3,   "Left-Cerebral-Cortex",    [205, 62, 78]],
        [4,   "Left-Lateral-Ventricle",    [120, 18, 134]],
        [5,   "Left-Inf-Lat-Ventricle",    [196, 58, 250]],
        [6,   "Left-Cerebellum-Exterior",    [0, 148, 0]],
        [7,   "Left-Cerebellum-White-Matter",    [220, 248, 164]],
        [8,   "Left-Cerebellum-Cortex",    [230, 148, 34]],
        [9,   "Left-Thalamus",    [0, 118, 14]],
        [10,  "Left-Thalamus-Proper",    [0, 118, 14]],
        [11,  "Left-Caudate",    [122, 186, 220]],
        [12,  "Left-Putamen",    [236, 13, 176]],
        [13,  "Left-Pallidum",    [12, 48, 255]],
        [14,  "3rd-Ventricle",    [204, 182, 142]],
        [15,  "4th-Ventricle",    [42, 204, 164]],
        [16,  "Brain-Stem",    [119, 159, 176]],
        [17,  "Left-Hippocampus",    [220, 216, 20]],
        [18,  "Left-Amygdala",    [103, 255, 255]],
        [19,  "Left-Insula",    [80, 196, 98]],
        [20,  "Left-Operculum",    [60, 58, 210]],
        [21,  "Line-1",    [60, 58, 210]],
        [22,  "Line-2",    [60, 58, 210]],
        [23,  "Line-3",    [60, 58, 210]],
        [24,  "CSF",    [60, 60, 60]],
        [25,  "Left-Lesion",    [255, 165, 0]],
        [26,  "Left-Accumbens-area",    [255, 165, 0]],
        [27,  "Left-Substancia-Nigra",    [0, 255, 127]],
        [28,  "Left-VentralDC",    [165, 42, 42]],
        [29,  "Left-undetermined",    [135, 206, 235]],
        [30,  "Left-vessel",    [160, 32, 240]],
        [31,  "Left-choroid-plexus",    [0, 200, 200]],
        [32,  "Left-F3orb",    [100, 50, 100]],
        [33,  "Left-lOg",    [135, 50, 74]],
        [34,  "Left-aOg",    [122, 135, 50]],
        [35,  "Left-mOg",    [51, 50, 135]],
        [36,  "Left-pOg",    [74, 155, 60]],
        [37,  "Left-Stellate",    [120, 62, 43]],
        [38,  "Left-Porg",    [74, 155, 60]],
        [39,  "Left-Aorg",    [122, 135, 50]],
        [40,  "Right-Cerebral-Exterior",    [70, 130, 180]],
        [41,  "Right-Cerebral-White-Matter",    [0, 225, 0]],
        [42,  "Right-Cerebral-Cortex",    [205, 62, 78]],
        [43,  "Right-Lateral-Ventricle",    [120, 18, 134]],
        [44,  "Right-Inf-Lat-Ventricle",    [196, 58, 250]],
        [45,  "Right-Cerebellum-Exterior",    [0, 148, 0]],
        [46,  "Right-Cerebellum-White-Matter",    [220, 248, 164]],
        [47,  "Right-Cerebellum-Cortex",    [230, 148, 34]],
        [48,  "Right-Thalamus",    [0, 118, 14]],
        [49,  "Right-Thalamus-Proper",    [0, 118, 14]],
        [50,  "Right-Caudate",    [122, 186, 220]],
        [51,  "Right-Putamen",    [236, 13, 176]],
        [52,  "Right-Pallidum",    [13, 48, 255]],
        [53,  "Right-Hippocampus",    [220, 216, 20]],
        [54,  "Right-Amygdala",    [103, 255, 255]],
        [55,  "Right-Insula",    [80, 196, 98]],
        [56,  "Right-Operculum",    [60, 58, 210]],
        [57,  "Right-Lesion",    [255, 165, 0]],
        [58,  "Right-Accumbens-area",    [255, 165, 0]],
        [59,  "Right-Substancia-Nigra",    [0, 255, 127]],
        [60,  "Right-VentralDC",    [165, 42, 42]],
        [61,  "Right-undetermined",    [135, 206, 235]],
        [62,  "Right-vessel",    [160, 32, 240]],
        [63,  "Right-choroid-plexus",    [0, 200, 221]],
        [64,  "Right-F3orb",    [100, 50, 100]],
        [65,  "Right-lOg",    [135, 50, 74]],
        [66,  "Right-aOg",    [122, 135, 50]],
        [67,  "Right-mOg",    [51, 50, 135]],
        [68,  "Right-pOg",    [74, 155, 60]],
        [69,  "Right-Stellate",    [120, 62, 43]],
        [70,  "Right-Porg",    [74, 155, 60]],
        [71,  "Right-Aorg",    [122, 135, 50]],
        [72,  "5th-Ventricle",    [120, 190, 150]],
        [73,  "Left-Interior",    [122, 135, 50]],
        [74,  "Right-Interior",    [122, 135, 50]],
        [75,  "Left-Lateral-Ventricles",    [120, 18, 134]],
        [76,  "Right-Lateral-Ventricles",    [120, 18, 134]],
        [77,  "WM-hypointensities",    [200, 70, 255]],
        [78,  "Left-WM-hypointensities",    [255, 148, 10]],
        [79,  "Right-WM-hypointensities",    [255, 148, 10]],
        [80,  "non-WM-hypointensities",    [164, 108, 226]],
        [81,  "Left-non-WM-hypointensities",    [164, 108, 226]],
        [82,  "Right-non-WM-hypointensities",    [164, 108, 226]],
        [83,  "Left-F1",    [255, 218, 185]],
        [84,  "Right-F1",    [255, 218, 185]],
        [85,  "Optic-Chiasm",    [234, 169, 30]],
        [86,  "Corpus_Callosum",    [250, 255, 50]],
        [91,  "Left-basal-forebrain",    [250, 255, 100]],
        [92,  "Right-basal-forebrain",    [250, 255, 120]],
        [96,  "Left-Amygdala-Anterior",    [205, 10, 125]],
        [97,  "Right-Amygdala-Anterior",    [205, 10, 125]],
        [98,  "Dura",    [160, 32, 240]],
        [100, "Left-wm-intensity-abnormality",    [124, 140, 178]],
        [101, "Left-caudate-intensity-abnormality",    [125, 140, 178]],
        [102, "Left-putamen-intensity-abnormality",    [126, 140, 178]],
        [103, "Left-accumbens-intensity-abnormality",    [127, 140, 178]],
        [104, "Left-pallidum-intensity-abnormality",    [124, 141, 178]],
        [105, "Left-amygdala-intensity-abnormality",    [124, 142, 178]],
        [106, "Left-hippocampus-intensity-abnormality",    [124, 143, 178]],
        [107, "Left-thalamus-intensity-abnormality",    [124, 144, 178]],
        [108, "Left-VDC-intensity-abnormality",    [124, 140, 179]],
        [109, "Right-wm-intensity-abnormality",    [124, 140, 178]],
        [110, "Right-caudate-intensity-abnormality",    [125, 140, 178]],
        [111, "Right-putamen-intensity-abnormality",    [126, 140, 178]],
        [112, "Right-accumbens-intensity-abnormality",    [127, 140, 178]],
        [113, "Right-pallidum-intensity-abnormality",    [124, 141, 178]],
        [114, "Right-amygdala-intensity-abnormality",    [124, 142, 178]],
        [115, "Right-hippocampus-intensity-abnormality",    [124, 143, 178]],
        [116, "Right-thalamus-intensity-abnormality",    [124, 144, 178]],
        [117, "Right-VDC-intensity-abnormality",    [124, 140, 179]],
        [118, "Epidermis",    [255, 20, 147]],
        [119, "Conn-Tissue",    [205, 179, 139]],
        [120, "SC-Fat/Muscle",    [238, 238, 209]],
        [121, "Cranium",    [200, 200, 200]],
        [122, "CSF-SA",    [74, 255, 74]],
        [123, "Muscle",    [238, 0, 0]],
        [124, "Ear",    [0, 0, 139]],
        [125, "Adipose",    [173, 255, 47]],
        [126, "Spinal-Cord",    [133, 203, 229]],
        [127, "Soft-Tissue",    [26, 237, 57]],
        [128, "Nerve",    [34, 139, 34]],
        [129, "Bone",    [30, 144, 255]],
        [130, "Air",    [147, 19, 173]],
        [131, "Orbital-Fat",    [238, 59, 59]],
        [132, "Tongue",    [221, 39, 200]],
        [133, "Nasal-Structures",    [238, 174, 238]],
        [134, "Globe",    [255, 0, 0]],
        [135, "Teeth",    [72, 61, 139]],
        [136, "Left-Caudate/Putamen",    [21, 39, 132]],
        [137, "Right-Caudate/Putamen",    [21, 39, 132]],
        [138, "Left-Claustrum",    [65, 135, 20]],
        [139, "Right-Claustrum",    [65, 135, 20]],
        [140, "Cornea",    [134, 4, 160]],
        [142, "Diploe",    [221, 226, 68]],
        [143, "Vitreous-Humor",    [255, 255, 254]],
        [144, "Lens",    [52, 209, 226]],
        [145, "Aqueous-Humor",    [239, 160, 223]],
        [146, "Outer-Table",    [70, 130, 180]],
        [147, "Inner-Table",    [70, 130, 181]],
        [148, "Periosteum",    [139, 121, 94]],
        [149, "Endosteum",    [224, 224, 224]],
        [150, "R/C/S",    [255, 0, 0]],
        [151, "Iris",    [205, 205, 0]],
        [152, "SC-Adipose/Muscle",    [238, 238, 209]],
        [153, "SC-Tissue",    [139, 121, 94]],
        [154, "Orbital-Adipose",    [238, 59, 59]],
        [155, "Left-IntCapsule-Ant",    [238, 59, 59]],
        [156, "Right-IntCapsule-Ant",    [238, 59, 59]],
        [157, "Left-IntCapsule-Pos",    [62, 10, 205]],
        [158, "Right-IntCapsule-Pos",    [62, 10, 205]],
        [193, "Left-hippocampal_fissure",    [0, 196, 255]],
        [194, "Left-CADG-head",    [255, 164, 164]],
        [195, "Left-subiculum",    [196, 196, 0]],
        [196, "Left-fimbria",    [0, 100, 255]],
        [197, "Right-hippocampal_fissure",    [128, 196, 164]],
        [198, "Right-CADG-head",    [0, 126, 75]],
        [199, "Right-subiculum",    [128, 96, 64]],
        [200, "Right-fimbria",    [0, 50, 128]],
        [201, "alveus",    [255, 204, 153]],
        [202, "perforant_pathway",    [255, 128, 128]],
        [203, "parasubiculum",    [255, 255, 0]],
        [204, "presubiculum",    [64, 0, 64]],
        [205, "subiculum",    [0, 0, 255]],
        [206, "CA1",    [255, 0, 0]],
        [207, "CA2",    [128, 128, 255]],
        [208, "CA3",    [0, 128, 0]],
        [209, "CA4",    [196, 160, 128]],
        [210, "GC-DG",    [32, 200, 255]],
        [211, "HATA",    [128, 255, 128]],
        [212, "fimbria",    [204, 153, 204]],
        [213, "Lateral_ventricle",    [121, 17, 136]],
        [214, "molecular_layer_HP",    [128, 0, 0]],
        [215, "hippocampal_fissure",    [128, 32, 255]],
        [216, "entorhinal_cortex",    [255, 204, 102]],
        [217, "molecular_layer_subiculum",    [128, 128, 128]],
        [218, "Amygdala",    [104, 255, 255]],
        [219, "cerebrum_White_Matter",    [0, 226, 0]],
        [220, "cerebrum_Cortex",    [205, 63, 78]],
        [221, "Inf_Lat_Vent",    [197, 58, 250]],
        [222, "Perirhinal",    [33, 150, 250]],
        [223, "cerebrum_White_Matter_Edge",    [226, 0, 0]],
        [224, "Background",    [100, 100, 100]],
        [225, "Ectorhinal",    [197, 150, 250]],
        [250, "Fornix",    [255, 0, 0]],
        [251, "CC_Posterior",    [0, 0, 64]],
        [252, "CC_Mid_Posterior",    [0, 0, 112]],
        [253, "CC_Central",    [0, 0, 160]],
        [254, "CC_Mid_Anterior",    [0, 0, 208]],
        [255, "CC_Anterior",    [0, 0, 255]],
        [630, "Cerebellar-vermal-lobulesI-V",    [250, 255, 130]],
        [631, "Cerebellar-vermal-lobulesVI-VII",    [250, 255, 140]],
        [632, "Cerebellar-vermal-lobulesVIII-X",    [250, 255, 150]],
        [1000,    "ctx-lh-unknown",    [25, 5, 25]],
        [1001,    "ctx-lh-bankssts",    [25, 100, 40]],
        [1002,    "ctx-lh-caudalanteriorcingulate",    [125, 100, 160]],
        [1003,    "ctx-lh-caudalmiddlefrontal",    [100, 25, 0]],
        [1004,    "ctx-lh-corpuscallosum",    [120, 70, 50]],
        [1005,    "ctx-lh-cuneus",    [220, 20, 100]],
        [1006,    "ctx-lh-entorhinal",    [220, 20, 10]],
        [1007,    "ctx-lh-fusiform",    [180, 220, 140]],
        [1008,    "ctx-lh-inferiorparietal",    [220, 60, 220]],
        [1009,    "ctx-lh-inferiortemporal",    [180, 40, 120]],
        [1010,    "ctx-lh-isthmuscingulate",    [140, 20, 140]],
        [1011,    "ctx-lh-lateraloccipital",    [20, 30, 140]],
        [1012,    "ctx-lh-lateralorbitofrontal",    [35, 75, 50]],
        [1013,    "ctx-lh-lingual",    [225, 140, 140]],
        [1014,    "ctx-lh-medialorbitofrontal",    [200, 35, 75]],
        [1015,    "ctx-lh-middletemporal",    [160, 100, 50]],
        [1016,    "ctx-lh-parahippocampal",    [20, 220, 60]],
        [1017,    "ctx-lh-paracentral",    [60, 220, 60]],
        [1018,    "ctx-lh-parsopercularis",    [220, 180, 140]],
        [1019,    "ctx-lh-parsorbitalis",    [20, 100, 50]],
        [1020,    "ctx-lh-parstriangularis",    [220, 60, 20]],
        [1021,    "ctx-lh-pericalcarine",    [120, 100, 60]],
        [1022,    "ctx-lh-postcentral",    [220, 20, 20]],
        [1023,    "ctx-lh-posteriorcingulate",    [220, 180, 220]],
        [1024,    "ctx-lh-precentral",    [60, 20, 220]],
        [1025,    "ctx-lh-precuneus",    [160, 140, 180]],
        [1026,    "ctx-lh-rostralanteriorcingulate",    [80, 20, 140]],
        [1027,    "ctx-lh-rostralmiddlefrontal",    [75, 50, 125]],
        [1028,    "ctx-lh-superiorfrontal",    [20, 220, 160]],
        [1029,    "ctx-lh-superiorparietal",    [20, 180, 140]],
        [1030,    "ctx-lh-superiortemporal",    [140, 220, 220]],
        [1031,    "ctx-lh-supramarginal",    [80, 160, 20]],
        [1032,    "ctx-lh-frontalpole",    [100, 0, 100]],
        [1033,    "ctx-lh-temporalpole",    [70, 70, 70]],
        [1034,    "ctx-lh-transversetemporal",    [150, 150, 200]],
        [1035,    "ctx-lh-insula",    [255, 192, 32]],
        [2000,    "ctx-rh-unknown",    [25, 5, 25]],
        [2001,    "ctx-rh-bankssts",    [25, 100, 40]],
        [2002,    "ctx-rh-caudalanteriorcingulate",    [125, 100, 160]],
        [2003,    "ctx-rh-caudalmiddlefrontal",    [100, 25, 0]],
        [2004,    "ctx-rh-corpuscallosum",    [120, 70, 50]],
        [2005,    "ctx-rh-cuneus",    [220, 20, 100]],
        [2006,    "ctx-rh-entorhinal",    [220, 20, 10]],
        [2007,    "ctx-rh-fusiform",    [180, 220, 140]],
        [2008,    "ctx-rh-inferiorparietal",    [220, 60, 220]],
        [2009,    "ctx-rh-inferiortemporal",    [180, 40, 120]],
        [2010,    "ctx-rh-isthmuscingulate",    [140, 20, 140]],
        [2011,    "ctx-rh-lateraloccipital",    [20, 30, 140]],
        [2012,    "ctx-rh-lateralorbitofrontal",    [35, 75, 50]],
        [2013,    "ctx-rh-lingual",    [225, 140, 140]],
        [2014,    "ctx-rh-medialorbitofrontal",    [200, 35, 75]],
        [2015,    "ctx-rh-middletemporal",    [160, 100, 50]],
        [2016,    "ctx-rh-parahippocampal",    [20, 220, 60]],
        [2017,    "ctx-rh-paracentral",    [60, 220, 60]],
        [2018,    "ctx-rh-parsopercularis",    [220, 180, 140]],
        [2019,    "ctx-rh-parsorbitalis",    [20, 100, 50]],
        [2020,    "ctx-rh-parstriangularis",    [220, 60, 20]],
        [2021,    "ctx-rh-pericalcarine",    [120, 100, 60]],
        [2022,    "ctx-rh-postcentral",    [220, 20, 20]],
        [2023,    "ctx-rh-posteriorcingulate",    [220, 180, 220]],
        [2024,    "ctx-rh-precentral",    [60, 20, 220]],
        [2025,    "ctx-rh-precuneus",    [160, 140, 180]],
        [2026,    "ctx-rh-rostralanteriorcingulate",    [80, 20, 140]],
        [2027,    "ctx-rh-rostralmiddlefrontal",    [75, 50, 125]],
        [2028,    "ctx-rh-superiorfrontal",    [20, 220, 160]],
        [2029,    "ctx-rh-superiorparietal",    [20, 180, 140]],
        [2030,    "ctx-rh-superiortemporal",    [140, 220, 220]],
        [2031,    "ctx-rh-supramarginal",    [80, 160, 20]],
        [2032,    "ctx-rh-frontalpole",    [100, 0, 100]],
        [2033,    "ctx-rh-temporalpole",    [70, 70, 70]],
        [2034,    "ctx-rh-transversetemporal",    [150, 150, 200]],
        [2035,    "ctx-rh-insula",    [255, 192, 32]]]

    colormap = []
    colormap_normalized = []
    for x in number_name_rgb:
        label = x[0]
        r = x[2][0]
        g = x[2][1]
        b = x[2][2]
        r1 = r / 255.0
        g1 = g / 255.0
        b1 = b / 255.0
        colormap.append([label, 1, r, g, b])
        colormap_normalized.append([label, 1, r1, g1, b1])

def print_colormap(colormap):
    """
    Print colormap.

    Parameters
    ----------
    colormap : list of lists of string and floats
        label, 1, red, green, blue

    Examples
    --------
    >>> from mindboggle.LABELS import DKTprotocol, print_colormap
    >>> dkt = DKTprotocol()
    >>> colormap = dkt.colormap_normalized
    >>> print_colormap(colormap)

    """

    print('<ColorMap name="Mindboggle Colormap" space="RGB" indexedLookup="false">')
    print('  <NaN r="0" g="0" b="0"/>')
    print('  <Point x="-1" o="0"  r="0" g="0" b="0"/>')
    for row in colormap:
        print('  <Point x="{0}" o="1"  '
              'r="{1:1.2f}" g="{2:1.2f}" b="{3:1.2f}"/>'.
              format(row[0], row[2], row[3], row[4]))
    print('</ColorMap>')

"""
#-------------------------------------------------------------------------
# FreeSurferColorLUT label numbers, names, rgb values
# (after FreeSurferColorLUT 2012)
#-------------------------------------------------------------------------
FreeSurferColorLUT = [
    [0,   "Unknown",    [0, 0, 0]],
    [1,   "Left-Cerebral-Exterior",    [70, 130, 180]],
    [2,   "Left-Cerebral-White-Matter",    [245, 245, 245]],
    [3,   "Left-Cerebral-Cortex",    [205, 62, 78]],
    [4,   "Left-Lateral-Ventricle",    [120, 18, 134]],
    [5,   "Left-Inf-Lat-Ventricle",    [196, 58, 250]],
    [6,   "Left-Cerebellum-Exterior",    [0, 148, 0]],
    [7,   "Left-Cerebellum-White-Matter",    [220, 248, 164]],
    [8,   "Left-Cerebellum-Cortex",    [230, 148, 34]],
    [9,   "Left-Thalamus",    [0, 118, 14]],
    [10,  "Left-Thalamus-Proper",    [0, 118, 14]],
    [11,  "Left-Caudate",    [122, 186, 220]],
    [12,  "Left-Putamen",    [236, 13, 176]],
    [13,  "Left-Pallidum",    [12, 48, 255]],
    [14,  "3rd-Ventricle",    [204, 182, 142]],
    [15,  "4th-Ventricle",    [42, 204, 164]],
    [16,  "Brain-Stem",    [119, 159, 176]],
    [17,  "Left-Hippocampus",    [220, 216, 20]],
    [18,  "Left-Amygdala",    [103, 255, 255]],
    [19,  "Left-Insula",    [80, 196, 98]],
    [20,  "Left-Operculum",    [60, 58, 210]],
    [21,  "Line-1",    [60, 58, 210]],
    [22,  "Line-2",    [60, 58, 210]],
    [23,  "Line-3",    [60, 58, 210]],
    [24,  "CSF",    [60, 60, 60]],
    [25,  "Left-Lesion",    [255, 165, 0]],
    [26,  "Left-Accumbens-area",    [255, 165, 0]],
    [27,  "Left-Substancia-Nigra",    [0, 255, 127]],
    [28,  "Left-VentralDC",    [165, 42, 42]],
    [29,  "Left-undetermined",    [135, 206, 235]],
    [30,  "Left-vessel",    [160, 32, 240]],
    [31,  "Left-choroid-plexus",    [0, 200, 200]],
    [32,  "Left-F3orb",    [100, 50, 100]],
    [33,  "Left-lOg",    [135, 50, 74]],
    [34,  "Left-aOg",    [122, 135, 50]],
    [35,  "Left-mOg",    [51, 50, 135]],
    [36,  "Left-pOg",    [74, 155, 60]],
    [37,  "Left-Stellate",    [120, 62, 43]],
    [38,  "Left-Porg",    [74, 155, 60]],
    [39,  "Left-Aorg",    [122, 135, 50]],
    [40,  "Right-Cerebral-Exterior",    [70, 130, 180]],
    [41,  "Right-Cerebral-White-Matter",    [0, 225, 0]],
    [42,  "Right-Cerebral-Cortex",    [205, 62, 78]],
    [43,  "Right-Lateral-Ventricle",    [120, 18, 134]],
    [44,  "Right-Inf-Lat-Ventricle",    [196, 58, 250]],
    [45,  "Right-Cerebellum-Exterior",    [0, 148, 0]],
    [46,  "Right-Cerebellum-White-Matter",    [220, 248, 164]],
    [47,  "Right-Cerebellum-Cortex",    [230, 148, 34]],
    [48,  "Right-Thalamus",    [0, 118, 14]],
    [49,  "Right-Thalamus-Proper",    [0, 118, 14]],
    [50,  "Right-Caudate",    [122, 186, 220]],
    [51,  "Right-Putamen",    [236, 13, 176]],
    [52,  "Right-Pallidum",    [13, 48, 255]],
    [53,  "Right-Hippocampus",    [220, 216, 20]],
    [54,  "Right-Amygdala",    [103, 255, 255]],
    [55,  "Right-Insula",    [80, 196, 98]],
    [56,  "Right-Operculum",    [60, 58, 210]],
    [57,  "Right-Lesion",    [255, 165, 0]],
    [58,  "Right-Accumbens-area",    [255, 165, 0]],
    [59,  "Right-Substancia-Nigra",    [0, 255, 127]],
    [60,  "Right-VentralDC",    [165, 42, 42]],
    [61,  "Right-undetermined",    [135, 206, 235]],
    [62,  "Right-vessel",    [160, 32, 240]],
    [63,  "Right-choroid-plexus",    [0, 200, 221]],
    [64,  "Right-F3orb",    [100, 50, 100]],
    [65,  "Right-lOg",    [135, 50, 74]],
    [66,  "Right-aOg",    [122, 135, 50]],
    [67,  "Right-mOg",    [51, 50, 135]],
    [68,  "Right-pOg",    [74, 155, 60]],
    [69,  "Right-Stellate",    [120, 62, 43]],
    [70,  "Right-Porg",    [74, 155, 60]],
    [71,  "Right-Aorg",    [122, 135, 50]],
    [72,  "5th-Ventricle",    [120, 190, 150]],
    [73,  "Left-Interior",    [122, 135, 50]],
    [74,  "Right-Interior",    [122, 135, 50]],
    [75,  "Left-Lateral-Ventricles",    [120, 18, 134]],
    [76,  "Right-Lateral-Ventricles",    [120, 18, 134]],
    [77,  "WM-hypointensities",    [200, 70, 255]],
    [78,  "Left-WM-hypointensities",    [255, 148, 10]],
    [79,  "Right-WM-hypointensities",    [255, 148, 10]],
    [80,  "non-WM-hypointensities",    [164, 108, 226]],
    [81,  "Left-non-WM-hypointensities",    [164, 108, 226]],
    [82,  "Right-non-WM-hypointensities",    [164, 108, 226]],
    [83,  "Left-F1",    [255, 218, 185]],
    [84,  "Right-F1",    [255, 218, 185]],
    [85,  "Optic-Chiasm",    [234, 169, 30]],
    [86,  "Corpus_Callosum",    [250, 255, 50]],
    [96,  "Left-Amygdala-Anterior",    [205, 10, 125]],
    [97,  "Right-Amygdala-Anterior",    [205, 10, 125]],
    [98,  "Dura",    [160, 32, 240]],
    [100, "Left-wm-intensity-abnormality",    [124, 140, 178]],
    [101, "Left-caudate-intensity-abnormality",    [125, 140, 178]],
    [102, "Left-putamen-intensity-abnormality",    [126, 140, 178]],
    [103, "Left-accumbens-intensity-abnormality",    [127, 140, 178]],
    [104, "Left-pallidum-intensity-abnormality",    [124, 141, 178]],
    [105, "Left-amygdala-intensity-abnormality",    [124, 142, 178]],
    [106, "Left-hippocampus-intensity-abnormality",    [124, 143, 178]],
    [107, "Left-thalamus-intensity-abnormality",    [124, 144, 178]],
    [108, "Left-VDC-intensity-abnormality",    [124, 140, 179]],
    [109, "Right-wm-intensity-abnormality",    [124, 140, 178]],
    [110, "Right-caudate-intensity-abnormality",    [125, 140, 178]],
    [111, "Right-putamen-intensity-abnormality",    [126, 140, 178]],
    [112, "Right-accumbens-intensity-abnormality",    [127, 140, 178]],
    [113, "Right-pallidum-intensity-abnormality",    [124, 141, 178]],
    [114, "Right-amygdala-intensity-abnormality",    [124, 142, 178]],
    [115, "Right-hippocampus-intensity-abnormality",    [124, 143, 178]],
    [116, "Right-thalamus-intensity-abnormality",    [124, 144, 178]],
    [117, "Right-VDC-intensity-abnormality",    [124, 140, 179]],
    [118, "Epidermis",    [255, 20, 147]],
    [119, "Conn-Tissue",    [205, 179, 139]],
    [120, "SC-Fat/Muscle",    [238, 238, 209]],
    [121, "Cranium",    [200, 200, 200]],
    [122, "CSF-SA",    [74, 255, 74]],
    [123, "Muscle",    [238, 0, 0]],
    [124, "Ear",    [0, 0, 139]],
    [125, "Adipose",    [173, 255, 47]],
    [126, "Spinal-Cord",    [133, 203, 229]],
    [127, "Soft-Tissue",    [26, 237, 57]],
    [128, "Nerve",    [34, 139, 34]],
    [129, "Bone",    [30, 144, 255]],
    [130, "Air",    [147, 19, 173]],
    [131, "Orbital-Fat",    [238, 59, 59]],
    [132, "Tongue",    [221, 39, 200]],
    [133, "Nasal-Structures",    [238, 174, 238]],
    [134, "Globe",    [255, 0, 0]],
    [135, "Teeth",    [72, 61, 139]],
    [136, "Left-Caudate/Putamen",    [21, 39, 132]],
    [137, "Right-Caudate/Putamen",    [21, 39, 132]],
    [138, "Left-Claustrum",    [65, 135, 20]],
    [139, "Right-Claustrum",    [65, 135, 20]],
    [140, "Cornea",    [134, 4, 160]],
    [142, "Diploe",    [221, 226, 68]],
    [143, "Vitreous-Humor",    [255, 255, 254]],
    [144, "Lens",    [52, 209, 226]],
    [145, "Aqueous-Humor",    [239, 160, 223]],
    [146, "Outer-Table",    [70, 130, 180]],
    [147, "Inner-Table",    [70, 130, 181]],
    [148, "Periosteum",    [139, 121, 94]],
    [149, "Endosteum",    [224, 224, 224]],
    [150, "R/C/S",    [255, 0, 0]],
    [151, "Iris",    [205, 205, 0]],
    [152, "SC-Adipose/Muscle",    [238, 238, 209]],
    [153, "SC-Tissue",    [139, 121, 94]],
    [154, "Orbital-Adipose",    [238, 59, 59]],
    [155, "Left-IntCapsule-Ant",    [238, 59, 59]],
    [156, "Right-IntCapsule-Ant",    [238, 59, 59]],
    [157, "Left-IntCapsule-Pos",    [62, 10, 205]],
    [158, "Right-IntCapsule-Pos",    [62, 10, 205]],
    [193, "Left-hippocampal_fissure",    [0, 196, 255]],
    [194, "Left-CADG-head",    [255, 164, 164]],
    [195, "Left-subiculum",    [196, 196, 0]],
    [196, "Left-fimbria",    [0, 100, 255]],
    [197, "Right-hippocampal_fissure",    [128, 196, 164]],
    [198, "Right-CADG-head",    [0, 126, 75]],
    [199, "Right-subiculum",    [128, 96, 64]],
    [200, "Right-fimbria",    [0, 50, 128]],
    [201, "alveus",    [255, 204, 153]],
    [202, "perforant_pathway",    [255, 128, 128]],
    [203, "parasubiculum",    [255, 255, 0]],
    [204, "presubiculum",    [64, 0, 64]],
    [205, "subiculum",    [0, 0, 255]],
    [206, "CA1",    [255, 0, 0]],
    [207, "CA2",    [128, 128, 255]],
    [208, "CA3",    [0, 128, 0]],
    [209, "CA4",    [196, 160, 128]],
    [210, "GC-DG",    [32, 200, 255]],
    [211, "HATA",    [128, 255, 128]],
    [212, "fimbria",    [204, 153, 204]],
    [213, "Lateral_ventricle",    [121, 17, 136]],
    [214, "molecular_layer_HP",    [128, 0, 0]],
    [215, "hippocampal_fissure",    [128, 32, 255]],
    [216, "entorhinal_cortex",    [255, 204, 102]],
    [217, "molecular_layer_subiculum",    [128, 128, 128]],
    [218, "Amygdala",    [104, 255, 255]],
    [219, "cerebrum_White_Matter",    [0, 226, 0]],
    [220, "cerebrum_Cortex",    [205, 63, 78]],
    [221, "Inf_Lat_Vent",    [197, 58, 250]],
    [222, "Perirhinal",    [33, 150, 250]],
    [223, "cerebrum_White_Matter_Edge",    [226, 0, 0]],
    [224, "Background",    [100, 100, 100]],
    [225, "Ectorhinal",    [197, 150, 250]],
    [250, "Fornix",    [255, 0, 0]],
    [251, "CC_Posterior",    [0, 0, 64]],
    [252, "CC_Mid_Posterior",    [0, 0, 112]],
    [253, "CC_Central",    [0, 0, 160]],
    [254, "CC_Mid_Anterior",    [0, 0, 208]],
    [255, "CC_Anterior",    [0, 0, 255]],
    [1000,    "ctx-lh-unknown",    [25, 5, 25]],
    [1001,    "ctx-lh-bankssts",    [25, 100, 40]],
    [1002,    "ctx-lh-caudalanteriorcingulate",    [125, 100, 160]],
    [1003,    "ctx-lh-caudalmiddlefrontal",    [100, 25, 0]],
    [1004,    "ctx-lh-corpuscallosum",    [120, 70, 50]],
    [1005,    "ctx-lh-cuneus",    [220, 20, 100]],
    [1006,    "ctx-lh-entorhinal",    [220, 20, 10]],
    [1007,    "ctx-lh-fusiform",    [180, 220, 140]],
    [1008,    "ctx-lh-inferiorparietal",    [220, 60, 220]],
    [1009,    "ctx-lh-inferiortemporal",    [180, 40, 120]],
    [1010,    "ctx-lh-isthmuscingulate",    [140, 20, 140]],
    [1011,    "ctx-lh-lateraloccipital",    [20, 30, 140]],
    [1012,    "ctx-lh-lateralorbitofrontal",    [35, 75, 50]],
    [1013,    "ctx-lh-lingual",    [225, 140, 140]],
    [1014,    "ctx-lh-medialorbitofrontal",    [200, 35, 75]],
    [1015,    "ctx-lh-middletemporal",    [160, 100, 50]],
    [1016,    "ctx-lh-parahippocampal",    [20, 220, 60]],
    [1017,    "ctx-lh-paracentral",    [60, 220, 60]],
    [1018,    "ctx-lh-parsopercularis",    [220, 180, 140]],
    [1019,    "ctx-lh-parsorbitalis",    [20, 100, 50]],
    [1020,    "ctx-lh-parstriangularis",    [220, 60, 20]],
    [1021,    "ctx-lh-pericalcarine",    [120, 100, 60]],
    [1022,    "ctx-lh-postcentral",    [220, 20, 20]],
    [1023,    "ctx-lh-posteriorcingulate",    [220, 180, 220]],
    [1024,    "ctx-lh-precentral",    [60, 20, 220]],
    [1025,    "ctx-lh-precuneus",    [160, 140, 180]],
    [1026,    "ctx-lh-rostralanteriorcingulate",    [80, 20, 140]],
    [1027,    "ctx-lh-rostralmiddlefrontal",    [75, 50, 125]],
    [1028,    "ctx-lh-superiorfrontal",    [20, 220, 160]],
    [1029,    "ctx-lh-superiorparietal",    [20, 180, 140]],
    [1030,    "ctx-lh-superiortemporal",    [140, 220, 220]],
    [1031,    "ctx-lh-supramarginal",    [80, 160, 20]],
    [1032,    "ctx-lh-frontalpole",    [100, 0, 100]],
    [1033,    "ctx-lh-temporalpole",    [70, 70, 70]],
    [1034,    "ctx-lh-transversetemporal",    [150, 150, 200]],
    [1035,    "ctx-lh-insula",    [255, 192, 32]],
    [2000,    "ctx-rh-unknown",    [25, 5, 25]],
    [2001,    "ctx-rh-bankssts",    [25, 100, 40]],
    [2002,    "ctx-rh-caudalanteriorcingulate",    [125, 100, 160]],
    [2003,    "ctx-rh-caudalmiddlefrontal",    [100, 25, 0]],
    [2004,    "ctx-rh-corpuscallosum",    [120, 70, 50]],
    [2005,    "ctx-rh-cuneus",    [220, 20, 100]],
    [2006,    "ctx-rh-entorhinal",    [220, 20, 10]],
    [2007,    "ctx-rh-fusiform",    [180, 220, 140]],
    [2008,    "ctx-rh-inferiorparietal",    [220, 60, 220]],
    [2009,    "ctx-rh-inferiortemporal",    [180, 40, 120]],
    [2010,    "ctx-rh-isthmuscingulate",    [140, 20, 140]],
    [2011,    "ctx-rh-lateraloccipital",    [20, 30, 140]],
    [2012,    "ctx-rh-lateralorbitofrontal",    [35, 75, 50]],
    [2013,    "ctx-rh-lingual",    [225, 140, 140]],
    [2014,    "ctx-rh-medialorbitofrontal",    [200, 35, 75]],
    [2015,    "ctx-rh-middletemporal",    [160, 100, 50]],
    [2016,    "ctx-rh-parahippocampal",    [20, 220, 60]],
    [2017,    "ctx-rh-paracentral",    [60, 220, 60]],
    [2018,    "ctx-rh-parsopercularis",    [220, 180, 140]],
    [2019,    "ctx-rh-parsorbitalis",    [20, 100, 50]],
    [2020,    "ctx-rh-parstriangularis",    [220, 60, 20]],
    [2021,    "ctx-rh-pericalcarine",    [120, 100, 60]],
    [2022,    "ctx-rh-postcentral",    [220, 20, 20]],
    [2023,    "ctx-rh-posteriorcingulate",    [220, 180, 220]],
    [2024,    "ctx-rh-precentral",    [60, 20, 220]],
    [2025,    "ctx-rh-precuneus",    [160, 140, 180]],
    [2026,    "ctx-rh-rostralanteriorcingulate",    [80, 20, 140]],
    [2027,    "ctx-rh-rostralmiddlefrontal",    [75, 50, 125]],
    [2028,    "ctx-rh-superiorfrontal",    [20, 220, 160]],
    [2029,    "ctx-rh-superiorparietal",    [20, 180, 140]],
    [2030,    "ctx-rh-superiortemporal",    [140, 220, 220]],
    [2031,    "ctx-rh-supramarginal",    [80, 160, 20]],
    [2032,    "ctx-rh-frontalpole",    [100, 0, 100]],
    [2033,    "ctx-rh-temporalpole",    [70, 70, 70]],
    [2034,    "ctx-rh-transversetemporal",    [150, 150, 200]],
    [2035,    "ctx-rh-insula",    [255, 192, 32]]]

    FSnumbers = [x[0] for x in FreeSurferColorLUT]
    FSnames = [x[1] for x in FreeSurferColorLUT]
    #FScolors = [x[2] for x in FreeSurferColorLUT]
"""
