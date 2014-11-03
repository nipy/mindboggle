#!/usr/bin/env python
"""
Brain label numbers, names, and colormap related to the DKT labeling protocol
==============================================================================

For more information about the Desikan-Killiany-Tourville cortical labeling
protocol see http://mindboggle.info/data and the article:

http://www.frontiersin.org/Brain_Imaging_Methods/10.3389/fnins.2012.00171/full
"101 labeled brain images and a consistent human cortical labeling protocol"
Arno Klein, Jason Tourville. Frontiers in Brain Imaging Methods. 6:171.
DOI: 10.3389/fnins.2012.00171

Calls return_numbers_names_colors() to return numbers, names, and colors
extracted from FreeSurfer's FreeSurferColorLUT.txt lookup table file
representing anatomical brain regions.
(See extract_numbers_names_colors() in LUT.py for extraction code.)

------------------------------------------------------------------------------
Combined/eliminated regions:
------------------------------------------------------------------------------
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

Jason Tourville:  "Regarding the lack of a full, consistent sulcal anterior
boundary for the inferior frontal gyrus:
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
rely on the margins of sulcal banks."

(5) Parahippocampal + entorhinal cortex + and lingual gyrus?

------------------------------------------------------------------------------
Regions bounded by sulcal fundi:
------------------------------------------------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Lateral surface:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- frontomarginal sulcus: [12,28]
- superior frontal: [3,28],[27,28]
- inferior frontal: [3,18],[3,19],[3,20], [18,27],[19,27],[20,27]
- precentral: [24,28]*, [[3,24],[24,27]]*, [[18,24],[19,24],[20,24]]*
- central sulcus: [22,24]
- postcentral: [22,29],[22,31], not:[22,24]
- intraparietal: [29,31], [8,29]
- primary intermediate sulcus /
    1st segment of the posterior superior temporal sulcus: [8,31]*
- sylvian fissure: [30,31]*, not:[18,30] (see note #2)
- lateral occipital sulcus: [8,11]*,[11,29]*
- anterior occipital sulcus: [11,15]*,[9,11]
- superior temporal sulcus: [15,30]
- inferior temporal sulcus: [9,15]*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PeriSylvian area (folds within the Sylvian fissure):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- circular sulcus: [12,35],[30,35],[34,35], [2,35],[10,35],[23,35],[26,35],
                 [22,35], [24,35], [31,35]
- 1st transverse temporal sulcus: [30,34]
- Heschl's sulcus: [30,34]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Medial surface:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- cingulate sulcus: [2,14],[10,14],[14,23],[14,26] (see note #3),
                  [2,28],[10,28],[23,28],[26,28],
                  [2,17],[10,17],[17,23],[17,26], [17,25]
- paracentral sulcus: [17,28]*
- parietooccipital fissure: [5,25]
- calcarine fissure: [13,25], [2,13],[10,13],[13,23],[13,26] not:[5,13] (note #4)
- superior rostral sulcus: [14,28]
- callosal sulcus: [2,4],[4,10],[4,23],[4,26]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ventral surface:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- lateral H-shaped orbital sulcus: [3,12],[12,27], [12,18],[12,19],[12,20]
- olfactory sulcus: [12,14]
- occipitotemporal sulcus: [7,9],[7,11]
- collateral sulcus: [6,7], [7,13], [7,16]

------------------------------------------------------------------------------
What boundaries will NEVER be derived by fundi, but instead by curvature, etc.
------------------------------------------------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Regions bounded by sulcal margins:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- interhemispheric fissure, dorsal margin:
    [17,28],[17,24],[17,22],[25,29],[5,29],[5,11]
- calcarine sulcus, dorsal margin: [5,21]
- calcarine sulcus, ventral margin: [21,13]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Regions with additional non-sulcal boundaries with subcortical regions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[16,6,9,30,12,14]

------------------------------------------------------------------------------
Notes:
------------------------------------------------------------------------------
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
    [left, right]_sulcus_label_pair_lists : list of two lists of lists of integer pairs
    [left, right]_cerebrum_cortex_[numbers, names, colors]
    [left, right]_cerebrum_noncortex_[numbers, names, colors]
    medial_cerebrum_noncortex_[numbers, names, colors]
    [left, right]_ventricle_[numbers, names, colors]
    medial_ventricle_[numbers, names, colors]
    extra_[numbers, names, colors]
    [left, right]_cerebellum_cortex_[numbers, names, colors]
    [left, right]_cerebellum_noncortex_[numbers, names, colors]
    medial_cerebellum_noncortex_[numbers, names, colors]
    [left, right]_ventricle_[numbers, names, colors]
    medial_ventricle_[numbers, names, colors]
    extra_[numbers, names, colors]
    [left_, right_]cerebrum_cortex_[numbers, names, colors][_DKT31, _DKT25]
    cerebrum_noncortex_[numbers, names, colors]
    [left, right]_cerebrum_[numbers, names, colors]
    cerebrum_[numbers, names, colors]
    [left, right]_cerebellum_[numbers, names, colors]
    cerebellum_cortex_[numbers, names, colors]
    cerebellum_noncortex_[numbers, names, colors]
    cerebellum_[numbers, names, colors]
    ventricle_[numbers, names, colors]
    label_[numbers, names, colors]
    colormap : list of lists
    colormap_normalized : list of lists

    Examples
    --------
    >>> from mindboggle.LABELS import DKTprotocol
    >>> dkt = DKTprotocol()
    >>> dkt.left_cerebrum_names

    """
    from mindboggle.LUT import return_numbers_names_colors

    #-------------------------------------------------------------------------
    # Return numbers, names, colors extracted from FreeSurferColorLUT.txt:
    #-------------------------------------------------------------------------
    numbers, names, colors = return_numbers_names_colors()

    #-------------------------------------------------------------------------
    # Cerebral cortex label numbers (31 + duplicates):
    #-------------------------------------------------------------------------
    left_cerebrum_cortex_list = \
        [3, 19, 20] + range(1000, 1004) + range(1005, 1036)
    right_cerebrum_cortex_list = \
        [42, 55, 56] + range(2000, 2004) + range(2005, 2036)

    #-------------------------------------------------------------------------
    # Cerebral ventricle label numbers:
    #-------------------------------------------------------------------------
    left_ventricle_list = [4, 5]
    right_ventricle_list = [43, 44]
    medial_ventricle_list = [14, 15, 72]

    #-------------------------------------------------------------------------
    # Cerebral noncortex label numbers (including ventricles above):
    #-------------------------------------------------------------------------
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
    left_cerebrum_noncortex_list = \
        [2, 9, 10, 11, 12, 13, 17, 18, 25, 26, 27, 28, 30, 31, 78, 91, 96] + \
        range(100, 109) + [155, 157] + range(550, 559) + [1004] + \
        range(3000, 3036) + [5001] + left_ventricle_list
    right_cerebrum_noncortex_list = \
        [41, 48, 49, 50, 51, 52, 53, 54, 57, 58, 59, 60, 62, 63, 79, 92, 97] + \
        range(109, 118) + [156, 158] + range(500, 509) + [2004] + \
        range(4000, 4036) + [5002] + right_ventricle_list
    medial_cerebrum_noncortex_list = medial_ventricle_list + \
                                        [192] + range(250, 256)

    #-------------------------------------------------------------------------
    # Cerebellar label numbers:
    #-------------------------------------------------------------------------
    # These labels [71, 72, 73] in Neuromorphometrics BrainCOLOR subcortex
    # did not have counterparts in FreeSurferColorLUT.txt, and were
    # reassigned to unused numbers in FreeSurferColorLUT.txt:
    #   [630, "cerebellar vermal lobules I-V"],
    #   [631, "cerebellar vermal lobules VI-VII"],
    #   [632, "cerebellar vermal lobules VIII-X"],
    #-------------------------------------------------------------------------
    left_cerebellum_cortex_list = [6, 8]
    right_cerebellum_cortex_list = [45, 47]
    left_cerebellum_noncortex_list = [7]
    right_cerebellum_noncortex_list = [46]
    medial_cerebellum_noncortex_list = [630, 631, 632]

    #-------------------------------------------------------------------------
    # Extra label numbers:
    #-------------------------------------------------------------------------
    #  [16, "Brain stem"],
    #  [24, "CSF"],
    #  [85, "optic chiasm"]]
    #  170-175: brain stem
    #-------------------------------------------------------------------------
    brainstem_list = [16] + range(170, 176)
    extra_list = [24, 85]

    #-------------------------------------------------------------------------
    # Label numbers:
    #-------------------------------------------------------------------------
    left_cerebrum_cortex_numbers = []
    right_cerebrum_cortex_numbers = []
    left_ventricle_numbers = []
    right_ventricle_numbers = []
    medial_ventricle_numbers = []
    left_cerebrum_noncortex_numbers = []
    right_cerebrum_noncortex_numbers = []
    medial_cerebrum_noncortex_numbers = []
    left_ventricle_numbers = []
    right_ventricle_numbers = []
    medial_ventricle_numbers = []
    left_cerebellum_cortex_numbers = []
    right_cerebellum_cortex_numbers = []
    left_cerebellum_noncortex_numbers = []
    right_cerebellum_noncortex_numbers = []
    medial_cerebellum_noncortex_numbers = []
    brainstem_numbers = []
    extra_numbers = []
    misc_numbers = []
    #-------------------------------------------------------------------------
    # Names corresponding to label numbers:
    #-------------------------------------------------------------------------
    left_cerebrum_cortex_names = []
    right_cerebrum_cortex_names = []
    left_ventricle_names = []
    right_ventricle_names = []
    medial_ventricle_names = []
    left_cerebrum_noncortex_names = []
    right_cerebrum_noncortex_names = []
    medial_cerebrum_noncortex_names = []
    left_ventricle_names = []
    right_ventricle_names = []
    medial_ventricle_names = []
    left_cerebellum_cortex_names = []
    right_cerebellum_cortex_names = []
    left_cerebellum_noncortex_names = []
    right_cerebellum_noncortex_names = []
    medial_cerebellum_noncortex_names = []
    brainstem_names = []
    extra_names = []
    misc_names = []
    #-------------------------------------------------------------------------
    # Colors corresponding to label numbers:
    #-------------------------------------------------------------------------
    left_cerebrum_cortex_colors = []
    right_cerebrum_cortex_colors = []
    left_ventricle_colors = []
    right_ventricle_colors = []
    medial_ventricle_colors = []
    left_cerebrum_noncortex_colors = []
    right_cerebrum_noncortex_colors = []
    medial_cerebrum_noncortex_colors = []
    left_ventricle_colors = []
    right_ventricle_colors = []
    medial_ventricle_colors = []
    left_cerebellum_cortex_colors = []
    right_cerebellum_cortex_colors = []
    left_cerebellum_noncortex_colors = []
    right_cerebellum_noncortex_colors = []
    medial_cerebellum_noncortex_colors = []
    brainstem_colors = []
    extra_colors = []
    misc_colors = []

    #-------------------------------------------------------------------------
    # Lists of numbers, names, and colors:
    #-------------------------------------------------------------------------
    for i, n in enumerate(numbers):
        if n in left_cerebrum_cortex_list:
            left_cerebrum_cortex_numbers.append(numbers[i])
            left_cerebrum_cortex_names.append(names[i])
            left_cerebrum_cortex_colors.append(colors[i])
        elif n in right_cerebrum_cortex_list:
            right_cerebrum_cortex_numbers.append(numbers[i])
            right_cerebrum_cortex_names.append(names[i])
            right_cerebrum_cortex_colors.append(colors[i])
        elif n in left_ventricle_list:
            left_ventricle_numbers.append(numbers[i])
            left_cerebrum_noncortex_numbers.append(numbers[i])
            left_ventricle_names.append(names[i])
            left_cerebrum_noncortex_names.append(names[i])
            left_ventricle_colors.append(colors[i])
            left_cerebrum_noncortex_colors.append(colors[i])
        elif n in right_ventricle_list:
            right_ventricle_numbers.append(numbers[i])
            right_cerebrum_noncortex_numbers.append(numbers[i])
            right_ventricle_names.append(names[i])
            right_cerebrum_noncortex_names.append(names[i])
            right_ventricle_colors.append(colors[i])
            right_cerebrum_noncortex_colors.append(colors[i])
        elif n in medial_ventricle_list:
            medial_ventricle_numbers.append(numbers[i])
            medial_cerebrum_noncortex_numbers.append(numbers[i])
            medial_ventricle_names.append(names[i])
            medial_cerebrum_noncortex_names.append(names[i])
            medial_ventricle_colors.append(colors[i])
            medial_cerebrum_noncortex_colors.append(colors[i])
        elif n in left_cerebrum_noncortex_list:
            left_cerebrum_noncortex_numbers.append(numbers[i])
            left_cerebrum_noncortex_names.append(names[i])
            left_cerebrum_noncortex_colors.append(colors[i])
        elif n in right_cerebrum_noncortex_list:
            right_cerebrum_noncortex_numbers.append(numbers[i])
            right_cerebrum_noncortex_names.append(names[i])
            right_cerebrum_noncortex_colors.append(colors[i])
        elif n in medial_cerebrum_noncortex_list:
            medial_cerebrum_noncortex_numbers.append(numbers[i])
            medial_cerebrum_noncortex_names.append(names[i])
            medial_cerebrum_noncortex_colors.append(colors[i])
        elif n in left_ventricle_list:
            left_ventricle_numbers.append(numbers[i])
            left_ventricle_names.append(names[i])
            left_ventricle_colors.append(colors[i])
        elif n in right_ventricle_list:
            right_ventricle_numbers.append(numbers[i])
            right_ventricle_names.append(names[i])
            right_ventricle_colors.append(colors[i])
        elif n in medial_ventricle_list:
            medial_ventricle_numbers.append(numbers[i])
            medial_ventricle_names.append(names[i])
            medial_ventricle_colors.append(colors[i])
        elif n in left_cerebellum_cortex_list:
            left_cerebellum_cortex_numbers.append(numbers[i])
            left_cerebellum_cortex_names.append(names[i])
            left_cerebellum_cortex_colors.append(colors[i])
        elif n in right_cerebellum_cortex_list:
            right_cerebellum_cortex_numbers.append(numbers[i])
            right_cerebellum_cortex_names.append(names[i])
            right_cerebellum_cortex_colors.append(colors[i])
        elif n in left_cerebellum_noncortex_list:
            left_cerebellum_noncortex_numbers.append(numbers[i])
            left_cerebellum_noncortex_names.append(names[i])
            left_cerebellum_noncortex_colors.append(colors[i])
        elif n in right_cerebellum_noncortex_list:
            right_cerebellum_noncortex_numbers.append(numbers[i])
            right_cerebellum_noncortex_names.append(names[i])
            right_cerebellum_noncortex_colors.append(colors[i])
        elif n in medial_cerebellum_noncortex_list:
            medial_cerebellum_noncortex_numbers.append(numbers[i])
            medial_cerebellum_noncortex_names.append(names[i])
            medial_cerebellum_noncortex_colors.append(colors[i])
        elif n in brainstem_list:
            brainstem_numbers.append(numbers[i])
            brainstem_names.append(names[i])
            brainstem_colors.append(colors[i])
        elif n in extra_list:
            extra_numbers.append(numbers[i])
            extra_names.append(names[i])
            extra_colors.append(colors[i])
        else:
            misc_numbers.append(numbers[i])
            misc_names.append(names[i])
            misc_colors.append(colors[i])

    #-------------------------------------------------------------------------
    # Aggregate lists of numbers, names, and colors:
    #-------------------------------------------------------------------------
    ventricle_numbers = left_ventricle_numbers + right_ventricle_numbers + \
                        medial_ventricle_numbers
    ventricle_names = left_ventricle_names + right_ventricle_names + \
                      medial_ventricle_names
    ventricle_colors = left_ventricle_colors + right_ventricle_colors + \
                       medial_ventricle_colors

    cerebrum_cortex_numbers = left_cerebrum_cortex_numbers + \
                              right_cerebrum_cortex_numbers
    cerebrum_cortex_names = left_cerebrum_cortex_names + \
                            right_cerebrum_cortex_names
    cerebrum_cortex_colors = left_cerebrum_cortex_colors + \
                             right_cerebrum_cortex_colors
    cerebrum_noncortex_numbers = left_cerebrum_noncortex_numbers + \
                                 right_cerebrum_noncortex_numbers + \
                                 medial_cerebrum_noncortex_numbers + \
                                 ventricle_numbers
    cerebrum_noncortex_names = left_cerebrum_noncortex_names + \
                               right_cerebrum_noncortex_names + \
                               medial_cerebrum_noncortex_names + \
                               ventricle_names
    cerebrum_noncortex_colors = left_cerebrum_noncortex_colors + \
                                right_cerebrum_noncortex_colors + \
                                medial_cerebrum_noncortex_colors + \
                                ventricle_colors
    left_cerebrum_numbers = left_cerebrum_cortex_numbers + \
                            left_cerebrum_noncortex_numbers
    left_cerebrum_names = left_cerebrum_cortex_names + \
                          left_cerebrum_noncortex_names
    left_cerebrum_colors = left_cerebrum_cortex_colors + \
                           left_cerebrum_noncortex_colors
    right_cerebrum_numbers = right_cerebrum_cortex_numbers + \
                             right_cerebrum_noncortex_numbers
    right_cerebrum_names = right_cerebrum_cortex_names + \
                           right_cerebrum_noncortex_names
    right_cerebrum_colors = right_cerebrum_cortex_colors + \
                            right_cerebrum_noncortex_colors
    cerebrum_numbers = left_cerebrum_numbers + right_cerebrum_numbers
    cerebrum_names = left_cerebrum_names + right_cerebrum_names
    cerebrum_colors = left_cerebrum_colors + right_cerebrum_colors

    left_cerebellum_numbers = left_cerebellum_cortex_numbers + \
                              left_cerebellum_noncortex_numbers
    left_cerebellum_names = left_cerebellum_cortex_names + \
                            left_cerebellum_noncortex_names
    right_cerebellum_numbers = right_cerebellum_cortex_numbers + \
                               right_cerebellum_noncortex_numbers
    right_cerebellum_names = right_cerebellum_cortex_names + \
                             right_cerebellum_noncortex_names

    cerebellum_cortex_numbers = left_cerebellum_cortex_numbers + \
                                right_cerebellum_cortex_numbers
    cerebellum_cortex_names = left_cerebellum_cortex_names + \
                              right_cerebellum_cortex_names
    cerebellum_cortex_colors = left_cerebellum_cortex_colors + \
                               right_cerebellum_cortex_colors
    cerebellum_noncortex_numbers = left_cerebellum_noncortex_numbers + \
                                   right_cerebellum_noncortex_numbers + \
                                   medial_cerebellum_noncortex_numbers
    cerebellum_noncortex_names = left_cerebellum_noncortex_names + \
                                 right_cerebellum_noncortex_names + \
                                 medial_cerebellum_noncortex_names
    cerebellum_noncortex_colors = left_cerebellum_noncortex_colors + \
                                  right_cerebellum_noncortex_colors
    cerebellum_numbers = cerebellum_cortex_numbers + \
                         cerebellum_noncortex_numbers
    cerebellum_names = cerebellum_cortex_names + cerebellum_noncortex_names
    cerebellum_colors = cerebellum_cortex_colors + cerebellum_noncortex_colors

    label_numbers = cerebrum_numbers + cerebellum_numbers + \
                    brainstem_numbers + extra_numbers
    label_names = cerebrum_names + cerebellum_names + \
                  brainstem_names + extra_names
    label_colors = cerebrum_colors + cerebellum_colors + \
                   brainstem_colors + extra_colors

    #-------------------------------------------------------------------------
    # Colormap:
    #-------------------------------------------------------------------------
    colormap = []
    colormap_normalized = []
    for i, x in enumerate(label_colors):
        colormap.append([label_numbers[i], 1, x[0], x[1], x[2]])
        colormap_normalized.append([label_numbers[i], 1,
                                    x[0]/255.0, x[1]/255.0, x[2]/255.0])

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
    sulcus_numbers = range(len(sulcus_names))

    sulcus_names_abbr = [
        "fms",
        "sfrs",
        "ifrs",
        "prcs",
        "cs",
        "pocs",
        "itps",
        "pis/csts1",
        "ls",
        "locs",
        "aocs",
        "sts",
        "its",
        "crs",
        "ftts/hs",
        "cgs",
        "pcs",
        "pos",
        "ccs",
        "sros",
        "cas",
        "lhos",
        "olfs",
        "ots",
        "cos"]

    #-------------------------------------------------------------------------
    # Lists of label pairs that define sulcus boundaries (or fundi)
    # according to the DKT labeling protocol.
    # 1000 [lh] or 2000 [rh] are added to these numbers below.
    #-------------------------------------------------------------------------
    pair_lists = [
        [[12, 28]],
        [[3, 28], [27, 28]],
        [[3, 18], [3, 19], [3, 20], [18, 27], [19, 27], [20, 27]],
        [[24, 28], [3, 24], [24, 27], [18, 24], [19, 24], [20, 24]],
        [[22, 24]],
        [[22, 29], [22, 31]],
        [[29, 31], [8, 29]],
        [[8, 31]],
        [[30, 31]],
        [[8, 11], [11, 29]],
        [[11, 15], [9, 11]],
        [[15, 30]],
        [[9, 15]],
        [[12, 35], [30, 35], [34, 35], [2, 35], [10, 35], [23, 35], [26, 35],
         [22, 35], [24, 35], [31, 35]],
        [[30, 34]],
        [[2, 14], [10, 14], [14, 23], [14, 26], [2, 28], [10, 28], [23, 28],
         [26, 28],
         [2, 17], [10, 17], [17, 23], [17, 26], [17, 25]],
        [[17, 28]],
        [[5, 25]],
        [[13, 25], [2, 13], [10, 13], [13, 23], [13, 26]],
        [[14, 28]],
        [[2, 4], [4, 10], [4, 23], [4, 26]],
        [[3, 12], [12, 27], [12, 18], [12, 19], [12, 20]],
        [[12, 14]],
        [[7, 9], [7, 11]],
        [[6, 7], [7, 16], [7, 13]]]

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

    print(
        '<ColorMap name="Mindboggle Colormap" space="RGB" indexedLookup="false">')
    print('  <NaN r="0" g="0" b="0"/>')
    print('  <Point x="-1" o="0"  r="0" g="0" b="0"/>')
    for row in colormap:
        print('  <Point x="{0}" o="1"  '
              'r="{1:1.2f}" g="{2:1.2f}" b="{3:1.2f}"/>'.
              format(row[0], row[2], row[3], row[4]))
    print('</ColorMap>')


"""
    #-------------------------------------------------------------------------
    # DKT cerebral cortical labeling protocol -- 25 labels:
    #-------------------------------------------------------------------------
    # Region numbers:
    # DKT31 to DKT25: [[10,23,26,27,19,20], [2,2,2,3,18,18]]
    left_cerebrum_cortex_numbers_DKT25 = range(1002, 1036)
    for n in [1004, 1010, 1019, 1020, 1023, 1026, 1027, 1032, 1033]:
        left_cerebrum_cortex_numbers_DKT25.remove(n)
    right_cerebrum_cortex_numbers_DKT25 = range(2002, 2036)
    for n in [2004, 2010, 2019, 2020, 2023, 2026, 2027, 2032, 2033]:
        right_cerebrum_cortex_numbers_DKT25.remove(n)
    cerebrum_cortex_numbers_DKT25 = left_cerebrum_cortex_numbers_DKT25 + \
                                    right_cerebrum_cortex_numbers_DKT25
    # Consolidate region labels:
    left_cerebrum_cortex_names_DKT25 = []
    for ilabel, label_number in enumerate(left_cerebrum_cortex_numbers):
        if label_number in left_cerebrum_cortex_numbers_DKT25:
            if label_number == 1002:
                left_cerebrum_cortex_names_DKT25.append('left cingulate')
            elif label_number == 1003:
                left_cerebrum_cortex_names_DKT25.append('left middle frontal')
            elif label_number == 1018:
                left_cerebrum_cortex_names_DKT25.append('left inferior frontal')
            else:
                left_cerebrum_cortex_names_DKT25. \
                    append(left_cerebrum_cortex_names[ilabel])
    right_cerebrum_cortex_names_DKT25 = []
    for ilabel, label_number in enumerate(right_cerebrum_cortex_numbers):
        if label_number in right_cerebrum_cortex_numbers_DKT25:
            if label_number == 2002:
                right_cerebrum_cortex_names_DKT25.append('right cingulate')
            elif label_number == 2003:
                right_cerebrum_cortex_names_DKT25.append('right middle frontal')
            elif label_number == 2018:
                right_cerebrum_cortex_names_DKT25.append('right inferior frontal')
            else:
                right_cerebrum_cortex_names_DKT25. \
                    append(right_cerebrum_cortex_names[ilabel])
    cerebrum_cortex_names_DKT25 = \
        left_cerebrum_cortex_names_DKT25 + right_cerebrum_cortex_names_DKT25
"""
"""
    # 2-D cortical map (UNDER CONSTRUCTION):
    #
    # Each gyrus has one or more sulcus boundaries -- rough breakdown:
    #
    # labels:  28  3   18  30  15  9   7   16   17    12   14    29  8      24   22      5    25   6   21    2    35   31
    #
    # Gyrii:  sFG mFG iFG sTG mTG iTG Fus pHip pCen  lOrb mOrb  sPL iPL    prCG poCG    Cun prCun ent pCalc Cing Ins sMarg
    #         ---|---|---|---|---|---|---|----|----  ---- ----  ---|---  -|----|----|-  ---|----- --- ----- ---- --- -----
    # Sulci:    sFS iFS Syl sTS iTS TO? Coll?        OrbF OrbT    iPS    prCS  CS poCS                 Cing Ins
    #
    # indices:   1   2   3   4   5   6   7      8    9      10      11   12  13   14  15  16   17  18
"""