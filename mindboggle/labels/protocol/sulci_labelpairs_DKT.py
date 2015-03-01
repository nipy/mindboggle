"""
Lists of label pairs that define sulcus boundaries (or fundi)
according to the Desikan-Killiany-Tourville cortical labeling protocol.

Authors:
    - Jason Tourville  (jtour@bu.edu)
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
def sulcus_boundaries():

    label_pair_lists = [[[12,28]],
                        [[3,28]],
                        [[3,18]],
                        [[24,28], [3,24], [18,24]],
                        [[22, 24]],
                        [[22,29], [22,31]],
                        [[29,31], [8,29]],
                        [[8,31]],
                        [[30,31]],
                        [[8,11], [11,29]],
                        [[11,15], [9,11]],
                        [[15,30]],
                        [[9,15]],
                        [[30,35], [34,35], [12,35], [2,35], [24,35], [22,35], [31,35]],
                        [[30,34]],
                        [[2,14], [2,28], [2,17], [17,25]],
                        [[17,28]],
                        [[5,25]],
                        [[13,25], [2,13]],
                        [[14,28]],
                        [[2,4]],
                        [[12,18], [3,12]],
                        [[12,14]],
                        [[7,9], [7,11]],
                        [[6,7], [7,16], [7,13]]]

    return label_pair_lists
