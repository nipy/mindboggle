file = '/hd2/Archive/registration_evaluation_2011_output/Results/dice_jacc_overlaps_fundi_brain_visa.txt'
file = '/hd2/Archive/registration_evaluation_2011_output/Results/dice_jacc_overlaps.txt'
#file = '/hd2/Archive/registration_evaluation_2011_output/Results/dice_jacc_overlaps_pits_forrest_bao.txt'
#file = '/hd2/Archive/registration_evaluation_2011_output/Results/dice_jacc_overlaps_ribbons_brain_visa.txt'

#file = '/hd2/Archive/registration_evaluation_2011_output/Results/avg_min_distance_fundi_forrest_bao_to_pits_kiho_im.txt'

field = 1
#field = 4

import numpy as np

f = open(file,'r')
T = f.readlines()
avg = 0
for row in T:
    avg += np.float(row.split(' ')[field])

avg = avg/len(T)
print(avg)
