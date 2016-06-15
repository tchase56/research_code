'''
Tyler Chase
06/14/2016
align and average multiple scans together 
aligning by 6 innermost bragg peaks
'''

import numpy as np
# Functions for aligning using cross correlaton
from skimage.feature import register_translation
from scipy.ndimage.interpolation import shift
# User defined function
import Region_Of_Interest_2




''' Values to change for each run'''
save_address = r'C:\Users\tchase56\Documents\UED\TaS2\20160611\Data\timescan\\'
load_address = r'C:\Users\tchase56\Documents\UED\TaS2\20160611\Data\timescan\\'
# Assumes files follow naming convention scan1 for example and the scan folders are in the same directory
scans = [9, 10, 11, 12, 13, 14, 15, 16]
pixel_accuracy = 100
ROI_width = 32      # Make sure this is a power of 2










# Load aligned files from each scan number listed above
run = []
for i in scans:
    load_address_2 = load_address + 'scan' + str(i) + '\images-ANDOR1\\'
    run.append(np.load(load_address_2 + 'averaged_aligned.npy'))
    
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(run[0][0,:,:] ,6, halfLength=ROI_width/2, contrastFactor = 0.25)


'''
# Align them with respect to eachother and average 
run_aligned = []    
run_aligned.append(run[0])
for i in range(1, len(run)):
    offset = register_translation(run[0][0,:,:], run[i][0,:,:], pixel_accuracy)[0]
    for j in range(0, np.shape(run[0])[0]):
        run[i][j,:,:] = shift(run[i][j,:,:], offset)
    run_aligned.append(run[i])
runs_aligned_averaged = np.average(run_aligned,0)

np.save(save_address + 'averaged_runs',runs_aligned_averaged)
'''


# Align them with respect to eachother and average 
run_aligned = []    
run_aligned.append(run[0])
for i in range(1, len(run)):
    offset_list = [register_translation(run[0][0, peak_region[l][0]:peak_region[l][1], peak_region[l][2]: peak_region[l][3]], 
                                            run[i][0, peak_region[l][0]:peak_region[l][1], peak_region[l][2]: peak_region[l][3]], 
                                            pixel_accuracy)[0] for l in range(0,6)]
    offset = np.average(offset_list,0)
    for j in range(0, np.shape(run[0])[0]):
        run[i][j,:,:] = shift(run[i][j,:,:], offset)
    run_aligned.append(run[i])
runs_aligned_averaged = np.average(run_aligned,0)

np.save(save_address + 'averaged_runs',runs_aligned_averaged)

           
