# OLD DETECTOR
'''
Tyler Chase
06/30/2016
average multiple scans together, and align them with respect to eachother

Instructions
0. Make sure align_average.py has already been run for each scan
1. Change addresses
2. Change scan number values (this code is written with the naming convention "scan#")
3. Run the code
2. Click the nbpeaks innermost bragg peaks
3. Click next to the nbpeaks innermost bragg peaks (background)
'''

''' Values to change for each run'''
# Location of where you want averaged of scans saved
save_address = r'C:\Users\tchase56\Documents\UED\Ni\data\20160629\short\\'
# Location of list of scan folders
load_address = r'C:\Users\tchase56\Documents\UED\Ni\data\20160629\short\\'
scans = [59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 78, 80, 81, 82, 84, 85, 86]
# Number of innermost peaks that you would like to align with respect to 
nbpeaks = 4
pixel_accuracy = 100
ROI_width = 64      # Make sure this is a power of 2









import numpy as np
# Functions for aligning using cross correlaton
from skimage.feature import register_translation
from scipy.ndimage.interpolation import shift
# User defined function
import Region_Of_Interest_2

# Load aligned files from each scan number listed above
run = []
for i in scans:
    load_address_2 = load_address + 'scan' + str(i) + '\images-ANDOR1\\'
    run.append(np.load(load_address_2 + 'averaged_aligned.npy'))
# Choose ROIs for nbpeaks lowest order peaks for alignment
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(run[0][0,:,:] ,nbpeaks, halfLength=ROI_width/2, contrastFactor = 0.25)

# look at runs averaged for performance checking of alignment 
runs_averaged = np.average(run,0)

# Align them with respect to eachother and average 
run_aligned = []    
run_aligned.append(run[0])
for i in range(1, len(run)):
    offset_list = [register_translation(run[0][0, peak_region[l][0]:peak_region[l][1], peak_region[l][2]: peak_region[l][3]], 
                                            run[i][0, peak_region[l][0]:peak_region[l][1], peak_region[l][2]: peak_region[l][3]], 
                                            pixel_accuracy)[0] for l in range(0,nbpeaks)]
    offset = np.average(offset_list,0)
    for j in range(0, np.shape(run[0])[0]):
        run[i][j,:,:] = shift(run[i][j,:,:], offset)
    run_aligned.append(run[i])
runs_aligned_averaged = np.average(run_aligned,0)

# Save the averaged and aligned runs
np.save(save_address + 'averaged_runs',runs_aligned_averaged)

'''
plt.figure(figsize=(20,20))
for i in range(np.shape(runs_averaged)[0]):
    plt.subplot(121)
    plt.title(str(i))
    plt.imshow(runs_averaged[i][360:430, 390:460])
    plt.subplot(122)
    plt.title(str(i))
    plt.imshow(runs_aligned_averaged[i][360:430, 390:460])
    plt.pause(0.5)
'''
           
