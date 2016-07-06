# OLD DETECTOR
'''
Tyler Chase
06/30/2016
This module reads in (.tif) data from a scan from the UED beamline. 
Then it aligns and averages for each time delay, also aligning between time delays.
It saves an aligned and unaligned numpy array, averaged.npy and averaged_aligned.npy respectively.
saves numpy_array[delay_stage_index][pixel_rows][pixel_columns]

Instructions
1. Change addresses
2. Change delay stage values
3. Run the code
4. Click the nbpeaks innermost bragg peaks
5. Click next to the nbpeaks innermost bragg peaks (background)
'''

''' Values to change for each run'''
load_address = r'E:\tchase56\TaS2\20160705\TaS2_150K\scan5\images-ANDOR1\\'
save_address = r'E:\tchase56\TaS2\20160705\TaS2_150K\scan5\images-ANDOR1\\'
# Delay stage settings
delayStage_start, delayStage_end, delayStage_step  = [67.025, 67.4, 0.0075]
# Number of innermost peaks that you would like to align with respect to 
nbpeaks = 6
# Pixel accuracy with wich to align
pixel_accuracy = 100      # Align with accuracy of 1/(pixel_accuracy)
ROI_width = 32      # Make sure this is a power of 2 (should be slightly faster)










import numpy as np
import glob
import matplotlib.pyplot as plt
# Functions for aligning using cross correlaton
from skimage.feature import register_translation
from scipy.ndimage.interpolation import shift
# User defined function
import Region_Of_Interest_2

# Make vector of delay stage values
delay = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)

# Choose ROIs for nbpeaks lowest order peaks for alignment
temp = plt.imread(glob.glob(load_address + '*.tif')[0])
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(temp , nbpeaks, halfLength=ROI_width/2, contrastFactor = 0.25)

# Average unaligned scans and save as numpy arrays
scan_times = []
scan_times_aligned = []
for i in delay:
    scan = []
    # Load sum of images
    for j in glob.glob(load_address + '*0' + "%0.4f" % i + '*.tif'):
        scan.append(np.array(plt.imread(j)))
    print('delay stage ' + str(i))
    scan_times.append(np.average(np.array(scan),0))     # creating list of unaligned averaged images

    # Align Images and print percentage complete for each delay stage value
    percent_1 = 0  
    scan_2 = []
    scan_2.append(scan[0])
    for k in range(1,len(scan)):
        offset_list = [register_translation(scan[0][peak_region[l][0]:peak_region[l][1], peak_region[l][2]: peak_region[l][3]], 
                                            scan[k][peak_region[l][0]:peak_region[l][1], peak_region[l][2]: peak_region[l][3]], 
                                            pixel_accuracy)[0] for l in range(0,nbpeaks)]
        offset = np.average(offset_list,0)
        scan_2.append(shift(scan[k], offset))
        percent_2 = round(k/len(scan)*100)             
        if percent_2 > percent_1:        
            print(str(round(k/len(scan)*100)) + '%')
            percent_1 = percent_2
    scan_times_aligned.append(np.average(np.array(scan_2),0))     # creating list of aligned averaged images

# If there is more than one time delay align the time delays
if len(scan_times_aligned) > 1:
    scan_times_aligned_2 = []
    scan_times_aligned_2.append(scan_times_aligned[0])
    for l in range(1,len(scan_times_aligned)):
        offset_list = [register_translation(scan_times_aligned[0][peak_region[m][0]:peak_region[m][1], peak_region[m][2]: peak_region[m][3]], 
                                            scan_times_aligned[l][peak_region[m][0]:peak_region[m][1], peak_region[m][2]: peak_region[m][3]], 
                                            pixel_accuracy)[0] for m in range(0,nbpeaks)]
        offset = np.average(offset_list,0)
        scan_times_aligned_2.append(shift(scan_times_aligned[l], offset))
    # Save the aligned averages
    np.save(save_address + 'averaged_aligned', scan_times_aligned_2)
else:
    np.save(save_address + 'averaged_aligned', scan_times_aligned)

# Save the unaligned averages and the aligned averages (averaged[delay_stage][row pixels][column pixels])     
np.save(save_address + 'averaged', scan_times)
