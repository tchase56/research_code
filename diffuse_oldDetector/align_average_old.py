# Align old data

import numpy as np
# Functions for aligning using cross correlaton
from skimage.feature import register_translation
from scipy.ndimage.interpolation import shift
# User defined function
import Region_Of_Interest_2

# Load and save addresses
load_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\old_data\\'
save_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\old_data\\'
ROI_width = 64      # Make sure this is a power of 2 (should be slightly faster)
pixel_accuracy = 100      # Align with accuracy of 1/(pixel_accuracy)

# Load pre-aligned sif files
scan = []
scan_static = []
for i in range(211,218):
    scan.append(np.loadtxt(load_address + str(i) + '\imageSumPump'))
    scan_static.append(np.loadtxt(load_address + str(i) + '\imageSumStatic'))
scan = np.array(scan)
scan_static = np.array(scan_static)


# Save unaligned files as a numpy array    
np.save(save_address + 'averaged', scan)

# Choose ROIs for 4 lowest order peaks for alignment 
peak_region, background_region = Region_Of_Interest_2.GetRegionOfInterest(scan[0], 4, halfLength = ROI_width/2, contrastFactor = 0.25)

scan_aligned = []
scan_aligned_static = []
scan_aligned.append(scan[0])
scan_aligned_static.append(scan_static[0])
for k in range(1,np.shape(scan)[0]):
    offset_list = [register_translation(scan[0][int(peak_region[l][0]):int(peak_region[l][1]), int(peak_region[l][2]): int(peak_region[l][3])], 
                                        scan[k][int(peak_region[l][0]):int(peak_region[l][1]), int(peak_region[l][2]): int(peak_region[l][3])], 
                                        pixel_accuracy)[0] for l in range(0,4)]
    offset = np.average(offset_list, 0)
    scan_aligned.append(shift(scan[k], offset))
    scan_aligned_static.append(shift(scan_static[k], offset))
scan_aligned = np.array(scan_aligned)
scan_aligned_static = np.array(scan_aligned_static)
    
# Save aligned files as a numpy array
np.save(save_address + 'averaged_aligned', scan_aligned)
np.save(save_address + 'averaged_aligned_static', scan_aligned_static)

