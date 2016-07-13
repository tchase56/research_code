'''
Tyler Chase
06/30/2016
Center the diffraction patterns
Rotate the image to align with screen (just looks better)
Symmetrize by rotation

Instructions
0. Make sure you have already run align_average.py or align_different_scans.py
   (set different_scans flag accordingly)
1. Change addresses
2. Run the code
3. Click nbpeaks innermost peaks for centering the images
4. Click background for nbpeaks innermost peaks
5. Click two peaks that you would like to be rotated to have a horizontal line
   for the vector connecting them for rotating the images (mostly for show)
6. Click background for those two peaks
'''

import matplotlib.pyplot as plt
import numpy as np
import Region_Of_Interest_2
import Fit_Peaks_2
from scipy.ndimage.interpolation import shift
from scipy.ndimage.interpolation import rotate





''' Values to change for each run''' 
# Pixel accuracy with wich to align
pixel_accuracy = 100      # Align with accuracy of 1/(pixel_accuracy)
ROI_width = 64      # Make sure this is a power of 2
# Addresses for loading and saving
load_address = r'C:\Users\tchase56\Documents\UED\Ni\data\20160628\long\\'
save_address = r'C:\Users\tchase56\Documents\UED\Ni\data\20160628\long\\'
# Is this the output of align_average or the output of align_different_scans
different_scans = 1
# Number of innermost peaks that you would like to align with respect to 
nbpeaks = 4
# What fold of rotation symmetry is the pattern (ex. 4-fold symmetry for FCC)
fold = 4










# Function for symmetrizing an image with six fold rotation symmetry
def symmetrize(image, fold):
    rotated = []
    for i in range(fold):
        rotated.append(image)
        image = rotate(image, 360/fold, reshape = False)
    symmetrized = np.average(rotated, 0)
    return(symmetrized)

# Read in the timescan you want to symmetrize
if different_scans == 0:
    temp = np.load(load_address + 'averaged_aligned.npy')
else :
    temp = np.load(load_address + 'averaged_runs.npy')
        
# Load in the runs
runs = []
for i in range(0, np.shape(temp)[0]):
    runs.append(temp[i,:,:])

# Center the image (using the average postion of the six innermost peaks)
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs[0] ,nbpeaks, halfLength=(ROI_width/2), contrastFactor = 0.25)
image_center = np.shape(temp[0,:,:])[0]/2 + 0.5 -1
diffraction_shift = []
diffraction_center = []
for i in range(0,np.shape(runs)[0]):
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs[i], peak_region, background_region, halfLength=ROI_width/2)
    #y_diffraction = np.mean([peak_fit_pump[0][0][1], peak_fit_pump[1][0][1], peak_fit_pump[2][0][1], peak_fit_pump[3][0][1]])
    #x_diffraction = np.mean([peak_fit_pump[0][1][1], peak_fit_pump[1][1][1], peak_fit_pump[2][1][1], peak_fit_pump[3][1][1]])
    y_diffraction = np.mean([peak_fit_pump[l][0][1] for l in range(nbpeaks)])
    x_diffraction = np.mean([peak_fit_pump[l][1][1] for l in range(nbpeaks)])
    diffraction_shift.append([(image_center-x_diffraction),(image_center-y_diffraction)])
    diffraction_center.append([x_diffraction,y_diffraction])
runs_centered = []    
for i in range(0,np.shape(runs)[0]):
    runs_centered.append(shift(runs[i],diffraction_shift[i]))    

# Rotate the image
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs_centered[0] ,2, halfLength=(ROI_width/2), contrastFactor = 0.25)
[peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs_centered[0], peak_region, background_region, halfLength=ROI_width/2)
deltaY = peak_fit_pump[1][1][1] - peak_fit_pump[0][1][1]
deltaX = peak_fit_pump[1][0][1] - peak_fit_pump[0][0][1]
angle = (360/(2*np.pi))*np.arctan(deltaY/deltaX)
rotated = []
for i in range(0,np.shape(runs)[0]):
    rotated.append(rotate(runs_centered[i], angle))
    
# Symmetrize image
symmetric = []
for i in range(np.shape(rotated)[0]):
    symmetric.append(np.array(symmetrize(rotated[i], fold)))

# Plot symmetrized images    
for i in range(np.shape(symmetric)[0]):
    plt.figure()
    plt.imshow(symmetric[i])
    plt.show()
 
# Save symmetrized images   
np.save(save_address + 'symmetrized.npy', symmetric)