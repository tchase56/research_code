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
from expand_contract import peak_expansion_f
from expand_contract import expand





''' Values to change for each run''' 
# Pixel accuracy with wich to align
pixel_accuracy = 100      # Align with accuracy of 1/(pixel_accuracy)
ROI_width = 64      # Make sure this is a power of 2
# Addresses for loading and saving
load_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
load_address_2 = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\long_scan\\'
save_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
save_address_2 = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\long_scan\\'
# Is this the output of align_average or the output of align_different_scans
different_scans = 1
# Number of innermost peaks that you would like to align with respect to 
nbpeaks = 4
# What fold of rotation symmetry is the pattern (ex. 4-fold symmetry for FCC)
fold = 4
# Flag to account for peak expansion in diffraction pattern (adds expansion interpolation)
expansion_interpolation = 1
norm_flag = 1










# Function for symmetrizing an image with six fold rotation symmetry
def symmetrize(image, fold):
    rotated = []
    for i in range(fold):
        rotated.append(image)
        image = rotate(image, 360/fold, reshape = False)
    symmetrized = np.average(rotated, 0)
    return(symmetrized)
    
# Old method of symmetrization with flipping
def symmetrize_flip(image):
    temp = np.flipud(image)
    temp_2 = np.average([image, temp],0)
    temp_3 = np.fliplr(temp_2)
    symmetrized = np.average([temp_3,temp_2],0)
    return symmetrized
    
# Create mask and mask image
def make_mask(image, outside, inside):
    mask = np.zeros(np.shape(image))
    mask[outside:inside,outside:-1-outside] = 1
    mask[-1-inside:-1-outside,outside:-1-outside] = 1
    mask[outside:-1-outside, outside:inside] = 1
    mask[outside:-1-outside, -1-inside:-1-outside] = 1
    masked = image*mask
    return masked

# Read in the timescan you want to symmetrize
if different_scans == 0:
    temp = np.load(load_address + 'averaged_aligned.npy')
else :
    temp = np.load(load_address + 'averaged_runs.npy')
if different_scans == 0:
    temp_2 = np.load(load_address_2 + 'averaged_aligned.npy')
else :
    temp_2 = np.load(load_address_2 + 'averaged_runs.npy')
        
# Load in the runs in chronological order
# (-25, -5, 20, 45, 70, 95) + (-6, -2, 2, 6, 10, 14, 18) -> (-25, -6, -5, -2, 2, 6, 10, 14, 18, 20, 45, 70, 95)
runs = []
runs.append(temp_2[0,:,:])
runs.append(temp[0,:,:])
runs.append(temp_2[1,:,:])
for i in range(1, np.shape(temp)[0]):
    runs.append(temp[i,:,:])
for i in range(2, np.shape(temp_2)[0]):
    runs.append(temp_2[i,:,:])

# Center the image (using the average postion of the nbpeaks innermost peaks)
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs[0] ,nbpeaks, 
    halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Select Innermost Peaks for centering')
image_center = np.shape(temp[0,:,:])[0]/2 + 0.5 -1
diffraction_shift = []
diffraction_center = []
for i in range(0,np.shape(runs)[0]):
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs[i], peak_region, background_region, halfLength=ROI_width/2)
    y_diffraction = np.mean([peak_fit_pump[l][0][1] for l in range(nbpeaks)])
    x_diffraction = np.mean([peak_fit_pump[l][1][1] for l in range(nbpeaks)])
    diffraction_shift.append([(image_center-x_diffraction),(image_center-y_diffraction)])
    diffraction_center.append([x_diffraction,y_diffraction])
runs_centered = []    
for i in range(0,np.shape(runs)[0]):
    runs_centered.append(shift(runs[i],diffraction_shift[i])) 
    
# Rotate the image
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs_centered[0] ,2, 
    halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Select 2 peaks that you would like to lie flat in rotated image')
[peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs_centered[0], peak_region, background_region, halfLength=ROI_width/2)
deltaY = peak_fit_pump[1][1][1] - peak_fit_pump[0][1][1]
deltaX = peak_fit_pump[1][0][1] - peak_fit_pump[0][0][1]
angle = (360/(2*np.pi))*np.arctan(deltaY/deltaX)
rotated = []
for i in range(0,np.shape(runs)[0]):
    rotated.append(rotate(runs_centered[i], angle, reshape=False))
    
if norm_flag == 1:           
    # Normalize Screen Intensity
    masked = []
    norm = []
    screen_intensity = []
    screen_intensity_norm = []
    norm.append(rotated[0])
    for i in range(np.shape(rotated)[0]):
        masked.append(make_mask(rotated[i], 133, int(1024/2-66)))
    sum_0 = np.sum(masked[0])
    screen_intensity.append(sum_0)
    screen_intensity_norm.append(sum_0)
    for i in range(1,np.shape(rotated)[0]):
        sum_i = np.sum(masked[i])
        screen_intensity.append(sum_i)
        factor = sum_0/sum_i
        screen_intensity_norm.append(sum_i*factor)
        norm.append(rotated[i]*factor)
    rotated = norm #line added to not break following code after adding normalization
    
# Account for peak expansion with interpolation and renormalize
if expansion_interpolation == 1:
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(rotated[0] ,nbpeaks, 
        halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Select Innermost Peaks for expansion interpolation')
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(rotated[0], peak_region, background_region, halfLength=ROI_width/2)
    peak_distance_ref = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])
    peak_distance = []
    peak_distance.append(peak_distance_ref)
    peak_distance_expanded = []
    peak_distance_expanded.append(peak_distance_ref)
    runs_expanded = []
    runs_expanded.append(rotated[0])
    for i in range(1,np.shape(runs)[0]):
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(rotated[i], peak_region, background_region, halfLength=ROI_width/2)
        peak_distance.append(peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1]))
        factor = peak_distance[0]/peak_distance[i]
        expanded = expand(rotated[i],factor)
        runs_expanded.append(expanded)
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(expanded, peak_region, background_region, halfLength=ROI_width/2)
        peak_distance_expanded.append(peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1]))
    rotated = runs_expanded #line added to not break following code after adding expansion interpolation
      
# Symmetrize image
symmetric_0 = []
for i in range(np.shape(rotated)[0]):
    symmetric_0.append(np.array(symmetrize(rotated[i], fold)))
    
# Complete symmetrization
symmetric = []
for i in range(np.shape(symmetric_0)[0]):
    symmetric.append(np.array(symmetrize_flip(symmetric_0[i])))

# Plot symmetrized images    
for i in range(1,np.shape(symmetric)[0]):
    plt.figure()
    plt.imshow(symmetric[i]-symmetric[0]).set_clim(0,10)
    plt.show()
    
# Separate long and short runs for saving
symmetric_short = []
symmetric_long = []
rotated_short = []
rotated_long = []
symmetric_long.append(symmetric[0])
symmetric_long.append(symmetric[2])
symmetric_long.extend(symmetric[9:])
symmetric_short.append(symmetric[1])
symmetric_short.extend(symmetric[3:9])
rotated_long.append(rotated[0])
rotated_long.append(rotated[2])
rotated_long.extend(rotated[9:])
rotated_short.append(rotated[1])
rotated_short.extend(rotated[3:9])
 
if expansion_interpolation == 1 and norm_flag == 1:
    # Save symmetrized images   
    np.save(save_address + 'expanded_symmetrized_norm.npy', symmetric_short)
    np.save(save_address + 'expanded_norm.npy', rotated_short)
    np.save(save_address_2 + 'expanded_symmetrized_norm.npy', symmetric_long)
    np.save(save_address_2 + 'expanded_norm.npy', rotated_long)
    plt.figure()
    plt.title('Expansion Interpolation Test')
    plt.plot(peak_distance, label = 'original')
    plt.plot(peak_distance_expanded, label = 'expanded')
    #plt.ylim(308.251, 308.6)
    plt.legend()
    plt.show()
    plt.figure()
    plt.plot(screen_intensity)
    plt.plot(screen_intensity_norm)
    plt.title('Intensity Normalization')
    plt.show()
if expansion_interpolation == 1 and norm_flag == 0:
    # Save symmetrized images   
    np.save(save_address + 'expanded_symmetrized.npy', symmetric_short)
    np.save(save_address + 'expanded.npy', rotated_short)
    np.save(save_address_2 + 'expanded_symmetrized.npy', symmetric_long)
    np.save(save_address_2 + 'expanded.npy', rotated_long)
    plt.figure()
    plt.title('Expansion Interpolation Test')
    plt.plot(peak_distance, label = 'original')
    plt.plot(peak_distance_expanded, label = 'expanded')
    #plt.ylim(308.251, 308.6)
    plt.legend()
    plt.show()
if expansion_interpolation == 0  and norm_flag == 1:
    # Save symmetrized images   
    np.save(save_address + 'symmetrized_norm.npy', symmetric_short)
    np.save(save_address + 'norm.npy', rotated_short)
    np.save(save_address_2 + 'symmetrized_norm.npy', symmetric_long)
    np.save(save_address_2 + 'norm.npy', rotated_long)    
    plt.figure()
    plt.plot(screen_intensity)
    plt.plot(screen_intensity_norm)
    plt.title('Intensity Normalization')
    plt.show()
if expansion_interpolation == 0 and norm_flag == 0:
    # Save symmetrized images   
    np.save(save_address + 'symmetrized.npy', symmetric_short)
    np.save(save_address + 'centered_rotated.npy', rotated_short)
    np.save(save_address_2 + 'symmetrized.npy', symmetric_long)
    np.save(save_address_2 + 'centered_rotated.npy', rotated_long)   