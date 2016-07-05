'''
Tyler Chase
04/28/2016
Group two time series, Rotate diffraction pattern, normalize, center, symmetrize
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
ROI_width = 15      # Make sure this is a power of 2
# Addresses for loading and saving
load_address = r'E:\tchase56\TaS2\20160611\TimeScan\scan6\images-ANDOR1\\'
save_address = r'E:\tchase56\TaS2\20160611\TimeScan\scan6\images-ANDOR1\\'
temp = np.load(load_address + 'averaged_aligned.npy')
peaks_saved = 0


def mask_1(image, radius_inner, radius_outer, ave_x, ave_y):
    mask = np.zeros(np.shape(image))
    for i in range(0, np.shape(image)[0]):
        for j in range(0, np.shape(image)[1]):
            if np.sqrt((i-ave_y)**2 + (j-ave_x)**2)>radius_inner:
                mask[i,j] = 1
            if np.sqrt((i-ave_y)**2 + (j-ave_x)**2)>radius_outer:
                mask[i,j] = 0
    return(mask)
    
def mask_3(image, radius_inner, ave_x, ave_y):
    mask = np.zeros(np.shape(image))
    for i in range(0, np.shape(image)[0]):
        for j in range(0, np.shape(image)[1]):
            if np.sqrt((i-ave_y)**2 + (j-ave_x)**2)<radius_inner:
                mask[i,j] = 1
    return(mask)
    
def mask_4(image, radius_outer, ave_x, ave_y):
    mask = np.ones(np.shape(image))
    for i in range(0, np.shape(image)[0]):
        for j in range(0, np.shape(image)[1]):
            if np.sqrt((i-ave_y)**2 + (j-ave_x)**2)<radius_outer:
                mask[i,j] = 0
    return(mask)
'''   
def mask_2(image, peak_radius, list_x, list_y):
    mask = np.ones(np.shape(image))
    for k in range(0, np.shape(list_x)[0]):
        for i in range(0, np.shape(image)[0]):
            for j in range(0, np.shape(image)[1]):
                if np.sqrt((i-list_y[k])**2 + (j-list_x[k])**2)<peak_radius:
                    mask[i,j] = 0
    return(mask)
'''    
def mask_2(image, peak_radius, list_x, list_y):
    list_x = np.array(list_x)
    list_y = np.array(list_y)
    mask = np.ones(np.shape(image))
    for i in range(0, np.shape(image)[0]):
        for j in range(0, np.shape(image)[1]):
            if np.any(((i-list_y)**2 + (j-list_x)**2)<peak_radius**2):
                mask[i,j] = 0
    return(mask)


# Load in the funs and choose ROIs for centering
runs = []
for i in range(0, np.shape(temp)[0]):
    runs.append(temp[i,:,:])


rotated = runs

'''
# Center the image
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs[0] ,6, halfLength=(ROI_width/2), contrastFactor = 0.25)
image_center = np.shape(temp[0,:,:])[0]/2 + 0.5 -1
diffraction_shift = []
diffraction_center = []
for i in range(0,np.shape(runs)[0]):
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs[i], peak_region, background_region, halfLength=ROI_width/2)
    y_diffraction = np.mean([peak_fit_pump[0][0][1], peak_fit_pump[1][0][1], peak_fit_pump[2][0][1], peak_fit_pump[3][0][1],
                             peak_fit_pump[4][0][1], peak_fit_pump[5][0][1]])
    x_diffraction = np.mean([peak_fit_pump[0][1][1], peak_fit_pump[1][1][1], peak_fit_pump[2][1][1], peak_fit_pump[3][1][1],
                             peak_fit_pump[4][1][1], peak_fit_pump[5][1][1]])
    diffraction_shift.append([-1*(x_diffraction-image_center),-1*(y_diffraction-image_center)])
    diffraction_center.append([x_diffraction,y_diffraction])
# Shift is about 5 pixels off
runs_centered = []    
for i in range(0,np.shape(runs)[0]):
    runs_centered.append(shift(runs[i],diffraction_shift[i]))
'''    
'''
# Rotate the image
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs_centered[0] ,2, halfLength=(ROI_width/2), contrastFactor = 0.25)
[peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs_centered[0], peak_region, background_region, halfLength=ROI_width/2)
deltaY = peak_fit_pump[1][1][1] - peak_fit_pump[0][1][1]
deltaX = peak_fit_pump[1][0][1] - peak_fit_pump[0][0][1]
angle = (360/(2*np.pi))*np.arctan(deltaY/deltaX)
rotated = []
for i in range(0,np.shape(runs)[0]):
    rotated.append(rotate(runs_centered[i], angle))
'''    
# Normalize the image
radius_inner = 60
radius_outer = 395
mask = np.ones(np.shape(rotated[0]))
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(rotated[0] ,6, halfLength=(ROI_width/2), contrastFactor = 0.05)
[peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(rotated[0], peak_region, background_region, halfLength=ROI_width/2)
ave_x = np.average([peak_fit_pump[0][0][1], peak_fit_pump[1][0][1], peak_fit_pump[2][0][1], peak_fit_pump[3][0][1], peak_fit_pump[4][0][1], peak_fit_pump[5][0][1]])
ave_y = np.average([peak_fit_pump[0][1][1], peak_fit_pump[1][1][1], peak_fit_pump[2][1][1], peak_fit_pump[3][1][1], peak_fit_pump[4][1][1], peak_fit_pump[5][1][1]])
mask = mask_1(rotated[0], radius_inner, radius_outer, ave_x, ave_y)
mask_throughbeam = mask_3(rotated[0], radius_inner, ave_x, ave_y)

peak_radius = 15
list_x = []
list_y = []
if peaks_saved == 0:
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(mask*rotated[0] ,5, halfLength=(ROI_width/2), contrastFactor = 0.05)
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(rotated[0], peak_region, background_region, halfLength=ROI_width/2)
    np.save(save_address + 'peak_region.npy', peak_region)
    np.save(save_address + 'peak_fit_pump.npy', peak_fit_pump)
else:
    peak_region = np.load(save_address + 'peak_region.npy')
    peak_fit_pump = np.load(save_address + 'peak_fit_pump.npy')

for i in range(0,np.shape(peak_region)[0]):
    list_x.append(peak_fit_pump[i][0][1])
    list_y.append(peak_fit_pump[i][1][1])
mask_peaks = mask_2(mask*rotated[0], peak_radius, list_x, list_y)

mask_detectorbkgnd = mask_4(rotated[0], 600, ave_x, ave_y)



# Save Mask
np.save(save_address + 'mask_detector_background.npy', mask_detectorbkgnd)
np.save(save_address + 'boundary_mask.npy', mask)
np.save(save_address + 'throughbeam_mask.npy', mask_throughbeam)
np.save(save_address + 'peak_mask.npy', mask_peaks)
np.save(save_address + 'aligned_rotated.npy', rotated)

plt.figure()
plt.subplot(121)
plt.imshow(mask_peaks*mask*rotated[0]).set_clim(0,2500)
plt.subplot(122)
plt.imshow(rotated[0]).set_clim(0,2500)
plt.show()






'''
mask[:int(high_q_min_y-40), :] = 0
mask[int(high_q_max_y+40):, :] = 0
mask[:, :int(high_q_min_x-40)] = 0
mask[:, int(high_q_max_x+40):] = 0
mask[int(low_q_min_y+40):int(low_q_max_y-40), int(low_q_min_x+40):int(low_q_max_x-40)] = 0
norm = []
norm.append(rotated[0])
reference = np.sum(mask*rotated[0])
for i in range(1,np.shape(runs)[0]):
    norm_factor = reference/np.sum(mask*rotated[i])
    norm.append(rotated[i]*norm_factor)
    


# Save grouped, normalized, centered, symmetrized, and un-symmetrized
np.save(save_address + 'runs_grouped_rotated_normalized_centered.npy', norm)
'''

