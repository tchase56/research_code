'''
Tyler Chase
07/22/2016
This module reads in the output files and plots both a diffuse linescan and
shows the region of interest that creates this linescan

Instructions
0. Make sure you have already run align_average.py or align_different_scans.py
   (set different_scans flag accordingly [different_scans flag])
1. Change address of load and save
2. Change delay stage values and make sure time zero is correct
3. Choose weather you want to plot a horizontal ROI in addition to vertical
   (if yes set x_flag to 1)
'''







''' Values to change for each run'''
load_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
save_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
# Delay stage settings
#delayStage_start, delayStage_end, delayStage_step  = [37.8, 56.55, 3.75]
delayStage_start, delayStage_end, delayStage_step  = [41.4, 45.0, 0.6]
#delayStage_start, delayStage_end, delayStage_step  = [64.65, 65.25, 0.1]
#delayStage_start, delayStage_end, delayStage_step  = [64.65, 65.25, 0.1]
#delayStage_start, delayStage_end, delayStage_step  = [63.7, 64.3, 0.1]
time_zero = 42.3
#time_zero = 64.67
#time_zero = 63.8
ROI_width = 50      # Make sure this is a power of 2
# Flag for forming linescan along x
x_flag = 0
# Flag for different file outputs
# (0, alignment output)
# (1, averaging runs output)
# (2, symmetrized expanded normalized)
# (3, symmetrized expanded)
# (4, symmetrized normalized)
# (5, symmetrized)
different_scans = 3











import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import Region_Of_Interest_2
import Fit_Peaks_2

# Delay Stage Settings
delay_stage = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
speed_of_light = 299792458 # meters/second
delay = (np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12 # Delay in ps


# Load data that you would like to see diffuse scattering in
if different_scans == 0:
    runs = np.load(load_address + 'averaged_aligned.npy')
elif different_scans == 1:
    runs = np.load(load_address + 'averaged_runs.npy')
elif different_scans == 2:
    runs = np.load(load_address + 'expanded_symmetrized_norm.npy')
elif different_scans == 3:
    runs = np.load(load_address + 'expanded_symmetrized.npy')
elif different_scans == 4:
    runs = np.load(load_address + 'symmetrized_norm.npy')
else:
    runs = np.load(load_address + 'symmetrized.npy')


# Form region of interest from two peaks along the vertical
[peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs[0,:,:] ,2, 
    halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Choose two peaks along the vertical direction to form linescan_y')
[peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs[0,:,:], peak_region, background_region, halfLength=ROI_width/2)
min_y = min([peak_fit_pump[1][1][1], peak_fit_pump[0][1][1]])
max_y = max([peak_fit_pump[1][1][1], peak_fit_pump[0][1][1]])
ROI_x_min = peak_fit_pump[0][0][1]-30
ROI_x_max = peak_fit_pump[0][0][1]+30
ROI_y_min = min_y-10
ROI_y_max = max_y+10
height = ROI_y_max - ROI_y_min
width = ROI_x_max - ROI_x_min

if x_flag ==1:
    # Form region of interest from two peaks along the horizontal
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs[0,:,:] ,2, 
        halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Choose two peaks along the horizontal direction to form linescan_x')
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs[0,:,:], peak_region, background_region, halfLength=ROI_width/2)
    min_x_2 = min([peak_fit_pump[1][0][1], peak_fit_pump[0][0][1]])
    max_x_2 = max([peak_fit_pump[1][0][1], peak_fit_pump[0][0][1]])
    ROI_y_min_2 = peak_fit_pump[0][1][1]-int(30*old_new_conversion)
    ROI_y_max_2 = peak_fit_pump[0][1][1]+int(30*old_new_conversion)
    ROI_x_min_2 = min_x_2-int(10*old_new_conversion)
    ROI_x_max_2 = max_x_2+int(10*old_new_conversion)
    height_2 = ROI_y_max_2 - ROI_y_min_2
    width_2 = ROI_x_max_2 - ROI_x_min_2

# Form difference images
reference = runs[0]
runs_2 =[]
for i in range(1, np.shape(runs)[0]):
    runs_2.append(runs[i]-reference)

# Form linescan of difference image along y
diffuse = []
for i in range(0,np.shape(runs_2)[0]):
    subimage = runs_2[i][ROI_y_min:ROI_y_max,ROI_x_min:ROI_x_max]
    diffuse.append(np.sum(subimage,1))    
if load_address == r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\long_scan\\':
    if different_scans == 2:
        np.save(save_address + 'diffuse_long_expanded_norm.npy', diffuse)
    elif different_scans == 3:
        np.save(save_address + 'diffuse_long_expanded.npy', diffuse)
    elif different_scans == 4:
        np.save(save_address + 'diffuse_long_norm.npy', diffuse)
    elif different_scans == 5:
        np.save(save_address + 'diffuse_long.npy', diffuse)
    else:
        np.save(save_address + 'diffuse_l.npy', diffuse)    
else:
    if different_scans == 2:
        np.save(save_address + 'diffuse_short_expanded_norm.npy', diffuse)
    if different_scans == 3:
        np.save(save_address + 'diffuse_short_expanded.npy', diffuse)
    if different_scans == 4:
        np.save(save_address + 'diffuse_short_norm.npy', diffuse)
    elif different_scans == 5:
        np.save(save_address + 'diffuse_short.npy', diffuse)
    else:
        np.save(save_address + 'diffuse_s.npy', diffuse)   
    
if x_flag == 1:
    # Form linescan of difference image along x
    diffuse_2 = []
    for i in range(0,np.shape(runs_2)[0]):
        subimage = runs_2[i][ROI_y_min_2:ROI_y_max_2,ROI_x_min_2:ROI_x_max_2]
        diffuse_2.append(np.sum(subimage,0)) 
    np.save(save_address + 'diffuse_horizontal.npy', diffuse_2)

# Plot diffuse linescan for each time delay along y
plt.figure()
for i in range(1, np.shape(runs)[0]):
    plt.plot(diffuse[i-1], label = str(round(delay[i],1)) + ' ps')
plt.legend(loc='upper center')
plt.title('Diffuse Linescan Along Vertical')
#plt.xlim(20, 218)
#plt.ylim(-10,325)

if x_flag == 1:
    # Plot diffuse linescan for each time delay along x
    plt.figure()
    for i in range(1, np.shape(runs)[0]):
        plt.plot(diffuse_2[i-1], label = str(int(round(delay[i]))) + ' ps')
    plt.legend(loc='upper center')
    plt.title('Diffuse Linescan Along Horizontal')
    #plt.xlim(20, 218)
    #plt.ylim(-10,325)
    
# Plot the region of interest we are integrating to form the linescan
# Region of interest for vertical
plt.figure()
ax1 = plt.subplot(111)
ax1.imshow(runs[-1]-reference).set_clim(0, 10)
plt.title('Vertical Linescan Region of Interest')
ax1.add_patch(
    patches.Rectangle((ROI_x_min, ROI_y_min), width, height, fill = False, edgecolor = 'white')
    )
    
if x_flag == 1:
    plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(runs[-1]-reference).set_clim(0, 10)
    plt.title('Horizontal Linescan Region of Interest')
    ax1.add_patch(
        patches.Rectangle((ROI_x_min_2, ROI_y_min_2), width_2, height_2, fill = False, edgecolor = 'white')
        )
        
plt.show()