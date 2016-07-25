# OLD DETECTOR
'''
Tyler Chase
06/27/2016
This module reads in the output files from align_average.py or align_different_scans.py
and subtracts a before time zero reference to see diffuse scattering

Instructions
0. Make sure you have already run align_average.py or align_different_scans.py
   (set different_scans flag accordingly)
1. Change address of output file you want to load from step 0 
2. Change delay stage values and make sure time zero is correct
3. Decide if you want to oversaturate the image and what you want the maximum to be
   (Do this with the clim_flag and clim_value)
4. Which run in the run list do you want as the reference? (reference_num)
'''

''' Values to change for each run'''
load_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
# Delay stage settings
delayStage_start, delayStage_end, delayStage_step  = [41.4, 45.0, 0.6]
time_zero = 42.3
# Flags
clim_flag = 1     # Do you want to saturate large values on colorbar when plotting?
clim_value = 10     # maximum value on colorbar when saturating large intensities
# Which image to use for reference?
reference_num = 0
# Flag for different file outputs
# (0, alignment output)
# (1, averaging runs output)
# (2, centered rotated)
# (3, symmetrized)
# (4, centered rotated expanded)
# (5, centered rotated expanded symmetrized)
different_scans = 5










import numpy as np
import matplotlib.pyplot as plt 

# Delay Stage Settings
delay_stage = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
speed_of_light = 299792458 # meters/second
delay = (np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12 # Delay in ps

# Load data that you would like to see diffuse scattering in
if different_scans == 0:
    pumped = np.load(load_address + 'averaged_aligned.npy')
    reference = pumped[reference_num,:,:]
elif different_scans == 1:
    pumped = np.load(load_address + 'averaged_runs.npy')
    reference = pumped[reference_num,:,:]
elif different_scans == 2:
    pumped = np.load(load_address + 'centered_rotated.npy')
    reference = pumped[reference_num,:,:]
elif different_scans == 3:
    pumped = np.load(load_address + 'symmetrized.npy')
    reference = pumped[reference_num,:,:]
elif different_scans == 4:
    pumped = np.load(load_address + 'expanded.npy')
    reference = pumped[reference_num,:,:]
else:
    pumped = np.load(load_address + 'expanded_symmetrized.npy')
    reference = pumped[reference_num,:,:]

# Plot diffuse scattering
if clim_flag == 1:
    for i in range(0,np.shape(pumped)[0]):
        plt.figure()
        plt.imshow(pumped[i]-reference, interpolation='none').set_clim(0,clim_value)
        plt.colorbar()
        plt.title('Time Delay {:g} ps'.format(delay[i]))
else:
    for i in range(0,np.shape(pumped)[0]):
        plt.figure()
        plt.imshow(pumped[i]-reference, interpolation='none')
        plt.colorbar()
        plt.title('Time Delay {:g} ps'.format(delay[i]))
