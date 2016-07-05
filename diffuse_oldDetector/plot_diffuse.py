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
3. 
'''

''' Values to change for each run'''
load_address = r'E:\Klaus\20160630\\'
# Delay stage settings
#delayStage_start, delayStage_end, delayStage_step  = [64.3, 65.5, 0.2] 
delayStage_start, delayStage_end, delayStage_step  = [64.3, 65.7, 0.2]

time_zero = 64.47
# Flags
clim_flag = 1     # Do you want to saturate large values on colorbar when plotting?
clim_value = 50     # maximum value on colorbar when saturating large intensities
# Which image to use for reference?
reference_num = 0
# Is this the output of align_average or the output of align_different_scans
different_scans = 2










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
    pumped = np.load(load_address + 'averaged_runs_highFluence.npy')
    reference = pumped[reference_num,:,:]
else:
    pumped = np.load(load_address + 'symmetrized_highFluence.npy')
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
