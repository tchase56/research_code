'''
Tyler Chase
06/14/2016
average multiple time delays together and form difference images to see changes
after the pump
'''


import numpy as np
import matplotlib.pyplot as plt 





''' Values to change for each run'''
load_address = r'C:\Users\tchase56\Documents\UED\TaS2\20160611\Data\timescan\\'
save_address = r'C:\Users\tchase56\Documents\UED\TaS2\20160611\Data\timescan\\'
# Delay stage settings
time_zero = 63.2
delayStage_start, delayStage_end, delayStage_step  = [63.035, 63.485, 0.0075] 
# How many time delays would you like to average together to improve statistics?
group = 10
















# Determine delay in ps from delay stage
delay_stage = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
speed_of_light = 299792458 # meters/second
delay = np.round((np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12,2) # Delay in ps


# Load Mask and Images
'''
boundary_mask = np.load(load_address + 'boundary_mask.npy')
peak_mask = np.load(load_address + 'peak_mask.npy')
peak_mask_inv = -1.0*(peak_mask-1.0)
throughbeam_mask = np.load(load_address + 'throughbeam_mask.npy')
background_mask = np.load(load_address + 'mask_detector_background.npy')
'''
images = np.load(load_address + 'averaged_runs.npy')

grouped = []
times = []
for i in range(0, np.shape(images)[0]//group):
    grouped.append(np.average(images[group*(i):group*(i+1)], axis = 0))
    times.append(np.average(delay[group*(i):group*(i+1)]))

difference = []
for i in range(0, np.shape(grouped)[0]):
    difference.append(grouped[i]-(grouped[0]))
    
plt.figure()
plt.imshow(grouped[0]).set_clim(300,1750)
plt.colorbar()
plt.show()
    
for i in range(0, np.shape(difference)[0]):
    plt.figure()
    plt.imshow(difference[i], cmap = 'coolwarm').set_clim(-35,35)
    plt.title(str(times[i]) + ' PS')
    plt.colorbar()

'''
for i in range(0, np.shape(grouped)[0]):
    plt.figure()
    plt.imshow(grouped[i])
    plt.title(str(times[i]))
'''
    
plt.show()