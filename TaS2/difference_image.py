'''
Tyler Chase
06/14/2016
average multiple time delays together and form difference images to see changes
after the pump

Instructions
0. Make sure you have already run align_average.py or align_different_scans.py or symmetrize
   (set different_scans flag accordingly)
1. Change addresses
2. Change delay stage values and make sure time zero is correct
3. Choose how many time delays you would like averaged together with group variable 
   (only intended to improve statistics)
4. Run the code
5. May need to change the set_clim for contrast line 77 or line 70
'''

import numpy as np
import matplotlib.pyplot as plt 





''' Values to change for each run'''
load_address = r'C:\Users\tchase56\Documents\UED\TaS2\20160611\Data\timescan\\'
save_address = r'C:\Users\tchase56\Documents\UED\TaS2\20160611\Data\timescan\\'
# Is this the output of align_average, the output of align_different_scans, or the output of symmetrize (0, 1, or 2)
different_scans = 2
# Delay stage settings
time_zero = 63.2
delayStage_start, delayStage_end, delayStage_step  = [63.035, 63.485, 0.0075] 
# How many time delays would you like to average together to improve statistics?
group = 10










# Determine delay in ps from delay stage
delay_stage = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
speed_of_light = 299792458 # meters/second
delay = np.round((np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12,2) # Delay in ps

# Read in output of align_average or align_different_scans and group
if different_scans == 0:
    images = np.load(load_address + 'averaged_aligned.npy')
elif different_scans == 1:
    images = np.load(load_address + 'averaged_runs.npy')
else:
    images = np.load(load_address + 'symmetrized.npy')
grouped = []
times = []
for i in range(0, np.shape(images)[0]//group):
    grouped.append(np.average(images[group*(i):group*(i+1)], axis = 0))
    times.append(np.average(delay[group*(i):group*(i+1)]))

# Take the difference between the time_delay(t) and time_delay(0) of grouped
difference = []
for i in range(0, np.shape(grouped)[0]):
    difference.append(grouped[i]-(grouped[0]))

# Plot pre time zero diffraction pattern    
plt.figure()
plt.imshow(grouped[0]).set_clim(300,1750)
plt.colorbar()
plt.show()
    
# Plot differences
for i in range(0, np.shape(difference)[0]):
    plt.figure()
    plt.imshow(difference[i], cmap = 'coolwarm').set_clim(-35,35)
    plt.title(str(times[i]) + ' PS')
    plt.colorbar()   
plt.show()