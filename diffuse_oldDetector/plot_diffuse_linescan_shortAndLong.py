'''
Tyler Chase
07/22/2016
This module takes in the diffuse lists from plot_diffuse_linescan.py and plots
both the short and the long scans together on the same graph

Instructions
0. make sure you have already run plot_diffuse_linescan.py
1. Change address of load and save for both long and short scans
2. Change delay stage values for both long and short timescan 
   and make sure time zero is correct
'''







''' Values to change for each run'''
load_address_1 = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
save_address_1 = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
load_address_2 = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\long_scan\\'
save_address_2 = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\long_scan\\'
# Delay stage settings
delayStage_start, delayStage_end, delayStage_step  = [41.4, 45.0, 0.6]
delayStage_start_2, delayStage_end_2, delayStage_step_2  = [37.8, 56.55, 3.75]
time_zero = 42.3












import numpy as np
import matplotlib.pyplot as plt 

def ps_convert(start, stop, step):
    # Delay Stage Settings
    delay_stage = np.arange(start, stop + (step)/2, step)
    speed_of_light = 299792458 # meters/second
    delay = (np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12 # Delay in ps
    return(delay)
    
delay_1 = ps_convert(delayStage_start, delayStage_end, delayStage_step)
delay_2 = ps_convert(delayStage_start_2, delayStage_end_2, delayStage_step_2)

diffuse_1 = np.load(load_address_1 + 'diffuse_vertical.npy')
diffuse_2 = np.load(load_address_2 + 'diffuse_vertical.npy')

# Plot diffuse linescan for each time delay along y
plt.figure()
for i in range(0,np.shape(diffuse_1)[0]):
    plt.plot(diffuse_1[i], label = str(int(round(delay_1[i+1]))) + ' ps')
    
# Plot diffuse linescan for each time delay along y
for i in range(0,np.shape(diffuse_2)[0]):
    plt.plot(diffuse_2[i], label = str(int(round(delay_2[i+1]))) + ' ps')
plt.legend(loc='upper center')
plt.title('Diffuse Linescan Along Vertical')
plt.show()
