# Old Gold Paper Timescales

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from scipy.interpolate import interp1d 

def ps_convert(start, stop, step, time_zero):
    # Delay Stage Settings
    delay_stage = np.arange(start, stop + (step)/2, step)
    speed_of_light = 299792458 # meters/second
    delay = (np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12 # Delay in ps
    return delay
    
def combine_short_long(short, long):
    combined = []
    combined.append(long[0])
    for i in range(len(short)):
        combined.append(short[i])
    for i in range(1,len(long)):
        combined.append(long[i])
    return combined
    
# Used to normalize debye-waller time scans with multiple values below time zero
def Scale(graph, points):
    points = list(range(points))
    norm = np.average(np.array(graph[points]))
    graph = graph / norm
    return graph

def exp_debye(t, tau, A):
    return (1-A)*np.exp((-1*t+0)/tau)+A  
def exp_debye_fit(x, data):
    tau = 5
    A = 0.9
    popt, pcov = opt.curve_fit(exp_debye, x, data, p0 = (tau, A))
    return popt, pcov

# Diffuse intensity Fits no offset
def exp_diffuse(t, tau, A):
    return A * (1-np.exp((-1*t+0)/tau))
def exp_diffuse_fit(x, data):
    tau = 5
    A = 300
    popt, pcov = opt.curve_fit(exp_diffuse, x, data, p0 = (tau, A,))
    return popt, pcov
    


# Addresses 
debyeWaller_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\timescan\scan32\images-ANDOR1\\'
lineScan_short_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
lineScan_long_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\long_scan\\'

# Load Data
peak_intensity = np.load(debyeWaller_address + 'peak_intensity_gaussian.npy')
lineScan_short = np.load(lineScan_short_address + 'diffuse_short_expanded_norm.npy')
lineScan_long = np.load(lineScan_long_address + 'diffuse_long_expanded_norm.npy')

# Note delay stage values and timezero
time_zero = 42.45
delayStage_start_timescan, delayStage_end_timescan, delayStage_step_timescan  = [39.3, 57.3, 0.15]
delayStage_start_short, delayStage_end_short, delayStage_step_short  = [41.4, 45.0, 0.6]
delayStage_start_long, delayStage_end_long, delayStage_step_long  = [37.8, 56.55, 3.75]

# Convert delay stage values to delay in picoseconds
times_timescan = ps_convert(delayStage_start_timescan, delayStage_end_timescan, delayStage_step_timescan, time_zero)
times_short = ps_convert(delayStage_start_short, delayStage_end_short, delayStage_step_short, time_zero)
times_long = ps_convert(delayStage_start_long, delayStage_end_long, delayStage_step_long, time_zero)
times_combined = combine_short_long(times_short[1:], times_long[1:])

# Normalize and scale debye-waller data
peak_intensity_norm = []
for i in range(1,np.shape(peak_intensity)[0]):
    peak_intensity_norm.append(Scale(peak_intensity[i,:]/peak_intensity[0,:],18))
peak_intensity_norm = np.array(peak_intensity_norm)
    


# Determine regions for determining changes in diffuse intensity across linescan
length = 35 # needs to be even number with how coded currently
lower_1 = 29
upper_1 = lower_1 + length
upper_5 = 204
lower_5 = upper_5 - length
lower_3 = int((upper_1+lower_5)/2-(length/2))
upper_3 = int((upper_1+lower_5)/2+(length/2))
lower_2 = int((upper_1+lower_3)/2-(length/2))
upper_2 = int((upper_1+lower_3)/2+(length/2))
lower_4 = int((lower_5+upper_3)/2-(length/2))
upper_4 = int((lower_5+upper_3)/2+(length/2))
lower_0 = lower_1
upper_0 = int(lower_0 + (length/2))
lower_1 = int(upper_1 - (length/2))

# Sum over regions determined in previous step
sums_short = []
for i in range(np.shape(np.array(lineScan_short))[0]):
    sums_short.append([])
    sums_short[i].append(np.sum(lineScan_short[i,lower_0:upper_0])/float(upper_0-lower_0))
    sums_short[i].append(np.sum(lineScan_short[i,lower_1:upper_1])/float(upper_1-lower_1))
    sums_short[i].append(np.sum(lineScan_short[i,lower_2:upper_2])/float(upper_2-lower_2))
    sums_short[i].append(np.sum(lineScan_short[i,lower_3:upper_3])/float(upper_3-lower_3))
    sums_short[i].append(np.sum(lineScan_short[i,lower_4:upper_4])/float(upper_4-lower_4))
    sums_short[i].append(np.sum(lineScan_short[i,lower_5:upper_5])/float(upper_5-lower_5))
offset = np.shape(np.array(lineScan_short))[0] - 1
sums_long = []
for i in range(np.shape(np.array(lineScan_long))[0]):
    sums_long.append([])
    sums_long[i].append(np.sum(lineScan_long[i,lower_0:upper_0])/float(upper_0-lower_0))
    sums_long[i].append(np.sum(lineScan_long[i,lower_1:upper_1])/float(upper_1-lower_1))
    sums_long[i].append(np.sum(lineScan_long[i,lower_2:upper_2])/float(upper_2-lower_2))
    sums_long[i].append(np.sum(lineScan_long[i,lower_3:upper_3])/float(upper_3-lower_3))
    sums_long[i].append(np.sum(lineScan_long[i,lower_4:upper_4])/float(upper_4-lower_4))
    sums_long[i].append(np.sum(lineScan_long[i,lower_5:upper_5])/float(upper_5-lower_5))
sums = combine_short_long(sums_short, sums_long)
sums = np.array(sums)

'''
# Normalize the end of these sums to 1 so you can see the time constants easier
for i in range(np.shape(sums)[1]):
    sums[:,i] = sums[:,i]/sums[-1,i]
'''

# Fit the debye-waller data
fits_debye = []
error_debye = []
debye_region_bottom = 21
#debye_region_top = -1
debye_region_top = np.shape(peak_intensity_norm)[1]
for i in range(np.shape(peak_intensity_norm)[0]):
    fits_debye.append(exp_debye_fit(times_timescan[debye_region_bottom:debye_region_top], peak_intensity_norm[i,debye_region_bottom:debye_region_top])[0])
    error_debye.append(np.sqrt(exp_debye_fit(times_timescan[debye_region_bottom:debye_region_top], peak_intensity_norm[i,debye_region_bottom:debye_region_top])[1][0,0]))

# Fit the diffuse data
fits_diffuse = []
error_diffuse = []
diffuse_region_bottom = 2
#diffuse_region_top = -1
diffuse_region_top = np.shape(sums)[0]
for i in range(np.shape(sums)[1]):
    fits_diffuse.append(exp_diffuse_fit(times_combined[diffuse_region_bottom:diffuse_region_top], sums[diffuse_region_bottom:diffuse_region_top,i])[0])
    error_diffuse.append(np.sqrt(exp_diffuse_fit(times_combined[diffuse_region_bottom:diffuse_region_top], sums[diffuse_region_bottom:diffuse_region_top,i])[1][0,0]))



# Plot figure similar to paper
f, (ax1, ax2) = plt.subplots(2, 1, sharex = True)

# Subplot with debye-waller
colors_0 = ['k', 'b', 'r', 'g', 'c']
for i in range(np.shape(peak_intensity_norm)[0]):
    ax1.scatter(times_timescan, peak_intensity_norm[i], c= colors_0[i])
    
# Plot vertical lines showing fit region
debye_region_bottom_ps = times_timescan[debye_region_bottom]
debye_region_top_ps = times_timescan[debye_region_top-1]
ax1.axvline(debye_region_bottom_ps, linestyle = ':', c = 'k')
ax1.axvline(debye_region_top_ps, linestyle = ':', c = 'k')
ax1.set_xlim(-23,101)

# debye-waller fits
for i in range(np.shape(peak_intensity_norm)[0]):
    ax1.plot(np.arange(debye_region_bottom_ps, debye_region_top_ps, 0.5), exp_debye(np.arange(debye_region_bottom_ps, debye_region_top_ps, 0.5), *fits_debye[i]), c = colors_0[i])
    
# Subplot with diffuse regions
colors = ['y', 'm', 'c', 'r', 'g', 'b']
for i in range(np.shape(sums)[1]):
    ax2.scatter(times_combined, sums[:,i], label = str(i), c = colors[i])
    
# Plot vertical line showing fit region
diffuse_region_bottom_ps = times_combined[diffuse_region_bottom]
diffuse_region_top_ps = times_combined[diffuse_region_top-1]
ax2.axvline(diffuse_region_bottom_ps, linestyle = ':', c = 'k')
ax2.axvline(diffuse_region_top_ps, linestyle = ':', c = 'k')
ax2.set_ylim(-25,375)


# diffuse fits
for i in range(np.shape(sums)[1]):
    ax2.plot(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.5), exp_diffuse(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.5), *fits_diffuse[i]), c = colors[i])
plt.show()







'''
plt.figure()
for i in range(np.shape(sums)[1]):
    plt.plot(times_combined, sums[:,i], label = str(i))
plt.legend(loc = 'lower right')
plt.show()
'''    

plt.figure()
ax1 = plt.subplot(111)
for i in range(np.shape(np.array(lineScan_short))[0]):
    ax1.plot(lineScan_short[i], label = str(int(round(times_short[i+1]))) + ' ps')    
ax1.axvline(lower_0, linestyle = '--', c = 'k')
ax1.axvline(upper_0, linestyle = '--', c = 'k')
ax1.axvline(lower_1, linestyle = '--', c = 'k')
ax1.axvline(upper_1, linestyle = '--', c = 'k')
ax1.axvline(lower_2, linestyle = '--', c = 'k')
ax1.axvline(upper_2, linestyle = '--', c = 'k')
ax1.axvline(lower_3, linestyle = '--', c = 'k')
ax1.axvline(upper_3, linestyle = '--', c = 'k')
ax1.axvline(lower_4, linestyle = '--', c = 'k')
ax1.axvline(upper_4, linestyle = '--', c = 'k')
ax1.axvline(lower_5, linestyle = '--', c = 'k')
ax1.axvline(upper_5, linestyle = '--', c = 'k')
ax1.set_ylim(-50,600)

for i in [2, 4]:
#for i in range(1,np.shape(np.array(lineScan_long))[0]):
    plt.plot(lineScan_long[i], label = str(int(round(times_long[i+1]))) + ' ps')
ax1.legend(loc = 'upper center')

ax2 = ax1.twinx()
ax2.bar(lower_0, fits_diffuse[0][0], upper_0-lower_0, color = colors[0], yerr = error_diffuse[0], alpha = 0.30, ecolor = 'k')
ax2.bar(lower_1, fits_diffuse[1][0], upper_1-lower_1, color = colors[1], yerr = error_diffuse[1], alpha = 0.20, ecolor = 'k')
ax2.bar(lower_2, fits_diffuse[2][0], upper_2-lower_2, color = colors[2], yerr = error_diffuse[2], alpha = 0.20, ecolor = 'k')
ax2.bar(lower_3, fits_diffuse[3][0], upper_3-lower_3, color = colors[3], yerr = error_diffuse[3], alpha = 0.20, ecolor = 'k')
ax2.bar(lower_4, fits_diffuse[4][0], upper_4-lower_4, color = colors[4], yerr = error_diffuse[4], alpha = 0.20, ecolor = 'k')
ax2.bar(lower_5, fits_diffuse[5][0], upper_5-lower_5, color = colors[5], yerr = error_diffuse[5], alpha = 0.20, ecolor = 'k')
ax2.set_xlim(16,215)
ax2.set_ylim(0,7)
plt.show()
