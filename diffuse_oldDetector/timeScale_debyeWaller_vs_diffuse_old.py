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
    tau = 10
    A = 0.9
    popt, pcov = opt.curve_fit(exp_debye, x, data, p0 = (tau, A))
    return popt, pcov

# Diffuse intensity Fits no offset
def exp_diffuse(t, tau, A):
    return A * (1-np.exp((-1*t+0.93)/tau))
def exp_diffuse_fit(x, data):
    tau = 10
    A = 100
    popt, pcov = opt.curve_fit(exp_diffuse, x, data, p0 = (tau, A))
    return popt, pcov
    
def line(t, A, B):
    return A*t+B
def line_fit(x, data):
    A = 1
    B = 0
    popt, pcov = opt.curve_fit(line, x, data, p0 = (A, B))
    return popt, pcov
    


# Addresses 
debyeWaller_address_integrated = r"C:\Users\tchase56\Documents\UED\20150121 Au\Normalization_Experimenting\normalize day_2 sum of regions 1 through 3 same intensity version 2\NewDataTo50ps\Data\\"
debyeWaller_address = r"C:\Users\tchase56\Documents\UED\20150121 Au\Normalization_Experimenting\\"
lineScan_short_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\old_data\\'

# Load Data
peak_intensity_integrated = np.loadtxt(debyeWaller_address_integrated + 'Static_Pumped_PeakIntensities_Norm')
peak_intensity = np.loadtxt(debyeWaller_address + 'peak_intensities_timescan')
lineScan_short = np.load(lineScan_short_address + 'diffuse_short_norm.npy')


# Convert delay stage values to delay in picoseconds
times_timescan = np.arange(-2.25,17.75,20.0/50.0)-.5
temp = np.ndarray.tolist(times_timescan)
temp.append(50)
times_timescan = np.array(temp)
times_short = np.array([-2,2,4,6,8,10,50])

# Normalize and scale debye-waller data
peak_intensity_integrated_norm = []
for i in range(5,7):
    peak_intensity_integrated_norm.append(Scale(peak_intensity_integrated[i,:]/peak_intensity_integrated[4,:],1))
peak_intensity_integrated_norm = np.array(peak_intensity_integrated_norm) 
peak_intensity_norm = []
for i in range(6,8):
    peak_intensity_norm.append(Scale(peak_intensity[i,:]/peak_intensity[5,:],5))
    
for i in range(np.shape(peak_intensity_norm)[0]):
    peak_intensity_norm[i] = np.ndarray.tolist(peak_intensity_norm[i])
    peak_intensity_norm[i].append(peak_intensity_integrated_norm[i][-1])
peak_intensity_norm = np.array(peak_intensity_norm)
# Append 50ps point to timescan

# Determine regions for determining changes in diffuse intensity across linescan
length = 49 # needs to be even number with how coded currently
lower_1 = 44
upper_1 = lower_1 + length
upper_5 = 289
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

# Subtract linear fit to background(sub timezero)
lineScan_short_2 = []
popt, pcov = line_fit(np.arange(lower_0, upper_5), lineScan_short[0,lower_0:upper_5])
fit = line(np.arange(np.shape(lineScan_short)[1]), *popt)
'''
plt.figure()
x = np.arange(np.shape(lineScan_short)[1])
plt.plot(line(x,*popt))
plt.plot(lineScan_short[0,:])
plt.axvline(lower_0, linestyle = ':', c = 'k')
plt.axvline(upper_5, linestyle = ':', c = 'k')
plt.show()
'''
for i in range(np.shape(np.array(lineScan_short))[0]):
    lineScan_short_2.append(lineScan_short[i,:]-fit)
lineScan_short = np.array(lineScan_short_2)
    
    

# Sum over regions determined in previous step
sums = []
for i in range(np.shape(np.array(lineScan_short))[0]):
    sums.append([])
    sums[i].append(np.sum(lineScan_short[i,lower_0:upper_0])/float(upper_0-lower_0))
    sums[i].append(np.sum(lineScan_short[i,lower_1:upper_1])/float(upper_1-lower_1))
    sums[i].append(np.sum(lineScan_short[i,lower_2:upper_2])/float(upper_2-lower_2))
    sums[i].append(np.sum(lineScan_short[i,lower_3:upper_3])/float(upper_3-lower_3))
    sums[i].append(np.sum(lineScan_short[i,lower_4:upper_4])/float(upper_4-lower_4))
    sums[i].append(np.sum(lineScan_short[i,lower_5:upper_5])/float(upper_5-lower_5))
sums = np.array(sums)

# Fit the debye-waller data
fits_debye = []
error_debye = []
debye_region_bottom = 7
debye_region_top = np.shape(peak_intensity_norm)[1]
#debye_region_top = -1
for i in range(np.shape(peak_intensity_norm)[0]):
    fits_debye.append(exp_debye_fit(times_timescan[debye_region_bottom:debye_region_top], peak_intensity_norm[i,debye_region_bottom:debye_region_top])[0])
    error_debye.append(exp_debye_fit(times_timescan[debye_region_bottom:debye_region_top], peak_intensity_norm[i,debye_region_bottom:debye_region_top])[1][0,0])

# Fit the diffuse data
fits_diffuse = []
error_diffuse = []
diffuse_region_bottom = 1
diffuse_region_top = np.shape(sums)[0]
#diffuse_region_top = -1
for i in range(np.shape(sums)[1]):
    fits_diffuse.append(exp_diffuse_fit(times_short[diffuse_region_bottom:diffuse_region_top], sums[diffuse_region_bottom:diffuse_region_top,i])[0])
    error_diffuse.append(np.sqrt(exp_diffuse_fit(times_short[diffuse_region_bottom:diffuse_region_top], sums[diffuse_region_bottom:diffuse_region_top,i])[1][0,0]))



# Plot figure similar to paper
f, (ax1, ax2) = plt.subplots(2, 1, sharex = True)

# Subplot with debye-waller
colors_0 = ['k', 'b', 'r', 'g', 'c']
for i in range(np.shape(peak_intensity_norm)[0]):
    ax1.scatter(times_timescan, peak_intensity_norm[i], c= colors_0[i])
    #ax1.scatter(times_short, peak_intensity_integrated_norm[i], c = colors_0[i])
    
# Plot vertical lines showing fit region
debye_region_bottom_ps = times_timescan[debye_region_bottom]
debye_region_top_ps = times_timescan[debye_region_top-1]
ax1.axvline(debye_region_bottom_ps, linestyle = ':', c = 'k')
ax1.axvline(debye_region_top_ps, linestyle = ':', c = 'k')
ax1.set_xlim(-5,55)
ax1.set_ylim(0.9, 1.01)

# debye-waller fits
for i in range(np.shape(peak_intensity_norm)[0]):
    ax1.plot(np.arange(debye_region_bottom_ps, debye_region_top_ps, 0.5), exp_debye(np.arange(debye_region_bottom_ps, debye_region_top_ps, 0.5), *fits_debye[i]), c = colors_0[i])
    
# Subplot with diffuse regions
colors = ['y', 'm', 'c', 'r', 'g', 'b']
for i in range(np.shape(sums)[1]):
    ax2.scatter(times_short, sums[:,i], label = str(i), c = colors[i])
    
# Plot vertical line showing fit region
diffuse_region_bottom_ps = times_short[diffuse_region_bottom]
diffuse_region_top_ps = times_short[diffuse_region_top-1]
ax2.axvline(diffuse_region_bottom_ps, linestyle = ':', c = 'k')
ax2.axvline(diffuse_region_top_ps, linestyle = ':', c = 'k')

# diffuse fits
for i in range(np.shape(sums)[1]):
    ax2.plot(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.5), exp_diffuse(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.5), *fits_diffuse[i]), c = colors[i])
ax2.set_ylim(-25,300)
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
    ax1.plot(lineScan_short[i], label = str(round(times_short[i])))    
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
ax1.legend(loc = 'upper center')
ax2 = ax1.twinx()
ax2.bar(lower_0, fits_diffuse[0][0], upper_0-lower_0, color = colors[0], yerr = error_diffuse[0], alpha = 0.30, ecolor = 'k')
ax2.bar(lower_1, fits_diffuse[1][0], upper_1-lower_1, color = colors[1], yerr = error_diffuse[1], alpha = 0.20, ecolor = 'k')
ax2.bar(lower_2, fits_diffuse[2][0], upper_2-lower_2, color = colors[2], yerr = error_diffuse[2], alpha = 0.20, ecolor = 'k')
ax2.bar(lower_3, fits_diffuse[3][0], upper_3-lower_3, color = colors[3], yerr = error_diffuse[3], alpha = 0.20, ecolor = 'k')
ax2.bar(lower_4, fits_diffuse[4][0], upper_4-lower_4, color = colors[4], yerr = error_diffuse[4], alpha = 0.20, ecolor = 'k')
ax2.bar(lower_5, fits_diffuse[5][0], upper_5-lower_5, color = colors[5], yerr = error_diffuse[5], alpha = 0.20, ecolor = 'k')
ax2.set_xlim(20,312)
ax2.set_ylim(0,7)
plt.show()
