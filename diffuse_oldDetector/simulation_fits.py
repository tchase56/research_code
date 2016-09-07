import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from scipy.interpolate import interp1d 
import scipy.io
import Region_Of_Interest_2
import Fit_Peaks_2
import matplotlib.patches as patches

# Flags
peaks_saved = 1
plot_simulation = 0
plot_data = 0
plot_simulation_ROIs = 0
plot_data_ROIs = 1
plot_short_region1 = 0
plot_short_region2 = 0
plot_long_region1 = 0
plot_long_region2 = 0
plot_fitParams = 0
# Region widths
ROI_width = 50
halfPipeLength_S = 11
# Delay stage values for timescans
delayStage_start_short, delayStage_end_short, delayStage_step_short  = [41.4, 45.0, 0.6]
delayStage_start_long, delayStage_end_long, delayStage_step_long  = [37.8, 56.55, 3.75]
#time_zero = 42.3
time_zero = 42.45
# Graph clim values
simulation_clim = 4e10
data_clim = 10





    


# Longitudinal Linescan From Theory
def theory_L(x, A, B):
    return A*(spline_T1_2(x)+spline_T2_2(x)+spline_L_2(x))+B     
# Longitudinal Linescan From Theory
def theory_L_fit(x, data):
    fit_func_L = lambda z, A, B: theory_L(z, A, B) 
    A = 2e-9
    B = 100
    popt, pcov = opt.curve_fit(fit_func_L, x, data, p0=(A, B))
    return popt    
    
# Transverse Linescan From Theory
def theory_T(x, A, B):
    return A*(spline_T1_1(x)+spline_T2_1(x)+spline_L_1(x))+B    
# Transverse Linescan From Theory
def theory_T_fit(x, data):
    fit_func_T = lambda z, A, B: theory_T(z, A, B)
    A = 2e-9
    B = 100
    popt, pcov = opt.curve_fit(fit_func_T, x, data, p0=(A, B))
    return popt

# Take in ROI information as well as phonon branch and return the phonon linescan value    
def ROI(xmin, xmax, ymin, ymax, branch, axis = 0, plot = 0, clim = 0, lim_start = 0):
    xmin = int(round(xmin))
    xmax = int(round(xmax))
    ymin = int(round(ymin))
    ymax = int(round(ymax))
    lim_start = int(round(lim_start))
    subimage = branch[xmin:xmax,ymin:ymax]
    linescan = np.sum(subimage, axis = axis)
    if plot == 1:
        plt.figure()
        ax1 = plt.subplot(111)
        if clim == 0:
            ax1.imshow(branch)
        else:
            ax1.imshow(branch).set_clim(0,clim)
        ax1.add_patch(
            patches.Rectangle((xmin, ymin), (xmax-xmin+1), (ymax-ymin+1), fill = False, edgecolor = 'white')
            )
        plt.xlim(lim_start, np.shape(branch)[0]-lim_start-1)
        plt.ylim(lim_start, np.shape(branch)[0]-lim_start-1)
        plt.show()
    return linescan
    
def ps_convert(start, stop, step):
    # Delay Stage Settings
    delay_stage = np.arange(start, stop + (step)/2, step)
    speed_of_light = 299792458 # meters/second
    delay = (np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12 # Delay in ps
    return(delay)
    
def combine_short_long(short, long):
    combined = []
    combined.append(long[0])
    for i in short:
        combined.append(i)
    for i in long[1:]:
        combined.append(i)
    return(combined)
    
def fancy_plot(t, A, B, title):
    fig, ax1 = plt.subplots()
    ax1.plot(t, A, 'bo-')
    ax1.set_xlabel('Time Delay(ps)')
    ax1.set_ylabel ('Fit Parameter A', color = 'b')
    for i in ax1.get_yticklabels():
        i.set_color('b')
    ax2= ax1.twinx()
    ax2.plot(t, B, 'ro-')
    ax2.set_ylabel('Fit Parameter B', color = 'r')
    for i in ax2.get_yticklabels():
        i.set_color('r')
    plt.title(title)
    plt.show
    
# Physical Constants
hbar = 1.05457173e-34
k = 1.3806488e-23
T_1 = 300
A = 5e-12
y_lim = 650

# Load Mariano's Simulation Data
simulation_path = r'C:\Users\tchase56\Documents\UED\20150121 Au\Normalization_Experimenting\q_exchange_for_pixels\\'
simulation = scipy.io.loadmat(simulation_path + 'Au_Sep17b')
simulation_old = scipy.io.loadmat(simulation_path + 'Au_Sep1') 
simulation_T2 = simulation['intens'][0,:,:]
simulation_T1 = simulation['intens'][1,:,:]
simulation_L = simulation['intens'][2,:,:]
simulation_total = simulation['intens'][0,:,:] + simulation['intens'][1,:,:] + simulation['intens'][2,:,:]

# Name Peaks Simulation
peak_S = []
peak_S.append([])
peak_S[0].append([86, 86])
peak_S[0].append([169, 86])
peak_S[0].append([169, 169])
peak_S[0].append([86, 169])
peak_S.append([])
peak_S[1].append([127, 45])
peak_S[1].append([210, 127])
peak_S[1].append([127, 210])
peak_S[1].append([45, 127])
peak_S.append([])
peak_S[2].append([45, 45])
peak_S[2].append([210, 45])
peak_S[2].append([210, 210])
peak_S[2].append([45, 210])
x_center = np.average(peak_S[0], axis = 0)

# Load Experimental Data
data_path_short = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\short_scan\\'
data_path_long = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\long_scan\\'
save_address = r'C:\Users\tchase56\Documents\UED\Gold_Diffuse_03212016\20160512\\'
data_short = np.array(np.load(data_path_short + 'expanded_symmetrized.npy'))
data_long = np.array(np.load(data_path_long + 'expanded_symmetrized.npy'))

# Name Peaks Data
if peaks_saved == 0:
    peak_D = []
    for i in range(3):
        peak_D.append([])
        plt.figure()
        plt.imshow(simulation_total).set_clim(0, 3e10)
        color = ['b', 'g', 'r', 'k']
        for j in range(4):
            plt.scatter(peak_S[i][j][0], peak_S[i][j][1], c = color[j] )
        [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(data_long[-1,:,:] ,4, 
            halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Choose four peaks in order of dots in figure one (blue, green, red, black)')
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(data_long[-1,:,:], peak_region, background_region, halfLength=ROI_width/2)
        #plt.figure()
        #plt.imshow(data_long[-1,:,:])
        for j in range(4):        
            peak_D[i].append([peak_fit_pump[j][0][1], peak_fit_pump[j][1][1]])
            #plt.scatter(peak_D[i][j][0], peak_D[i][j][1], c = color[j])
    np.save(save_address + 'data_peaks.npy', peak_D)
else:
    peak_D = np.load(save_address + 'data_peaks.npy')
    
# Determing conversion factor between simulation and data
length_D = peak_D[1][1][0] - peak_D[2][0][0] 
length_S = peak_S[1][1][0] - peak_S[2][0][0] 
factor = length_D/length_S
halfPipeLength_D = int(round(halfPipeLength_S * factor))
# Determine cutoff values for ROI plots
lim_start_S = peak_S[2][0][0]-20
lim_start_D = peak_D[2][0][0]-int(round(factor*20))

if plot_simulation == 1:
    # Plot simulation of phonon modes
    plt.figure()
    plt.imshow(simulation_T1).set_clim(0,simulation_clim)
    plt.title('Transverse 1')
    plt.xlim(lim_start_S, np.shape(simulation_total)[0]-lim_start_S-1)
    plt.ylim(lim_start_S, np.shape(simulation_total)[0]-lim_start_S-1)
    plt.figure()
    plt.imshow(simulation_T2).set_clim(0,simulation_clim)
    plt.title('Transverse 2')
    plt.xlim(lim_start_S, np.shape(simulation_total)[0]-lim_start_S-1)
    plt.ylim(lim_start_S, np.shape(simulation_total)[0]-lim_start_S-1)
    plt.figure()
    plt.imshow(simulation_L).set_clim(0,simulation_clim)
    plt.title('Longitudinal')
    plt.xlim(lim_start_S, np.shape(simulation_total)[0]-lim_start_S-1)
    plt.ylim(lim_start_S, np.shape(simulation_total)[0]-lim_start_S-1)
    plt.figure()
    plt.imshow(simulation_total).set_clim(0,simulation_clim)
    plt.title('All Three Branches')
    plt.xlim(lim_start_S, np.shape(simulation_total)[0]-lim_start_S-1)
    plt.ylim(lim_start_S, np.shape(simulation_total)[0]-lim_start_S-1)
    plt.show()

# Determine the time delays for the long and short scans
delay_short = ps_convert(delayStage_start_short, delayStage_end_short, delayStage_step_short)
delay_long = ps_convert(delayStage_start_long, delayStage_end_long, delayStage_step_long)


if plot_data == 1:
    # Plot short scan diffuse scattering
    for i in range(1,np.shape(delay_short)[0]):
        plt.figure()
        plt.imshow(data_short[i,:,:]-data_short[0,:,:]).set_clim(0,data_clim)
        plt.title(str(int(round(delay_short[i]))) + 'PS')
        plt.xlim(lim_start_D, np.shape(data_long)[1]-lim_start_D-1)
        plt.ylim(lim_start_D, np.shape(data_long)[1]-lim_start_D-1)
        plt.show()
    # Plot long scan diffuse scattering
    for i in range(2, np.shape(delay_long)[0]):
        plt.figure()
        plt.imshow(data_long[i,:,:]-data_long[0,:,:]).set_clim(0,data_clim)
        plt.title(str(int(round(delay_long[i]))) + 'PS')
        plt.xlim(lim_start_D, np.shape(data_long)[1]-lim_start_D-1)
        plt.ylim(lim_start_D, np.shape(data_long)[1]-lim_start_D-1)
        plt.show()


    
# Determine Linescans Simulation
line_T1_1 = ROI(peak_S[2][1][0]-halfPipeLength_S, peak_S[2][1][0]+halfPipeLength_S, peak_S[2][1][1], peak_S[1][1][1], simulation_T1)
line_T2_1 = ROI(peak_S[2][1][0]-halfPipeLength_S, peak_S[2][1][0]+halfPipeLength_S, peak_S[2][1][1], peak_S[1][1][1], simulation_T2)
line_L_1 = ROI(peak_S[2][1][0]-halfPipeLength_S, peak_S[2][1][0]+halfPipeLength_S, peak_S[2][1][1], peak_S[1][1][1], simulation_L)
line_T1_2 = ROI(peak_S[1][2][0]-halfPipeLength_S, peak_S[1][2][0]+halfPipeLength_S, peak_S[0][2][1], peak_S[1][2][1], simulation_T1)
line_T2_2 = ROI(peak_S[1][2][0]-halfPipeLength_S, peak_S[1][2][0]+halfPipeLength_S, peak_S[0][2][1], peak_S[1][2][1], simulation_T1)
line_L_2 = ROI(peak_S[1][2][0]-halfPipeLength_S, peak_S[1][2][0]+halfPipeLength_S, peak_S[0][2][1], peak_S[1][2][1], simulation_T1)
# Plot simulation
if plot_simulation_ROIs == 1:
    ROI(peak_S[2][1][0]-halfPipeLength_S, peak_S[2][1][0]+halfPipeLength_S, peak_S[2][1][1], peak_S[1][1][1], simulation_total, plot = 1, clim = simulation_clim, lim_start = lim_start_S)
    ROI(peak_S[1][2][0]-halfPipeLength_S, peak_S[1][2][0]+halfPipeLength_S, peak_S[0][2][1], peak_S[1][2][1], simulation_total, plot = 1, clim = simulation_clim, lim_start = lim_start_S)

# Determine Linescans Data
line_D_1_short = []
line_D_2_short = []
for i in range(1,np.shape(data_short)[0]):
    line_D_1_short.append(ROI(peak_D[2][1][0]-halfPipeLength_D, peak_D[2][1][0]+halfPipeLength_D, peak_D[2][1][1], peak_D[1][1][1], data_short[i]-data_short[0]))
    line_D_2_short.append(ROI(peak_D[1][2][0]-halfPipeLength_D, peak_D[1][2][0]+halfPipeLength_D, peak_D[0][2][1], peak_D[1][2][1], data_short[i]-data_short[0]))
line_D_1_long = []
line_D_2_long = []
for i in range(1,np.shape(data_long)[0]):
    line_D_1_long.append(ROI(peak_D[2][1][0]-halfPipeLength_D, peak_D[2][1][0]+halfPipeLength_D, peak_D[2][1][1], peak_D[1][1][1], data_long[i]-data_long[0]))
    line_D_2_long.append(ROI(peak_D[1][2][0]-halfPipeLength_D, peak_D[1][2][0]+halfPipeLength_D, peak_D[0][2][1], peak_D[1][2][1], data_long[i]-data_long[0]))  
# Plot data
if plot_data_ROIs == 1:
    ROI(peak_D[2][1][0]-halfPipeLength_D, peak_D[2][1][0]+halfPipeLength_D, peak_D[2][1][1], peak_D[1][1][1], data_long[-1]-data_long[0], plot = 1, clim = data_clim, lim_start = lim_start_D)
    ROI(peak_D[1][2][0]-halfPipeLength_D, peak_D[1][2][0]+halfPipeLength_D, peak_D[0][2][1], peak_D[1][2][1], data_long[-1]-data_long[0], plot = 1, clim = data_clim, lim_start = lim_start_D)
    ROI(peak_D[0][2][0]-halfPipeLength_D, peak_D[0][2][0]+halfPipeLength_D, peak_D[0][2][1], peak_D[1][1][1], data_long[-1]-data_long[0], plot = 1, clim = data_clim, lim_start = lim_start_D)
    
# Save the linescan data
np.save(save_address + 'linescan_region1_short.npy', line_D_1_short)
np.save(save_address + 'linescan_region2_short.npy', line_D_2_short)
np.save(save_address + 'linescan_region1_long.npy', line_D_1_long)
np.save(save_address + 'linescan_region2_long.npy', line_D_2_long)

# Make simulation linescan the same length as the data linescan
x = np.linspace(0, np.shape(line_D_1_short)[1]-1, np.shape(line_T1_1)[0])
spline_T1_1 = interp1d(x, line_T1_1, kind = 'cubic')
spline_T2_1 = interp1d(x, line_T2_1, kind = 'cubic') 
spline_L_1 = interp1d(x, line_L_1, kind = 'cubic')
x = np.linspace(0, np.shape(line_D_2_short)[1]-1, np.shape(line_T1_2)[0])
spline_T1_2 = interp1d(x, line_T1_2, kind = 'cubic') 
spline_T2_2 = interp1d(x, line_T1_2, kind = 'cubic') 
spline_L_2 = interp1d(x, line_T1_2, kind = 'cubic') 

if plot_short_region1 == 1:
    # Plot short scan and fits for region 1
    fitParamA_region1_short = []
    fitParamB_region1_short = []
    #fit_lower = 13
    fit_lower = 18
    fit_higher = 192
    for i in range(np.shape(line_D_1_short)[0]):
        plt.figure()
        ax1 = plt.subplot(211)
        plt.plot(line_D_1_short[i], label = 'Data')
        fit_variables = theory_T_fit(range(fit_lower, fit_higher+1), line_D_1_short[i][fit_lower:fit_higher+1])
        
        fitParamA_region1_short.append(fit_variables[0])
        fitParamB_region1_short.append(fit_variables[1])
        fit = theory_T(range(np.shape(line_D_1_short)[1]), *fit_variables)
        plt.plot(fit, linestyle = '--', label = 'Fit (A = ' + '{:.2e}'.format(fit_variables[0]) + '  B = ' + str(round(fit_variables[1], 2)) + ')')
        '''
        fitParamA_region1_short.append(fit_variables[0])
        fit = theory_T(range(np.shape(line_D_1_short)[1]), *fit_variables)
        plt.plot(fit, linestyle = '--', label = 'Fit (A = ' + '{:.2e}'.format(fit_variables[0]))   
        '''
        plt.title(str(int(round(delay_short[i+1]))) + 'PS')
        plt.legend(loc = 'upper center')
        plt.axvline(fit_lower, linestyle = ':', c = 'k')
        plt.axvline(fit_higher, linestyle = ':', c = 'k')
        plt.ylim(-5, 600)
        plt.subplot(212, sharex=ax1)
        plt.title('Simulation-Data')
        plt.plot(range(fit_lower, fit_higher+1),fit[fit_lower:fit_higher+1]-line_D_1_short[i][fit_lower:fit_higher+1])
        plt.axvline(fit_lower, linestyle = ':', c = 'k')
        plt.axvline(fit_higher, linestyle = ':', c = 'k')
        plt.xlim(0, 210)
        plt.ylim(-45,120)
    plt.show()
    np.save(save_address + 'fitParamA_region1_short.npy', fitParamA_region1_short)
    np.save(save_address + 'fitParamB_region1_short.npy', fitParamB_region1_short)

if plot_long_region1 == 1:
    # Plot long scan and fits for region 1
    fitParamA_region1_long = []
    fitParamB_region1_long = []
    fit_lower = 15
    fit_higher = 195
    for i in range(np.shape(line_D_1_long)[0]):
        plt.figure()
        ax1 = plt.subplot(211)
        plt.plot(line_D_1_long[i], label = 'Data')
        fit_variables = theory_T_fit(range(fit_lower, fit_higher+1), line_D_1_long[i][fit_lower:fit_higher+1])
        
        fitParamA_region1_long.append(fit_variables[0])
        fitParamB_region1_long.append(fit_variables[1])
        fit = theory_T(range(np.shape(line_D_1_long)[1]), *fit_variables)
        plt.plot(fit, linestyle = '--', label = 'Fit (A = ' + '{:.2e}'.format(fit_variables[0]) + '  B = ' + str(round(fit_variables[1], 2)) + ')')
        '''
        fitParamA_region1_long.append(fit_variables[0])
        fit = theory_T(range(np.shape(line_D_1_long)[1]), *fit_variables)
        plt.plot(fit, linestyle = '--', label = 'Fit (A = ' + '{:.2e}'.format(fit_variables[0]))        
        '''
        plt.title(str(int(round(delay_long[i+1]))) + 'PS')
        plt.legend(loc = 'upper center')
        plt.axvline(fit_lower, linestyle = ':', c = 'k')
        plt.axvline(fit_higher, linestyle = ':', c = 'k')
        plt.ylim(-5, 600)
        plt.subplot(212, sharex=ax1)
        plt.title('Simulation-Data')
        plt.plot(range(fit_lower, fit_higher+1),fit[fit_lower:fit_higher+1]-line_D_1_long[i][fit_lower:fit_higher+1])
        plt.axvline(fit_lower, linestyle = ':', c = 'k')
        plt.axvline(fit_higher, linestyle = ':', c = 'k')
        plt.xlim(0, 210)
        plt.ylim(-45,120)
    plt.show()
    np.save(save_address + 'fitParamA_region1_long.npy', fitParamA_region1_long)
    np.save(save_address + 'fitParamB_region1_long.npy', fitParamB_region1_long)

if plot_short_region2 == 1:
    # Plot short scan and fits for region 2
    fitParamA_region2_short = []
    fitParamB_region2_short = []
    fit_lower = 18
    fit_higher = 95
    for i in range(np.shape(line_D_2_short)[0]):
        plt.figure()
        ax1 = plt.subplot(211)
        plt.plot(line_D_2_short[i], label = 'Data')
        fit_variables = theory_L_fit(range(fit_lower, fit_higher+1), line_D_2_short[i][fit_lower:fit_higher+1])
        
        fitParamA_region2_short.append(fit_variables[0])
        fitParamB_region2_short.append(fit_variables[1])
        fit = theory_L(range(np.shape(line_D_2_short)[1]), *fit_variables)
        plt.plot(fit, linestyle = '--', label = 'Fit (A = ' + '{:.2e}'.format(fit_variables[0]) + '  B = ' + str(round(fit_variables[1], 2)) + ')')
        '''
        fitParamA_region2_short.append(fit_variables[0])
        fit = theory_L(range(np.shape(line_D_2_short)[1]), *fit_variables)
        plt.plot(fit, linestyle = '--', label = 'Fit (A = ' + '{:.2e}'.format(fit_variables[0]))
        '''
        plt.title(str(int(round(delay_short[i+1]))) + 'PS')
        plt.legend(loc = 'upper center')
        plt.axvline(fit_lower, linestyle = ':', c = 'k')
        plt.axvline(fit_higher, linestyle = ':', c = 'k')
        plt.ylim(-5, 800)
        plt.subplot(212, sharex = ax1)
        plt.title('Simulation-Data')
        plt.plot(range(fit_lower, fit_higher+1),fit[fit_lower:fit_higher+1]-line_D_2_short[i][fit_lower:fit_higher+1])
        plt.axvline(fit_lower, linestyle = ':', c = 'k')
        plt.axvline(fit_higher, linestyle = ':', c = 'k')
        plt.xlim(0, 105)
        plt.ylim(-25,50)

    plt.show()   
    np.save(save_address + 'fitParamA_region2_short.npy', fitParamA_region2_short)
    np.save(save_address + 'fitParamB_region2_short.npy', fitParamB_region2_short)

if plot_long_region2 == 1:
    # Plot short scan and fits for region 2
    fitParamA_region2_long = []
    fitParamB_region2_long = []
    fit_lower = 18
    fit_higher = 95
    for i in range(np.shape(line_D_2_long)[0]):
        plt.figure()
        ax1 = plt.subplot(211)
        plt.plot(line_D_2_long[i], label = 'Data')
        fit_variables = theory_L_fit(range(fit_lower, fit_higher+1), line_D_2_long[i][fit_lower:fit_higher+1])
        
        fitParamA_region2_long.append(fit_variables[0])
        fitParamB_region2_long.append(fit_variables[1])
        fit = theory_L(range(np.shape(line_D_2_long)[1]), *fit_variables)
        plt.plot(fit, linestyle = '--', label = 'Fit (A = ' + '{:.2e}'.format(fit_variables[0]) + '  B = ' + str(round(fit_variables[1], 2)) + ')')
        '''
        fitParamA_region2_long.append(fit_variables[0])
        fit = theory_L(range(np.shape(line_D_2_long)[1]), *fit_variables)
        plt.plot(fit, linestyle = '--', label = 'Fit (A = ' + '{:.2e}'.format(fit_variables[0]))        
        '''
        plt.title(str(int(round(delay_long[i+1]))) + 'PS')
        plt.legend(loc = 'upper center')
        plt.axvline(fit_lower, linestyle = ':', c = 'k')
        plt.axvline(fit_higher, linestyle = ':', c = 'k')
        plt.ylim(-5, 800)
        plt.subplot(212, sharex = ax1)
        plt.title('Simulation-Data')
        plt.plot(range(fit_lower, fit_higher+1),fit[fit_lower:fit_higher+1]-line_D_2_long[i][fit_lower:fit_higher+1])
        plt.axvline(fit_lower, linestyle = ':', c = 'k')
        plt.axvline(fit_higher, linestyle = ':', c = 'k')
        plt.xlim(0, 105)
        plt.ylim(-25,50)
    plt.show()  
    np.save(save_address + 'fitParamA_region2_long.npy', fitParamA_region2_long)
    np.save(save_address + 'fitParamB_region2_long.npy', fitParamB_region2_long)
    
# Combine short and long scans in order of delay time
if plot_fitParams == 1:
    fitParamA_region1 = combine_short_long(fitParamA_region1_short, fitParamA_region1_long)
    fitParamB_region1 = combine_short_long(fitParamB_region1_short, fitParamB_region1_long)
    fitParamA_region2 = combine_short_long(fitParamA_region2_short, fitParamA_region2_long)
    fitParamB_region2 = combine_short_long(fitParamB_region2_short, fitParamB_region2_long)
    time_delay = combine_short_long(delay_short[1:], delay_long[1:])
    # Plot fit parameters as a function of time delay
    fancy_plot(time_delay, fitParamA_region1, fitParamB_region1, 'Region 1')
    fancy_plot(time_delay, fitParamA_region2, fitParamB_region2, 'Region 2')
    plt.figure()
    plt.plot(time_delay, fitParamA_region1, 'bo-', label = 'region 1')
    plt.plot(time_delay, fitParamA_region2, 'ro-', label = 'region 2')
    plt.legend(loc = 'lower right')
    plt.xlabel('Time Delay(ps)')
    plt.ylabel('Fit Parameter A')
    plt.figure()
    plt.plot(time_delay, fitParamB_region1, 'bo-', label = 'region 1')
    plt.plot(time_delay, fitParamB_region2, 'ro-', label = 'region 2')
    plt.legend(loc = 'lower right')
    plt.xlabel('Time Delay(ps)')
    plt.ylabel('Fit Parameter B')
'''
fig, ax1 = plt.subplots()
ax1.plot(time_delay, fitParamA_region1, 'bo-')
ax1.set_xlabel('Time Delay(ps)')
ax1.set_ylabel ('Fit Parameter A', color = 'b')
for i in ax1.get_yticklabels():
    i.set_color('b')
ax2= ax1.twinx()
ax2.plot(time_delay, fitParamB_region1, 'ro-')
ax2.set_ylabel('Fit Parameter B', color = 'r')
for i in ax2.get_yticklabels():
    i.set_color('r')
plt.show
'''    