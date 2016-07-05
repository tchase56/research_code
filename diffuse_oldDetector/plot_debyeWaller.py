# OLD DETECTOR
'''
Tyler Chase
06/27/2016
This module reads in the output files from align_average.py or align_different_scans.py
and analyzes the change in Bragg peak intensity as a function of time delay

Instructions
0. Make sure you have already run align_average.py or align_different_scans.py
   (set different_scans flag accordingly)
1. Change addresses
2. Change delay stage values and make sure time zero is correct
3. Flag which Bragg peaks you would like to plot (normalized by innermost)
4. Ones you have chosen peaks once intensities_saved flag can be set to 1 to 
   avoid having to choose all of the peaks again
5. If you want in x axis picoseconds set ps_flag = 1 otherwise units of delay stage
6. y_min and y_max values set the y-scale of the plots for debye-waller
7. The drop is normalized so you see percent change
   points is the number of time delays you are taking average of and normalizing
   the average to 1
8. Run the code
9. Click on the innermost bragg peaks
10. Click next to the innermost bragg peaks (background)
11. Iteritively repeat steps 9 and 10 (for next innermost bragg peaks) 
until you are no longer prompted to choose more peaks
'''

''' Values to change for each run'''
load_address = r'E:\Klaus\20160701\Ni_scan12\images-ANDOR1\\'
save_address = r'E:\Klaus\20160701\Ni_scan12\images-ANDOR1\\'
# Delay stage settings
time_zero = 64.47
#delayStage_start, delayStage_end, delayStage_step  = [64.15, 66.15, 0.02]
#delayStage_start, delayStage_end, delayStage_step  = [64.3, 65.5, 0.02]  
#delayStage_start, delayStage_end, delayStage_step  = [62.0, 68.0, 0.1]
delayStage_start, delayStage_end, delayStage_step  = [64.4, 65.2, 0.02]

# How wide do you want the ROIs for chooseing bragg peaks?
ROI_width = 64      
# Which Debye-Waller peaks would you like Saved/Plotted
debye_flag_2 = 0    
debye_flag_3 = 0
debye_flag_4 = 0   
debye_flag_5 = 0    
debye_flag_6 = 0
# Have the peak intensities been saved? (Have you chosen ROI's Already?)
intensities_saved = 0
# Y-scale for debye-waller plots
y_min_debye_1 = 0.925
y_min_debye_2 = 0.95
y_min_debye_3 = 0.92
y_min_debye_4 = 0.90
y_min_debye_5 = 0.90
y_min_debye_6 = 0.90
y_max_debye = 1.015
# Plot in units of ps?
ps_flag = 0
# How many of the initial debye-waller points do you want to scale by?
points = 4
# Is this the output of align_average or the output of align_different_scans
different_scans = 0










import numpy as np
import matplotlib.pyplot as plt 
# User defined functions
import Region_Of_Interest_2
import Fit_Peaks_2

# Used to normalize time scans with multiple values below time zero
def Scale(graph, timeStepForNorm):
    norm = np.average(np.array(graph[timeStepForNorm]))
    graph = graph / norm
    return graph

# Determine delay in ps from delay stage
delay_stage = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
speed_of_light = 299792458 # meters/second
delay = np.round((np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12,2) # Delay in ps
#times_timescan = np.arange(-2.25,17.75,20.0/50.0)

# Load data 
if intensities_saved == 0:
    if different_scans == 0:
        pumped = np.load(load_address + 'averaged.npy')
    else:
        pumped = np.load(load_address + 'averaged_runs.npy')
    pumped_ave = np.average(pumped, 0)
    # Define peak ROIs
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25)
    [peak_region_2, background_region_2] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25)
    if debye_flag_2 == 1:
        [peak_region_3, background_region_3] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25)
        if debye_flag_3 == 1:
            [peak_region_4, background_region_4] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,8, halfLength=ROI_width/2, contrastFactor = 0.25)
            if debye_flag_4 == 1:
                [peak_region_5, background_region_5] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25)
                if debye_flag_5 == 1:
                    [peak_region_6, background_region_6] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25)
                    if debye_flag_6 == 1:
                        [peak_region_7, background_region_7] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,8, halfLength=ROI_width/2, contrastFactor = 0.25)
    # Number of rows to save in intensity file
    length = 2
    if debye_flag_2 == 1:
        length = 3
        if debye_flag_3 == 1:
            length = 4
            if debye_flag_4 == 1:
                length = 5
                if debye_flag_5 == 1:
                    length = 6
                    if debye_flag_6 == 1:
                        length = 7

    # Save peak intensity (both intensity of ROI and Gaussian fit amplitude)
    peak_intensity = np.zeros((length,np.shape(pumped)[0]))
    peak_intensity_gaussian = np.zeros((length, np.shape(pumped)[0]))
    for i in range(0,np.shape(pumped)[0]):    
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region, background_region, halfLength=ROI_width/2)
        peak_intensity[0,i] = np.average(peak_intensity_pump) - np.average(background_intensity)
        peak_intensity_gaussian[0,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_2, background_region_2, halfLength=ROI_width/2)
        peak_intensity[1,i] = np.average(peak_intensity_pump) - np.average(background_intensity)   
        peak_intensity_gaussian[1,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
        if debye_flag_2 == 1:
            [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_3, background_region_3, halfLength=ROI_width/2)
            peak_intensity[2,i] = np.average(peak_intensity_pump) - np.average(background_intensity)
            peak_intensity_gaussian[2,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
            if debye_flag_3 == 1:
                [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_4, background_region_4, halfLength=ROI_width/2)
                peak_intensity[3,i] = np.average(peak_intensity_pump) - np.average(background_intensity)  
                peak_intensity_gaussian[3,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
                if debye_flag_4 == 1:
                    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_5, background_region_5, halfLength=ROI_width/2)
                    peak_intensity[4,i] = np.average(peak_intensity_pump) - np.average(background_intensity) 
                    peak_intensity_gaussian[4,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
                    if debye_flag_5 == 1:
                        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_6, background_region_6, halfLength=ROI_width/2)
                        peak_intensity[5,i] = np.average(peak_intensity_pump) - np.average(background_intensity) 
                        peak_intensity_gaussian[5,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
                        if debye_flag_6 == 1:
                            [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_7, background_region_7, halfLength=ROI_width/2)
                            peak_intensity[6,i] = np.average(peak_intensity_pump) - np.average(background_intensity) 
                            peak_intensity_gaussian[0,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
    np.save(save_address + "peak_intensity", peak_intensity)
    np.save(save_address + "peak_intensity_gaussian", peak_intensity_gaussian)
else:
    peak_intensity = np.load(save_address + 'peak_intensity.npy')
    peak_intensity_gaussian = np.load(save_address + 'peak_intensity_gaussian.npy')

# Normalize peaks by innermost peak and scale them to 1 (for gaussian intensity)
new_debye_gaussian_1 = peak_intensity_gaussian[1,:]/peak_intensity_gaussian[0,:]
new_debye_gaussian_1 = Scale(new_debye_gaussian_1, list(range(points)))
if debye_flag_2 == 1:
    new_debye_gaussian_2 = peak_intensity_gaussian[2,:]/peak_intensity_gaussian[0,:]
    new_debye_gaussian_2 = Scale(new_debye_gaussian_2, list(range(points)))
    if debye_flag_3 == 1:
        new_debye_gaussian_3 = peak_intensity_gaussian[3,:]/peak_intensity_gaussian[0,:]
        new_debye_gaussian_3 = Scale(new_debye_gaussian_3, list(range(points)))
        if debye_flag_4 == 1:
            new_debye_gaussian_4 = peak_intensity_gaussian[4,:]/peak_intensity_gaussian[0,:]
            new_debye_gaussian_4 = Scale(new_debye_gaussian_4, list(range(points)))
            if debye_flag_5 == 1:
                new_debye_gaussian_5 = peak_intensity_gaussian[5,:]/peak_intensity_gaussian[0,:]
                new_debye_gaussian_5 = Scale(new_debye_gaussian_5, list(range(points)))
                if debye_flag_6 == 1:
                    new_debye_gaussian_6 = peak_intensity_gaussian[6,:]/peak_intensity_gaussian[0,:]
                    new_debye_gaussian_6 = Scale(new_debye_gaussian_6, list(range(points)))
                    
'''                   
# Normalize peaks by innermost peak and scale them to 1 (for sum of ROI intensity)
new_debye_gaussian_1 = peak_intensity[1,:]/peak_intensity[0,:]
new_debye_gaussian_1 = Scale(new_debye_gaussian_1, list(range(points)))
if debye_flag_2 == 1:
    new_debye_gaussian_2 = peak_intensity[2,:]/peak_intensity[0,:]
    new_debye_gaussian_2 = Scale(new_debye_gaussian_2, list(range(points)))
    if debye_flag_3 == 1:
        new_debye_gaussian_3 = peak_intensity[3,:]/peak_intensity[0,:]
        new_debye_gaussian_3 = Scale(new_debye_gaussian_3, list(range(points)))
        if debye_flag_4 == 1:
            new_debye_gaussian_4 = peak_intensity[4,:]/peak_intensity[0,:]
            new_debye_gaussian_4 = Scale(new_debye_gaussian_4, list(range(points)))
            if debye_flag_5 == 1:
                new_debye_gaussian_5 = peak_intensity[5,:]/peak_intensity[0,:]
                new_debye_gaussian_5 = Scale(new_debye_gaussian_5, list(range(points)))
                if debye_flag_6 == 1:
                    new_debye_gaussian_6 = peak_intensity[6,:]/peak_intensity[0,:]
                    new_debye_gaussian_6 = Scale(new_debye_gaussian_6, list(range(points)))
'''

# Plot debye-waller (in units of picoseconds)
if ps_flag == 1:
    plt.figure()
    plt.plot(delay, new_debye_gaussian_1, 'ro-', label = "New Data")
    plt.title('SecondOrder/FirstOrder')
    plt.xlabel('Time Delay (ps)')
    plt.legend()
    plt.ylim(y_min_debye_1,y_max_debye)
    plt.show()

    if debye_flag_2 == 1:
        plt.figure()
        plt.plot(delay, new_debye_gaussian_2, 'ro-', label = "New Data")
        plt.title('ThirdOrder/FirstOrder')
        plt.xlabel('Time Delay (ps)')
        plt.legend()
        plt.ylim(y_min_debye_2,y_max_debye)
        plt.show()

        if debye_flag_3 == 1:
            plt.figure()
            plt.plot(delay, new_debye_gaussian_3, 'ro-', label = "New Data")
            plt.title('FourthOrder/FirstOrder')
            plt.xlabel('Time Delay (ps)')
            plt.legend()
            plt.ylim(y_min_debye_3,y_max_debye)
            plt.show()
        
            if debye_flag_4 == 1:
                plt.figure()
                plt.plot(delay, new_debye_gaussian_4, 'ro-', label = "New Data")
                plt.title('FifthOrder/FirstOrder')
                plt.xlabel('Time Delay (ps)')
                plt.legend()
                plt.ylim(y_min_debye_4,y_max_debye)
                plt.show()
        
                if debye_flag_5 == 1:
                    plt.figure()
                    plt.plot(delay, new_debye_gaussian_5, 'ro-', label = "New Data")
                    plt.title('sat/FirstOrder')
                    plt.xlabel('Time Delay (ps)')
                    plt.legend()
                    plt.ylim(y_min_debye_5,y_max_debye)
                    plt.show()
                
                    if debye_flag_6 == 1:
                        plt.figure()
                        plt.plot(delay, new_debye_gaussian_6, 'ro-', label = "New Data")
                        plt.title('SeventhOrder/FirstOrder')
                        plt.xlabel('Time Delay (ps)')
                        plt.legend()
                        plt.ylim(y_min_debye_5,y_max_debye)
                        plt.show()
# Plot debye-waller (in units of delay stage distance)
else:
    plt.figure()
    plt.plot(delay_stage, new_debye_gaussian_1, 'ro-', label = "New Data")
    plt.title('SecondOrder/FirstOrder')
    plt.xlabel('Delay Stage')
    plt.legend()
    plt.ylim(y_min_debye_1,y_max_debye)
    plt.show()

    if debye_flag_2 == 1:
        plt.figure()
        plt.plot(delay_stage, new_debye_gaussian_2, 'ro-', label = "New Data")
        plt.title('ThirdOrder/FirstOrder')
        plt.xlabel('Delay Stage')
        plt.legend()
        plt.ylim(y_min_debye_2,y_max_debye)
        plt.show()

        if debye_flag_3 == 1:
            plt.figure()
            plt.plot(delay_stage, new_debye_gaussian_3, 'ro-', label = "New Data")
            plt.title('FourthOrder/FirstOrder')
            plt.xlabel('Delay Stage')
            plt.legend()
            plt.ylim(y_min_debye_3,y_max_debye)
            plt.show()
        
            if debye_flag_4 == 1:
                plt.figure()
                plt.plot(delay_stage, new_debye_gaussian_4, 'ro-', label = "New Data")
                plt.title('FifthOrder/FirstOrder')
                plt.xlabel('Delay Stage')
                plt.legend()
                plt.ylim(y_min_debye_4,y_max_debye)
                plt.show()
        
                if debye_flag_5 == 1:
                    plt.figure()
                    plt.plot(delay_stage, new_debye_gaussian_5, 'ro-', label = "New Data")
                    plt.title('sat/FirstOrder')
                    plt.xlabel('Delay Stage')
                    plt.legend()
                    plt.ylim(y_min_debye_5,y_max_debye)
                    plt.show()
                
                    if debye_flag_6 == 1:
                        plt.figure()
                        plt.plot(delay_stage, new_debye_gaussian_6, 'ro-', label = "New Data")
                        plt.title('SeventhOrder/FirstOrder')
                        plt.xlabel('Delay Stage')
                        plt.legend()
                        plt.ylim(y_min_debye_5,y_max_debye)
                        plt.show()   