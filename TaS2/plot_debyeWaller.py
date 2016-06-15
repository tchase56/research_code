import numpy as np
import matplotlib.pyplot as plt 
# User defined functions
import Region_Of_Interest_2
import Fit_Peaks_2

''' Values to change for each run'''
load_address = r'E:\tchase56\TaS2\20160611\TimeScan\scan3\images-ANDOR1\\'
save_address = r'E:\tchase56\TaS2\20160611\TimeScan\scan3\images-ANDOR1\\'
# Delay stage settings
time_zero = 63.2


delayStage_start, delayStage_end, delayStage_step  = [62.90, 63.80, 0.0150] 


# How wide do you want the ROIs for chooseing bragg peaks?
ROI_width = 50      
# Which Debye-Waller peaks would you like Saved/Plotted
debye_flag_2 = 1    
debye_flag_3 = 1
debye_flag_4 = 1   
debye_flag_5 = 0    
debye_flag_6 = 0
# Have the peak intensities been saved? (Have you chosen ROI's Already?)
intensities_saved = 0    
# Y-scale for debye-waller plots
y_min_debye_1 = 0.90
y_min_debye_2 = 0.90
y_min_debye_3 = 0.90
y_min_debye_4 = 0.90
y_min_debye_5 = 0.90
y_min_debye_6 = 0.90
y_max_debye = 1.05
# Plot in units of ps?
ps_flag = 1
# How many of the initial debye-waller points do you want to scale by?
points = 2





# Used to normalize time scans with multiple values below time zero
def Scale(graph, timeStepForNorm):
    norm = np.average(np.array(graph[timeStepForNorm]))
    graph = graph / norm
    return graph

# Paper peak intensity data
paper_address = r"E:\tchase56\Gold_Diffuse\old_data\\"

# Determine delay in ps from delay stage
delay_stage = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
speed_of_light = 299792458 # meters/second
delay = np.round((np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12,2) # Delay in ps
#times_timescan = np.arange(-2.25,17.75,20.0/50.0)

# Load data 
paper_intensity = np.loadtxt(paper_address + 'peak_intensities_timescan')
if intensities_saved == 0:
    pumped = np.load(load_address + 'averaged_aligned.npy')

    # Define peak ROIs
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(pumped[0] ,6, halfLength=ROI_width/2, contrastFactor = 0.25)
    [peak_region_2, background_region_2] = Region_Of_Interest_2.GetRegionOfInterest(pumped[0] ,6, halfLength=ROI_width/2, contrastFactor = 0.25)
    if debye_flag_2 == 1:
        [peak_region_3, background_region_3] = Region_Of_Interest_2.GetRegionOfInterest(pumped[0] ,6, halfLength=ROI_width/2, contrastFactor = 0.25)
        if debye_flag_3 == 1:
            [peak_region_4, background_region_4] = Region_Of_Interest_2.GetRegionOfInterest(pumped[0] ,12, halfLength=ROI_width/2, contrastFactor = 0.25)
            if debye_flag_4 == 1:
                [peak_region_5, background_region_5] = Region_Of_Interest_2.GetRegionOfInterest(pumped[0] ,6, halfLength=ROI_width/2, contrastFactor = 0.25)
                if debye_flag_5 == 1:
                    [peak_region_6, background_region_6] = Region_Of_Interest_2.GetRegionOfInterest(pumped[0] ,4, halfLength=ROI_width/2, contrastFactor = 0.25)
                    if debye_flag_6 == 1:
                        [peak_region_7, background_region_7] = Region_Of_Interest_2.GetRegionOfInterest(pumped[0] ,8, halfLength=ROI_width/2, contrastFactor = 0.25)
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

    # Save peak intensity 
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

# Normalize peaks by innermost peak and scale them to 1
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
# Normalize peaks by innermost peak and scale them to 1
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

old_debye_1 = paper_intensity[6,:]/paper_intensity[5]
old_debye_1 = Scale(old_debye_1, [0,1,2,3,4])
if debye_flag_2 == 1:
    old_debye_2 = paper_intensity[7,:]/paper_intensity[5]
    old_debye_2 = Scale(old_debye_2, [0,1,2,3,4])
    if debye_flag_3 == 1:
        old_debye_3 = paper_intensity[8,:]/paper_intensity[5]
        old_debye_3 = Scale(old_debye_3, [0,1,2,3,4])

# Plot debye-waller
if ps_flag == 1:
    plt.figure()
    plt.plot(delay, new_debye_gaussian_1, 'ro-', label = "New Data")
    plt.title('[220]/[200]')
    plt.xlabel('Time Delay (ps)')
    plt.legend()
    plt.ylim(y_min_debye_1,y_max_debye)
    plt.show()

    if debye_flag_2 == 1:
        plt.figure()
        plt.plot(delay, new_debye_gaussian_2, 'ro-', label = "New Data")
        plt.title('[400]/[200]')
        plt.xlabel('Time Delay (ps)')
        plt.legend()
        plt.ylim(y_min_debye_2,y_max_debye)
        plt.show()

        if debye_flag_3 == 1:
            plt.figure()
            plt.plot(delay, new_debye_gaussian_3, 'ro-', label = "New Data")
            plt.title('[420]/[200]')
            plt.xlabel('Time Delay (ps)')
            plt.legend()
            plt.ylim(y_min_debye_3,y_max_debye)
            plt.show()
        
            if debye_flag_4 == 1:
                plt.figure()
                plt.plot(delay, new_debye_gaussian_4, 'ro-', label = "New Data")
                plt.title('[440]/[200]')
                plt.xlabel('Time Delay (ps)')
                plt.legend()
                plt.ylim(y_min_debye_4,y_max_debye)
                plt.show()
        
                if debye_flag_5 == 1:
                    plt.figure()
                    plt.plot(delay, new_debye_gaussian_5, 'ro-', label = "New Data")
                    plt.title('[600]/[200]')
                    plt.xlabel('Time Delay (ps)')
                    plt.legend()
                    plt.ylim(y_min_debye_5,y_max_debye)
                    plt.show()
                
                    if debye_flag_6 == 1:
                        plt.figure()
                        plt.plot(delay, new_debye_gaussian_6, 'ro-', label = "New Data")
                        plt.title('[620]/[200]')
                        plt.xlabel('Time Delay (ps)')
                        plt.legend()
                        plt.ylim(y_min_debye_5,y_max_debye)
                        plt.show()
else:
    plt.figure()
    plt.plot(delay_stage, new_debye_gaussian_1, 'ro-', label = "New Data")
    plt.title('[220/200]')
    plt.xlabel('Delay Stage')
    plt.legend()
#    plt.ylim(y_min_debye_1,y_max_debye)
    plt.show()

    if debye_flag_2 == 1:
        plt.figure()
        plt.plot(delay_stage, new_debye_gaussian_2, 'ro-', label = "New Data")
        plt.title('[400]/[200]')
        plt.xlabel('Delay Stage')
        plt.legend()
#        plt.ylim(y_min_debye_2,y_max_debye)
        plt.show()

        if debye_flag_3 == 1:
            plt.figure()
            plt.plot(delay_stage, new_debye_gaussian_3, 'ro-', label = "New Data")
            plt.title('[420]/[200]')
            plt.xlabel('Delay Stage')
            plt.legend()
#            plt.ylim(y_min_debye_3,y_max_debye)
            plt.show()
        
            if debye_flag_4 == 1:
                plt.figure()
                plt.plot(delay_stage, new_debye_gaussian_4, 'ro-', label = "New Data")
                plt.title('[440]/[200]')
                plt.xlabel('Delay Stage')
                plt.legend()
#                plt.ylim(y_min_debye_4,y_max_debye)
                plt.show()
        
                if debye_flag_5 == 1:
                    plt.figure()
                    plt.plot(delay_stage, new_debye_gaussian_5, 'ro-', label = "New Data")
                    plt.title('[600]/[200]')
                    plt.xlabel('Delay Stage')
                    plt.legend()
#                    plt.ylim(y_min_debye_5,y_max_debye)
                    plt.show()
                
                    if debye_flag_6 == 1:
                        plt.figure()
                        plt.plot(delay_stage, new_debye_gaussian_6, 'ro-', label = "New Data")
                        plt.title('[620]/[200]')
                        plt.xlabel('Delay Stage')
                        plt.legend()
#                        plt.ylim(y_min_debye_5,y_max_debye)
                        plt.show()   