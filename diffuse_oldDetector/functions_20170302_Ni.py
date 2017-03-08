'''
Tyler Chase
09/07/2016
- functions for symmetrizing diffraction pattern
- function for masking out edges and center of phosphor screen
- function for forming linescan and plotting ROI
- function for converting from delay stage to picoseconds
- function for expanding using the 2-D interpolation
- function for calculating the peak separation for a given order
- function for loading and aligning timescan
'''





# Diffuse Tools
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
import matplotlib.patches as patches
from scipy import interpolate
import glob
import scipy.optimize as opt
import copy
# Functions for aligning using cross correlaton
from skimage.feature import register_translation
from scipy.ndimage.interpolation import shift
# User defined function
import Region_Of_Interest_2
import Fit_Peaks_2


# Function for symmetrizing an image with six fold rotation symmetry
def symmetrize(image, fold):
    rotated = []
    for i in range(fold):
        rotated.append(image)
        image = rotate(image, 360/fold, reshape = False)
    symmetrized = np.average(rotated, 0)
    return(symmetrized)
    
# Old method of symmetrization with flipping
def symmetrize_flip(image):
    temp = np.flipud(image)
    temp_2 = np.average([image, temp],0)
    temp_3 = np.fliplr(temp_2)
    symmetrized = np.average([temp_3,temp_2],0)
    return symmetrized
    
# Create mask and mask image
def make_mask(image, outside, inside):
    mask = np.zeros(np.shape(image))
    mask[outside:inside,outside:-1-outside] = 1
    mask[-1-inside:-1-outside,outside:-1-outside] = 1
    mask[outside:-1-outside, outside:inside] = 1
    mask[outside:-1-outside, -1-inside:-1-outside] = 1
    masked = mask
    return masked
    
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
'''   
def ps_convert(start, stop, step, time_zero):
    # Delay Stage Settings
    delay_stage = np.arange(start, stop + (step)/2, step)
    speed_of_light = 299792458 # meters/second
    delay = (np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12 # Delay in ps
    return(delay)
'''
def ps_convert(delay_stage, time_zero):
    # Delay Stage Settings
    speed_of_light = 299792458 # meters/second
    delay = (np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12 # Delay in ps
    return(delay)
    
def delay_convert(ps, time_zero):
    # Delay Stage Settings
    speed_of_light = 299792458 # meters/second
    ps = np.float32(ps)
    delay = ps * speed_of_light/(2*10**9)+time_zero
    return(delay)

    
'''expand() takes in an image and a factor, then returns an image expanded by 
the input factor (this returned image has the same pixel count as the input image)
'''
def expand(image, factor, kind = 'cubic'):
    # Form 2-D interpolation of the input image
    x_len = np.shape(image)[0]
    y_len = np.shape(image)[1]    
    x = np.arange(x_len)
    y = np.arange(y_len)
    f = interpolate.interp2d(x,y,image, kind = kind, fill_value = 0)    
    # Use the interpolation to set the appropriate pixel count for the expanded image
    x_start = (x_len-x_len/factor)/2
    x_stop = x_len-((x_len-x_len/factor)/2)-1
    y_start = (y_len-y_len/factor)/2
    y_stop = y_len-((y_len-y_len/factor)/2)-1
    x_2 = np.linspace(x_start, x_stop, x_len)
    y_2 = np.linspace(y_start, y_stop, y_len)
    z_2 = f(x_2, y_2)
    return(z_2)
    
# Calculate average peak expansion of bragg peaks given x and y positions 
def peak_expansion_f(array_x, array_y):
    length = len(array_x)
    distance = []
    for i in range(length//2):
        distance.append(np.sqrt((array_x[i]-array_x[i+length//2])**2 + (array_y[i]-array_y[i+length//2])**2))
    distance_ave = np.average(distance,0)
    return(distance_ave)

'''    
def load_timescan(load_address, delayStage_start, delayStage_end, delayStage_step, align = 0, center = 0, 
                  ROI_width = 64, pixel_accuracy = 100, nbpeaks = 4, peak_region = 0, background_region = 0):
    # Make vector of delay stage values
    delay = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
    
    if align == 0:
        # Average unaligned scans and save as numpy arrays
        scan_times = []
        for i in delay:
            scan = []
            # Load sum of images
            for j in glob.glob(load_address + "*%08.4f" % i + '*.tif'):
                scan.append(np.array(plt.imread(j)))
            print('delay stage ' + str(i))
            scan_times.append(np.average(np.array(scan),0))     # creating list of unaligned averaged images        
        return scan_times

    else:
        # Choose ROIs for nbpeaks lowest order peaks for alignment
        temp = plt.imread(glob.glob(load_address + '*.tif')[0])
        if peak_region == 0 and background_region == 0:
            [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(temp , nbpeaks, 
                                            halfLength=ROI_width/2, contrastFactor = 0.25, message = "choose innermost peaks for alignment")        
        
        # Average unaligned scans and save as numpy arrays
        scan_times = []
        scan_times_aligned = []
        for i in delay:
            scan = []
            # Load sum of images
            for j in glob.glob(load_address + "*%08.4f" % i + '*.tif'):
                scan.append(np.array(plt.imread(j)))
            print('delay stage ' + str(i))
            scan_times.append(np.average(np.array(scan),0))     # creating list of unaligned averaged images         
        
            # Align Images and print percentage complete for each delay stage value
            percent_1 = 0  
            scan_2 = []
            scan_2.append(scan[0])
            for k in range(1,len(scan)):
                offset_list = [register_translation(scan[0][int(peak_region[l][0]):int(peak_region[l][1]), int(peak_region[l][2]): int(peak_region[l][3])], 
                                                    scan[k][int(peak_region[l][0]):int(peak_region[l][1]), int(peak_region[l][2]): int(peak_region[l][3])], 
                                                    pixel_accuracy)[0] for l in range(0,nbpeaks)]
                offset = np.average(offset_list,0)
                scan_2.append(shift(scan[k], offset))
                percent_2 = round(k/len(scan)*100)             
                if percent_2 > percent_1:        
                    print(str(round(k/len(scan)*100)) + '%')
                    percent_1 = percent_2
            scan_times_aligned.append(np.average(np.array(scan_2),0))     # creating list of aligned averaged images

        # If there is more than one time delay align the time delays
        if len(scan_times_aligned) > 1:
            scan_times_aligned_2 = []
            scan_times_aligned_2.append(scan_times_aligned[0])
            for l in range(1,len(scan_times_aligned)):
                offset_list = [register_translation(scan_times_aligned[0][int(peak_region[m][0]):int(peak_region[m][1]), int(peak_region[m][2]): int(peak_region[m][3])], 
                                                scan_times_aligned[l][int(peak_region[m][0]):int(peak_region[m][1]), int(peak_region[m][2]): int(peak_region[m][3])], 
                                                pixel_accuracy)[0] for m in range(0,nbpeaks)]
                offset = np.average(offset_list,0)
                scan_times_aligned_2.append(shift(scan_times_aligned[l], offset))
            scan_times_aligned = scan_times_aligned_2
        
        if center == 1:
            image_center = np.shape(scan_times_aligned[0])[0]/2 + 0.5 -1
            diffraction_shift = []
            diffraction_center = []
            for i in range(0,np.shape(scan_times_aligned)[0]):
                [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(scan_times_aligned[i], peak_region, background_region, halfLength=ROI_width/2)
                y_diffraction = np.mean([peak_fit_pump[l][0][1] for l in range(nbpeaks)])
                x_diffraction = np.mean([peak_fit_pump[l][1][1] for l in range(nbpeaks)])
                diffraction_shift.append([(image_center-x_diffraction),(image_center-y_diffraction)])
                diffraction_center.append([x_diffraction,y_diffraction])
            runs_centered = []    
            for i in range(0,np.shape(scan_times_aligned)[0]):
                runs_centered.append(shift(scan_times_aligned[i],diffraction_shift[i])) 
            scan_times_aligned = runs_centered

        return scan_times_aligned, peak_region, background_region
'''

def load_timescan(load_address, delay, align = 0, center = 0, 
                  ROI_width = 64, pixel_accuracy = 100, nbpeaks = 4, peak_region = 0, background_region = 0):
    
    # Make vector of delay stage values
    if align == 0:
        # Average unaligned scans and save as numpy arrays
        scan_times = []
        for i in delay:
            scan = []
            # Load sum of images
            for j in glob.glob(load_address + "*%08.4f" % i + '*.tif'):
                scan.append(np.array(plt.imread(j)))
            print('delay stage ' + str(i))
            scan_times.append(np.average(np.array(scan),0))     # creating list of unaligned averaged images        
        return scan_times

    else:
        # Choose ROIs for nbpeaks lowest order peaks for alignment
        temp = plt.imread(glob.glob(load_address + '*.tif')[0])
        if peak_region == 0 and background_region == 0:
            [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(temp , nbpeaks, 
                                            halfLength=ROI_width/2, contrastFactor = 0.25, message = "choose innermost peaks for alignment")        
        
        # Average unaligned scans and save as numpy arrays
        scan_times = []
        scan_times_aligned = []
        for i in delay:
            scan = []
            # Load sum of images
            for j in glob.glob(load_address + "*%08.4f" % i + '*.tif'):
                scan.append(np.array(plt.imread(j)))
            print('delay stage ' + str(i))
            scan_times.append(np.average(np.array(scan),0))     # creating list of unaligned averaged images         
        
            # Align Images and print percentage complete for each delay stage value
            percent_1 = 0  
            scan_2 = []
            scan_2.append(scan[0])
            for k in range(1,len(scan)):
                offset_list = [register_translation(scan[0][int(peak_region[l][0]):int(peak_region[l][1]), int(peak_region[l][2]): int(peak_region[l][3])], 
                                                    scan[k][int(peak_region[l][0]):int(peak_region[l][1]), int(peak_region[l][2]): int(peak_region[l][3])], 
                                                    pixel_accuracy)[0] for l in range(0,nbpeaks)]
                offset = np.average(offset_list,0)
                scan_2.append(shift(scan[k], offset))
                percent_2 = round(k/len(scan)*100)             
                if percent_2 > percent_1:        
                    print(str(round(k/len(scan)*100)) + '%')
                    percent_1 = percent_2
            scan_times_aligned.append(np.average(np.array(scan_2),0))     # creating list of aligned averaged images

        # If there is more than one time delay align the time delays
        if len(scan_times_aligned) > 1:
            scan_times_aligned_2 = []
            scan_times_aligned_2.append(scan_times_aligned[0])
            for l in range(1,len(scan_times_aligned)):
                offset_list = [register_translation(scan_times_aligned[0][int(peak_region[m][0]):int(peak_region[m][1]), int(peak_region[m][2]): int(peak_region[m][3])], 
                                                scan_times_aligned[l][int(peak_region[m][0]):int(peak_region[m][1]), int(peak_region[m][2]): int(peak_region[m][3])], 
                                                pixel_accuracy)[0] for m in range(0,nbpeaks)]
                offset = np.average(offset_list,0)
                scan_times_aligned_2.append(shift(scan_times_aligned[l], offset))
            scan_times_aligned = scan_times_aligned_2
        
        if center == 1:
            image_center = np.shape(scan_times_aligned[0])[0]/2 + 0.5 -1
            diffraction_shift = []
            diffraction_center = []
            for i in range(0,np.shape(scan_times_aligned)[0]):
                [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(scan_times_aligned[i], peak_region, background_region, halfLength=ROI_width/2)
                y_diffraction = np.mean([peak_fit_pump[l][0][1] for l in range(nbpeaks)])
                x_diffraction = np.mean([peak_fit_pump[l][1][1] for l in range(nbpeaks)])
                diffraction_shift.append([(image_center-x_diffraction),(image_center-y_diffraction)])
                diffraction_center.append([x_diffraction,y_diffraction])
            runs_centered = []    
            for i in range(0,np.shape(scan_times_aligned)[0]):
                runs_centered.append(shift(scan_times_aligned[i],diffraction_shift[i])) 
            scan_times_aligned = runs_centered

        return scan_times_aligned, peak_region, background_region
        
def center_run(runs, ROI_width = 64, nbpeaks = 4):
    # Center the image (using the average postion of the nbpeaks innermost peaks)
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs[0] ,nbpeaks, 
        halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Select Innermost Peaks for centering')
    image_center = np.shape(runs[0])[0]/2 + 0.5 -1
    diffraction_shift = []
    diffraction_center = []
    for i in range(0,np.shape(runs)[0]):
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs[i], peak_region, background_region, halfLength=ROI_width/2)
        y_diffraction = np.mean([peak_fit_pump[l][0][1] for l in range(nbpeaks)])
        x_diffraction = np.mean([peak_fit_pump[l][1][1] for l in range(nbpeaks)])
        diffraction_shift.append([(image_center-x_diffraction),(image_center-y_diffraction)])
        diffraction_center.append([x_diffraction,y_diffraction])
    runs_centered = []    
    for i in range(0,np.shape(runs)[0]):
        runs_centered.append(shift(runs[i],diffraction_shift[i]))   
    return runs_centered
        
def rotate_square(run_aligned_centered, ROI_width = 64):
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(run_aligned_centered[0] ,2, 
        halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Select 2 peaks that you would like to lie flat in rotated image')
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(run_aligned_centered[0], peak_region, background_region, halfLength=ROI_width/2)
    deltaY = peak_fit_pump[1][1][1] - peak_fit_pump[0][1][1]
    deltaX = peak_fit_pump[1][0][1] - peak_fit_pump[0][0][1]
    angle = (360/(2*np.pi))*np.arctan(deltaY/deltaX)
    rotated = []
    for i in range(0,np.shape(run_aligned_centered)[0]):
        rotated.append(rotate(run_aligned_centered[i], angle, reshape=False))
    return rotated
    
# Mask for the detector background (allowing only detector background through)   
def radial_mask(image, radius_outer, ave_x = 511.5, ave_y = 511.5):
    mask = np.ones(np.shape(image))
    for i in range(0, np.shape(image)[0]):
        for j in range(0, np.shape(image)[1]):
            if np.sqrt((i-ave_y)**2 + (j-ave_x)**2)<radius_outer:
                mask[i,j] = 0
    return mask
    
# Load in the runs in chronological order
#  (-6, -2, 2, 6, 10, 14, 18) + (-25, -5, 20, 45, 70, 95) -> (-25, -6, -5, -2, 2, 6, 10, 14, 18, 20, 45, 70, 95)
# or 
#  (-2, 2, 6, 10, 14, 18) + (-5, 20, 45, 70, 95) -> (-5, -2, 2, 6, 10, 14, 18, 20, 45, 70, 95)
def combine_short_long(short, long):
    if len(short) == 7 and len(long) == 6:
        combined = []
        combined.append(long[0])
        combined.append(short[0])
        combined.append(long[1])
        for i in range(1,len(short)):
            combined.append(short[i])
        for i in range(2,len(long)):
            combined.append(long[i])
        return combined
    elif len(short) == 6 and len(long) == 5:
        combined = []
        combined.append(long[0])
        for i in range(0,len(short)):
            combined.append(short[i])
        for i in range(1,len(long)):
            combined.append(long[i])
        return combined
    else:
        raise NameError('inputs are wrong length')
        
        
#  (-25, -6, -5, -2, 2, 6, 10, 14, 18, 20, 45, 70, 95) -> (-6, -2, 2, 6, 10, 14, 18) + (-25, -5, 20, 45, 70, 95) 
# or 
#  (-5, -2, 2, 6, 10, 14, 18, 20, 45, 70, 95) -> (-2, 2, 6, 10, 14, 18) + (-5, 20, 45, 70, 95) 
def separate_short_long(combined):
    if len(combined) == 13:
        short = []
        short.append(combined[1])
        for i in range(3,9):
            short.append(combined[i])
        long = []
        long.append(combined[0])
        long.append(combined[2])
        for i in range(9,13):
            long.append(combined[i])
        return short, long
    if len(combined) == 11:
        short = []
        for i in range(1,7):
            short.append(combined[i])
        long = []
        long.append(combined[0])
        for i in range(7,11):
            long.append(combined[i])
        return short, long
    else:
        raise NameError('input is wrong length')
    
def norm_short_long(combined):
    short, long = separate_short_long(combined)
    intens_factor_s = np.average(short[0:2])
    intens_factor_l = np.average(long[0:2])
    short = np.array(short)/intens_factor_s
    short = np.ndarray.tolist(short)
    long = np.array(long)/intens_factor_l
    long = np.ndarray.tolist(long)
    combined_2 = combine_short_long(short, long)
    return combined_2, intens_factor_s, intens_factor_l
    
def norm_short_long_image(scans, scans_cen_rot):
    inten = screen_intensity(scans_cen_rot)
    inten_norm, factor_s, factor_l = norm_short_long(inten)
    # Normalize images
    scans_s, scans_l = separate_short_long(scans)
    for i in range(len(scans_s)):
        scans_s[i] = scans_s[i]/factor_s
    for i in range(len(scans_l)):
        scans_l[i] = scans_l[i]/factor_l
    scans = combine_short_long(scans_s, scans_l)
    scans_rot_s, scans_rot_l = separate_short_long(scans_cen_rot)
    for i in range(len(scans_rot_s)):
        scans_rot_s[i] = scans_rot_s[i]/factor_s
    for i in range(len(scans_rot_l)):
        scans_rot_l[i] = scans_rot_l[i]/factor_l
    scans_rot = combine_short_long(scans_rot_s, scans_rot_l)
    return scans, scans_rot    

def subtract_detector_background(scans, scans_cen_rot, radius_mask = 505):
    # determine background detector fluctuations per pixel
    mask = radial_mask(scans[0], radius_mask)
    mask_area = np.sum(mask)
    background_intensity = []
    for i in range(len(scans)):
        background_intensity.append(np.sum(scans[i]*mask)/mask_area)
    # Subtract background from images
    scans_2 = []
    scans_cen_rot_2 = []
    for i in range(len(scans_cen_rot)):
        scans_cen_rot_2.append(scans_cen_rot[i] - background_intensity[i])
        scans_2.append(scans[i] - background_intensity[i])
        
    f, (ax1, ax2) = plt.subplots(1,2)
    ax1.imshow(scans[0]*mask)
    ax2.imshow(scans[0]*( -1*(mask-1) ) )
    plt.title('background mask test')
    return scans_2, scans_cen_rot_2, background_intensity
        
def screen_intensity(scans_cen_rot, edge = 133, center = 66):
    # Determine screen intensity fluctuations per pixel
    screen_intensity = []
    mask_2 = make_mask(scans_cen_rot[0], edge, int(1024/2-center))
    mask_2_area = np.sum(mask_2)
    for i in range(len(scans_cen_rot)):
        screen_intensity.append(np.sum(mask_2 * scans_cen_rot[i])/mask_2_area)
    return screen_intensity
    
def norm_screen(runs, edge = 133, center = 66):
    # Normalize Screen Intensity
    screen_inten_norm = []
    norm = []
    screen_inten = screen_intensity(runs, edge = edge, center = center)
    sum_0 = screen_inten[0]
    norm.append(runs[0])
    screen_inten_norm.append(sum_0)
    for i in range(1,np.shape(runs)[0]):
        sum_i = screen_inten[i]
        factor = sum_0/sum_i
        screen_inten_norm.append(sum_i*factor)
        norm.append(runs[i]*factor)
    return screen_inten, screen_inten_norm, norm
    
def expansion_correction(runs_centered, ROI_width = 64, nbpeaks = 4):
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs_centered[0] ,nbpeaks, 
        halfLength=(ROI_width/2), contrastFactor = 0.25, message = 'Select Innermost Peaks for expansion interpolation')
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs_centered[0], peak_region, background_region, halfLength=ROI_width/2)
    peak_distance_ref = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])
    peak_distance = []
    peak_distance.append(peak_distance_ref)
    peak_distance_expanded = []
    peak_distance_expanded.append(peak_distance_ref)
    runs_expanded = []
    runs_expanded.append(runs_centered[0])
    for i in range(1,np.shape(runs_centered)[0]):
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs_centered[i], peak_region, background_region, halfLength=ROI_width/2)
        peak_distance.append(peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1]))
        factor = peak_distance[0]/peak_distance[i]
        expanded = expand(runs_centered[i],factor)
        runs_expanded.append(expanded)
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(expanded, peak_region, background_region, halfLength=ROI_width/2)
        peak_distance_expanded.append(peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1]))
    return peak_distance, peak_distance_expanded, runs_expanded
    
def recyprocal_lattice_points(pattern, ROI_width = 50, halfFlag = 0):
    # Take single diffraction pattern (should be before time zero)
    # Use user defined functions to choose two peaks in gamma, k, x, k, gamma direction
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(pattern ,2, 
        halfLength=(ROI_width/2), contrastFactor = 0.25, 
        message = 'choose two vertical peaks along $\Gamma, K, X, K, \Gamma$ direction')
    [peak_fit, peak_intensity, background_intensity] = Fit_Peaks_2.FitPeaks(pattern,
        peak_region, background_region, halfLength=ROI_width/2)
    if halfFlag == 0:
        upper_gamma = max([peak_fit[1][1][1], peak_fit[0][1][1]])
        lower_gamma = min([peak_fit[1][1][1], peak_fit[0][1][1]])
        gamma_minus_gamma = upper_gamma - lower_gamma
        side_length = (np.sqrt(2)/(2+2*np.sqrt(2)))*gamma_minus_gamma
        x_minus_gamma = gamma_minus_gamma/2
        x_minus_k = side_length/2
        k_minus_gamma = gamma_minus_gamma/2 - side_length/2
        # Cast the values to integers because they are pixel addresses
        gamma_minus_gamma = int(round(gamma_minus_gamma))
        x_minus_gamma = int(round(x_minus_gamma))
        x_minus_k = int(round(x_minus_k))
        x_minus_gamma = int(round(x_minus_gamma))
        k_minus_gamma = int(round(k_minus_gamma))    
        # Store the values in a dictionary and save that dictionary
        recyprocalLatticePoints = {"gamma_minus_gamma":gamma_minus_gamma, 
                                   "x_minus_gamma":x_minus_gamma,
                                   "x_minus_k":x_minus_k, 
                                   "x_minus_gamma":x_minus_gamma,
                                   "k_minus_gamma":k_minus_gamma,
                                   "gaussian_fits": peak_fit} 
    elif halfFlag == 1:
        lower_gamma = min([peak_fit[1][1][1], peak_fit[0][1][1]])
        x = max([peak_fit[1][1][1], peak_fit[0][1][1]])
        x_minus_gamma = x - lower_gamma
        gamma_minus_gamma = 2 * x_minus_gamma
        side_length = (np.sqrt(2)/(2+2*np.sqrt(2)))*gamma_minus_gamma
        x_minus_k = side_length/2
        k_minus_gamma = gamma_minus_gamma/2 - side_length/2
        # Cast the values to integers because they are pixel addresses
        gamma_minus_gamma = int(round(gamma_minus_gamma))
        x_minus_gamma = int(round(x_minus_gamma))
        x_minus_k = int(round(x_minus_k))
        x_minus_gamma = int(round(x_minus_gamma))
        k_minus_gamma = int(round(k_minus_gamma))
        recyprocalLatticePoints = {"gamma_minus_gamma":gamma_minus_gamma, 
                                   "x_minus_gamma":x_minus_gamma,
                                   "x_minus_k":x_minus_k, 
                                   "x_minus_gamma":x_minus_gamma,
                                   "k_minus_gamma":k_minus_gamma,
                                   "gaussian_fits": peak_fit}
    else:
        raise Exception("Invalid Flag Value")
    return(recyprocalLatticePoints)
    
def alignCenterAve_NiRuns_20170501():
    # inputs
    save_address = r'C:\Users\tchase56\Documents\UED\Ni\data\20160629\short\\'
    run_numbers = [59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 78, 80, 81, 82, 84, 85, 86]
    load_addresses = []
    delayStart, delayStop, delayStep = [64.65, 65.25, 0.1]
    delayStart_2, delayStop_2, delayStep_2 = [64.67, 65.27, 0.1]
    
    # form array of load addresses
    for i in range(len(run_numbers)):
        temp = save_address + r'scan' + str(run_numbers[i]) + '\images-ANDOR1\\'
        load_addresses.append(temp)
        
    # initialize run arrays
    runs = []
    runs_aligned_centered = [] 
    # attain first averaged and averaged/aligned run to form ROIs
    runs.append(load_timescan(load_addresses[0], delayStart, delayStop, delayStep))
    temp, peak_region, background_region = load_timescan(load_addresses[0], delayStart, delayStop, delayStep, align = 1, center = 1)
    runs_aligned_centered.append(temp)
    print('scan ' + str(run_numbers[0]))
    
    # Load runs, average, center and average
    # Two loops are necessary because there was an adjustment to account for drift
    for i in range(1,3):
        print('scan ' + str(run_numbers[i]))
        runs.append(load_timescan(load_addresses[i], delayStart, delayStop, delayStep))
        temp, _, _ = load_timescan(load_addresses[i], delayStart, delayStop, delayStep, align = 1, center = 1, peak_region = peak_region, background_region = background_region)
        runs_aligned_centered.append(temp)   
    for i in range(3,len(run_numbers)):
        print('scan ' + str(run_numbers[i]))
        runs.append(load_timescan(load_addresses[i], delayStart_2, delayStop_2, delayStep_2))
        temp, _, _ = load_timescan(load_addresses[i], delayStart_2, delayStop_2, delayStep_2, align = 1, center = 1, peak_region = peak_region, background_region = background_region)
        runs_aligned_centered.append(temp)
        
    # Average all runs together and save
    runs_ave = np.average(runs, 0)
    runs_aligned_centered_ave = np.average(runs_aligned_centered, 0)
    np.save(save_address + 'average_noAlign.npy', runs_ave)
    np.save(save_address + 'average_align_center.npy', runs_aligned_centered_ave) 
    return(runs_ave, runs_aligned_centered_ave)
    
def alignCenterAveRuns(save_address, run_numbers, delayStart, delayStop, delayStep):
    # inputs
    load_addresses = []
    
    # form array of load addresses
    for i in range(len(run_numbers)):
        temp = save_address + r'scan' + str(run_numbers[i]) + '\images-ANDOR1\\'
        load_addresses.append(temp)
        
    # initialize run arrays
    runs = []
    runs_aligned_centered = [] 
    # attain first averaged and averaged/aligned run to form ROIs
    runs.append(load_timescan(load_addresses[0], delayStart, delayStop, delayStep))
    temp, peak_region, background_region = load_timescan(load_addresses[0], delayStart, delayStop, delayStep, align = 1, center = 1)
    runs_aligned_centered.append(temp)
    print('scan ' + str(run_numbers[0]))
    
    # Load runs, average, center and average
    # Two loops are necessary because there was an adjustment to account for drift
    for i in range(len(run_numbers)):
        print('scan ' + str(run_numbers[i]))
        runs.append(load_timescan(load_addresses[i], delayStart, delayStop, delayStep))
        temp, _, _ = load_timescan(load_addresses[i], delayStart, delayStop, delayStep, align = 1, center = 1, peak_region = peak_region, background_region = background_region)
        runs_aligned_centered.append(temp)   
        
    # Average all runs together and save
    runs_ave = np.average(runs, 0)
    runs_aligned_centered_ave = np.average(runs_aligned_centered, 0)
    np.save(save_address + 'average_noAlign.npy', runs_ave)
    np.save(save_address + 'average_align_center.npy', runs_aligned_centered_ave) 
    return(runs_ave, runs_aligned_centered_ave)
    
    
def backgroundNormExpansionSym(save_address, runs_ave, runs_aligned_centered_ave, times, radius_mask, subtract_background = False, normalize = False, expansion_interpolation = False, symmetrize_image = False):

    # Rotate the image
    rotated = rotate_square(runs_aligned_centered_ave)  
    
    if subtract_background:
        # Subtract detector background
        scans_sub, scans_rot_sub, background = subtract_detector_background(runs_ave, rotated, radius_mask = radius_mask)
        # Calculate intensities of initial and subtracted
        initial_intensity = screen_intensity(rotated)
        subtracted_intensity = screen_intensity(scans_rot_sub)
        # Compare these two intensities by overlapping them
        fig, ax1 = plt.subplots()
        ax1.plot(times, initial_intensity - np.average(initial_intensity), label = 'Original')
        ax1.plot(times, subtracted_intensity - np.average(subtracted_intensity), label = 'Background Subtracted')
        ax1.set_title('Background Subtraction Results Overlapped')
        ax1.set_ylabel('Intensity (AU)')
        ax1.set_xlabel('Time Delay (PS)')
        plt.legend(loc = 'lower right')
        plt.show()
        # Plot the detector background
        fig, ax1 = plt.subplots()
        ax1.plot(times, background)
        ax1.set_title('Background')
        ax1.set_ylabel('Intensity (AU)')
        ax1.set_xlabel('Time Delay (PS)')
        plt.show()
        rotated = scans_rot_sub
    
    if normalize:           
        # Normalize Screen Intensity
        detector_intensity, detector_intensity_norm, norm = norm_screen(rotated)
        fig, ax1 = plt.subplots()
        ax1.plot(times, detector_intensity, label = 'Orginal')
        ax1.plot(times, detector_intensity_norm, label = 'Normalized')
        ax1.set_title('Intensity Normalization')
        ax1.set_ylabel('Intensity (AU)')
        ax1.set_xlabel('Time Delay (PS)')
        plt.legend(loc = 'upper left')
        plt.show()
        rotated = norm
    
    # Account for peak expansion with interpolation and renormalize
    if expansion_interpolation:
        peak_distance, peak_distance_expanded, runs_expanded = expansion_correction(rotated)
        fig, ax1 = plt.subplots()    
        ax1.set_title('Expansion Interpolation Test')
        ax1.set_ylabel('Intensity (AU)')
        ax1.set_xlabel('Time Delay (PS)')
        ax1.plot(times, peak_distance, label = 'original')
        ax1.plot(times, peak_distance_expanded, label = 'expanded')
        plt.legend()
        plt.show()
        rotated = runs_expanded
        
    if symmetrize_image:
        # Symmetrixe by rotating 
        symmetric_0 = []
        for i in range(np.shape(rotated)[0]):
            symmetric_0.append(np.array(symmetrize(rotated[i], fold = 4)))
    
        # Symmetrize by flipping
        symmetric = []
        for i in range(np.shape(symmetric_0)[0]):
            symmetric.append(np.array(symmetrize_flip(symmetric_0[i])))
        rotated = symmetric
        
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 1 and symmetrize_image == 1:
        np.save(save_address + 'expanded_norm_subBack_symmetrized.npy', rotated)
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 1 and symmetrize_image == 1:
        np.save(save_address + 'norm_subBack_symmetrized.npy', rotated)
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 1 and symmetrize_image == 1:
        np.save(save_address + 'expanded_subBack_symmetrized.npy', rotated)
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 0 and symmetrize_image == 1:
        np.save(save_address + 'expanded_norm_symmetrized.npy', rotated)
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 1 and symmetrize_image == 0:
        np.save(save_address + 'expanded_norm_subBack.npy', rotated)
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 1 and symmetrize_image == 1:
        np.save(save_address + 'subBack_symmetrized.npy', rotated)        
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 0 and symmetrize_image == 1:
        np.save(save_address + 'norm_symmetrized.npy', rotated)
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 1 and symmetrize_image == 0:
        np.save(save_address + 'norm_subBack.npy', rotated)
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 0 and symmetrize_image == 1:
        np.save(save_address + 'expanded_symmetrized.npy', rotated)
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 1 and symmetrize_image == 0:
        np.save(save_address + 'expanded_subBack.npy', rotated)
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 0 and symmetrize_image == 0:
        np.save(save_address + 'expanded_norm.npy', rotated)
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 0 and symmetrize_image == 1:
        np.save(save_address + 'symmetrized.npy', rotated)
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 1 and symmetrize_image == 0:
        np.save(save_address + 'subBack.npy', rotated)
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 0 and symmetrize_image == 0:
        np.save(save_address + 'norm.npy', rotated)
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 0 and symmetrize_image == 0:
        np.save(save_address + 'expanded.npy', rotated)
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 0 and symmetrize_image == 0:
        np.save(save_address + 'rotated.npy', rotated)
        
    # Plot returned image
    plt.figure()
    plt.imshow(rotated[-1]-rotated[0]).set_clim(0,12)
    plt.show() 
        
    return(rotated)
    
def difference_image(pattern):
    pattern = list(pattern)
    reference = pattern[0]
    difference = []
    for i in range(1, np.shape(pattern)[0]):
        difference.append(pattern[i]-reference)
    return(difference)
 
def formDiffuseLinescan(times, runs, runs_sub, 
                        ROI_width = 64, halfFlag = 0, plot = False, 
                        message = 'Choose two peaks along the vertical direction to form linescan_y'):
        
    # Form region of interest from two peaks along the vertical
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(runs[0] ,2, 
        halfLength=(ROI_width/2), contrastFactor = 0.25, message = message)
    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(runs[0], peak_region, background_region, halfLength=ROI_width/2)
    min_y = min([peak_fit_pump[1][1][1], peak_fit_pump[0][1][1]])
    max_y = max([peak_fit_pump[1][1][1], peak_fit_pump[0][1][1]])
    ROI_x_min = peak_fit_pump[0][0][1]-30
    ROI_x_max = peak_fit_pump[0][0][1]+30
    ROI_y_min = min_y
    ROI_y_max = max_y
    height = ROI_y_max - ROI_y_min
    width = ROI_x_max - ROI_x_min
    ROI = [ROI_x_min, ROI_y_min, width, height]
    
    # Form linescan of difference image along y
    diffuse = []
    for i in range(0,np.shape(runs_sub)[0]):
        subimage = runs_sub[i][ROI_y_min:ROI_y_max,ROI_x_min:ROI_x_max]
        diffuse.append(np.sum(subimage,1)) 
        
    if plot:            
        # Plot diffuse linescan for each time delay along y
        plt.figure()
        for i in range(np.shape(runs_sub)[0]):
            plt.plot(diffuse[i], label = str(round(times[i+1],1)) + ' ps')
        plt.legend(loc='upper center')
        plt.title('Diffuse Linescan Along Vertical')
   
        # Mask out the center and edges of diffraction pattern
        mask = make_mask(runs_sub[0], 250, int(1024/2-66))
   
        # Plot the region of interest we are integrating to form the linescan
        # Region of interest for vertical
        plt.figure()
        ax1 = plt.subplot(111)
        ax1.imshow(mask*runs_sub[-1], origin = 'upper').set_clim(0, 12)
        plt.title('ROI', fontsize = 20)
        ax1.add_patch(
            patches.Rectangle((ROI_x_min, ROI_y_min), width, height, fill = False, edgecolor = 'white', lw = 5)
        )
        ax1.set_xlim(250, 1024-250)
        ax1.set_ylim(1024-250, 250)
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)  
        plt.show()
        
    return(diffuse, ROI)
    
# Used to normalize time scans with multiple values below time zero
def Scale(graph, timeStepForNorm):
    norm = np.average(np.array(graph[timeStepForNorm]))
    graph = graph / norm
    return (graph, norm) 
    
def debyeWaller(save_address, delay_stage, pumped, points, numPeaks=5, ROI_width = 64, gaussian = True, plot = False):
    # Set the flags from the number of peaks
    debye_flag_2 = 0
    debye_flag_3 = 0
    debye_flag_4 = 0
    debye_flag_5 = 0
    debye_flag_6 = 0        
    if numPeaks > 2:
        debye_flag_2 = 1
    if numPeaks > 3:
        debye_flag_3 = 1
    if numPeaks > 4:
        debye_flag_4 = 1
    if numPeaks > 5:
        debye_flag_5 = 1
    if numPeaks > 6:
        debye_flag_6 = 1
    pumped_ave = np.average(pumped, 0)
    # Define peak ROIs
    [peak_region, background_region] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25, message = 'Choose 1st innermost peaks for Debye-Waller')
    [peak_region_2, background_region_2] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25, message = 'Choose 2nd innermost peaks for Debye-Waller')
    if debye_flag_2 == 1:
        [peak_region_3, background_region_3] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25, message = 'Choose 3rd innermost peaks for Debye-Waller')
        if debye_flag_3 == 1:
            [peak_region_4, background_region_4] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,8, halfLength=ROI_width/2, contrastFactor = 0.25, message = 'Choose 4th innermost peaks for Debye-Waller')
            if debye_flag_4 == 1:
                [peak_region_5, background_region_5] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25, message = 'Choose 5th innermost peaks for Debye-Waller')
                if debye_flag_5 == 1:
                    [peak_region_6, background_region_6] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,4, halfLength=ROI_width/2, contrastFactor = 0.25, message = 'Choose 6th innermost peaks for Debye-Waller')
                    if debye_flag_6 == 1:
                        [peak_region_7, background_region_7] = Region_Of_Interest_2.GetRegionOfInterest(pumped_ave ,8, halfLength=ROI_width/2, contrastFactor = 0.25, message = 'Choose 7th innermost peaks for Debye-Waller')
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

    # Save peak intensity (both intensity of ROI and Gaussian fit amplitude), as well as peak expansion
    peak_intensity = np.zeros((length,np.shape(pumped)[0]))
    peak_intensity_gaussian = np.zeros((length, np.shape(pumped)[0]))
    peak_expansion = np.zeros((length, np.shape(pumped)[0]))
    debye_intensity = np.zeros((length,np.shape(pumped)[0]))
    debye_gaussian = np.zeros((length-1, np.shape(pumped)[0]))
    for i in range(0,np.shape(pumped)[0]):    
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region, background_region, halfLength=ROI_width/2)
        peak_intensity[0,i] = np.average(peak_intensity_pump) - np.average(background_intensity)
        peak_intensity_gaussian[0,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
        peak_expansion[0,i] = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])
        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_2, background_region_2, halfLength=ROI_width/2)
        peak_intensity[1,i] = np.average(peak_intensity_pump) - np.average(background_intensity)   
        peak_intensity_gaussian[1,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
        peak_expansion[1,i] = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])
        if debye_flag_2 == 1:
            [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_3, background_region_3, halfLength=ROI_width/2)
            peak_intensity[2,i] = np.average(peak_intensity_pump) - np.average(background_intensity)
            peak_intensity_gaussian[2,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
            peak_expansion[2,i] = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])
            if debye_flag_3 == 1:
                [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_4, background_region_4, halfLength=ROI_width/2)
                peak_intensity[3,i] = np.average(peak_intensity_pump) - np.average(background_intensity)  
                peak_intensity_gaussian[3,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
                peak_expansion[3,i] = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])
                if debye_flag_4 == 1:
                    [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_5, background_region_5, halfLength=ROI_width/2)
                    peak_intensity[4,i] = np.average(peak_intensity_pump) - np.average(background_intensity) 
                    peak_intensity_gaussian[4,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
                    peak_expansion[4,i] = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])
                    if debye_flag_5 == 1:
                        [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_6, background_region_6, halfLength=ROI_width/2)
                        peak_intensity[5,i] = np.average(peak_intensity_pump) - np.average(background_intensity) 
                        peak_intensity_gaussian[5,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
                        peak_expansion[5,i] = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])
                        if debye_flag_6 == 1:
                            [peak_fit_pump, peak_intensity_pump, background_intensity] = Fit_Peaks_2.FitPeaks(pumped[i], peak_region_7, background_region_7, halfLength=ROI_width/2)
                            peak_intensity[6,i] = np.average(peak_intensity_pump) - np.average(background_intensity) 
                            peak_intensity_gaussian[6,i] = (np.average(np.array(peak_fit_pump)[:,0,0]) + np.average(np.array(peak_fit_pump)[:,1,0]))/2
                            peak_expansion[6,i] = peak_expansion_f(np.array(peak_fit_pump)[:,0,1], np.array(peak_fit_pump)[:,1,1])

    
    # Normalize peaks by innermost peak and scale them to 1 (for gaussian intensity)
    new_debye_gaussian_1 = peak_intensity_gaussian[1,:]/peak_intensity_gaussian[0,:]
    debye_gaussian[0,:], _ = Scale(new_debye_gaussian_1, list(range(points)))
    if debye_flag_2 == 1:
        new_debye_gaussian_2 = peak_intensity_gaussian[2,:]/peak_intensity_gaussian[0,:]
        debye_gaussian[1,:], _ = Scale(new_debye_gaussian_2, list(range(points)))
        if debye_flag_3 == 1:
            new_debye_gaussian_3 = peak_intensity_gaussian[3,:]/peak_intensity_gaussian[0,:]
            debye_gaussian[2,:], _ = Scale(new_debye_gaussian_3, list(range(points)))
            if debye_flag_4 == 1:
                new_debye_gaussian_4 = peak_intensity_gaussian[4,:]/peak_intensity_gaussian[0,:]
                debye_gaussian[3,:], _ = Scale(new_debye_gaussian_4, list(range(points)))
                if debye_flag_5 == 1:
                    new_debye_gaussian_5 = peak_intensity_gaussian[5,:]/peak_intensity_gaussian[0,:]
                    debye_gaussian[4,:], _ = Scale(new_debye_gaussian_5, list(range(points)))
                    if debye_flag_6 == 1:
                        new_debye_gaussian_6 = peak_intensity_gaussian[6,:]/peak_intensity_gaussian[0,:]
                        debye_gaussian[5,:], _ = Scale(new_debye_gaussian_6, list(range(points)))
                    
                   
    # Normalize peaks by innermost peak and scale them to 1 (for sum of ROI intensity)
    new_debye_intensity_1 = peak_intensity[1,:]/peak_intensity[0,:]
    debye_intensity[0,:], _ = Scale(new_debye_intensity_1, list(range(points)))
    if debye_flag_2 == 1:
        new_debye_intensity_2 = peak_intensity[2,:]/peak_intensity[0,:]
        debye_intensity[1,:], _ = Scale(new_debye_intensity_2, list(range(points)))
        if debye_flag_3 == 1:
            new_debye_intensity_3 = peak_intensity[3,:]/peak_intensity[0,:]
            debye_intensity[2,:], _ = Scale(new_debye_intensity_3, list(range(points)))
            if debye_flag_4 == 1:
                new_debye_intensity_4 = peak_intensity[4,:]/peak_intensity[0,:]
                debye_intensity[3,:], _ = Scale(new_debye_intensity_4, list(range(points)))
                if debye_flag_5 == 1:
                    new_debye_intensity_5 = peak_intensity[5,:]/peak_intensity[0,:]
                    debye_intensity[4,:], _ = Scale(new_debye_intensity_5, list(range(points)))
                    if debye_flag_6 == 1:
                        new_debye_intensity_6 = peak_intensity[6,:]/peak_intensity[0,:]
                        debye_intensity[5,:], _ = Scale(new_debye_intensity_6, list(range(points)))
 
    static_q = []
    # Normalize peak expansion by the first point
    static_q.append(peak_expansion[0,0]/2)
    peak_expansion_1 = peak_expansion[0,:]
    peak_expansion[0,:], _ = Scale(peak_expansion_1, list(range(points)))
    static_q.append(peak_expansion[1,0]/2)
    peak_expansion_2 = peak_expansion[1,:]
    peak_expansion[1,:], _ = Scale(peak_expansion_2, list(range(points)))
    if debye_flag_2 == 1:
        static_q.append(peak_expansion[2,0]/2)
        peak_expansion_3 = peak_expansion[2,:]
        peak_expansion[2,:], _ = Scale(peak_expansion_3, list(range(points)))
        if debye_flag_3 == 1:
            static_q.append(peak_expansion[3,0]/2)
            peak_expansion_4 = peak_expansion[3,:]
            peak_expansion[3,:], _ = Scale(peak_expansion_4, list(range(points)))
            if debye_flag_4 == 1:
                static_q.append(peak_expansion[4,0]/2)
                peak_expansion_5 = peak_expansion[4,:]
                peak_expansion[4,:], _ = Scale(peak_expansion_5, list(range(points)))
                if debye_flag_5 == 1:
                    static_q.append(peak_expansion[5,0]/2)
                    peak_expansion_6 = peak_expansion[5,:]
                    peak_expansion[5,:], _ = Scale(peak_expansion_6, list(range(points)))
                    if debye_flag_6 == 1:
                        static_q.append(peak_expansion[6,0]/2)
                        peak_expansion_7 = peak_expansion[6,:]
                        peak_expansion[6,:], _ = Scale(peak_expansion_7, list(range(points))) 
    np.save(save_address + "static_q", static_q)
    np.save(save_address + "debye_intensity", debye_intensity)
    np.save(save_address + "debye_intensity_gaussian", debye_gaussian)
    np.save(save_address + "peak_expansion", peak_expansion)   
                    
    if gaussian:
        debye = debye_gaussian
    else:
        debye = debye_intensity
        
    if plot:
        plt.figure()
        plt.plot(delay_stage, debye_intensity[0,:], label = "SecondOrder/FirstOrder")
        plt.title('Debye-Waller')
        plt.xlabel('Delay Stage')
        plt.show()

        if debye_flag_2 == 1:
            plt.plot(delay_stage, debye_intensity[2,:], label = "ThirdOrder/FirstOrder")
            plt.legend()
            plt.show()

            if debye_flag_3 == 1:
                plt.plot(delay_stage, debye_intensity[3,:], label = "FourthOrder/FirstOrder")
                plt.legend()
                plt.show()
        
                if debye_flag_4 == 1:
                    plt.plot(delay_stage, debye_intensity[4,:], label = "FifthOrder/FirstOrder")
                    plt.legend()
                    plt.show()
        
                    if debye_flag_5 == 1:
                        plt.plot(delay_stage, debye_intensity[5,:], label = "SixthOrder/FirstOrder")
                        #plt.legend(loc = 'lower left')
                        plt.legend()
                        plt.show()
                
                        if debye_flag_6 == 1:
                            plt.plot(delay_stage, debye_intensity[6,:], label = "SeventhOrder/FirstOrder")
                            plt.legend()
                            plt.show()   
                            
        plt.figure()
        plt.plot(delay_stage, peak_expansion[0,:], label = 'FirstOrder Expansion')
        plt.legend()
        plt.title('Lattice Expansion')
        plt.xlabel('Delay Stage')
        plt.show()
        plt.plot(delay_stage, peak_expansion[1,:], label = 'SecondOrder Expansion')
        plt.legend()
        plt.show()
        if debye_flag_2 == 1:
            plt.plot(delay_stage, peak_expansion[2,:], label = 'ThirdOrder Expansion')
            plt.legend()
            plt.show()

            if debye_flag_3 == 1:
                plt.plot(delay_stage, peak_expansion[3,:], label = 'FourthOrder Expansion')
                plt.legend()
                plt.show()
        
                if debye_flag_4 == 1:
                    plt.plot(delay_stage, peak_expansion[4,:], label = 'FifthOrder Expansion')
                    plt.legend()
                    plt.show()
        
                    if debye_flag_5 == 1:
                        plt.plot(delay_stage, peak_expansion[5,:], label = 'SixthOrder Expansion')
                        #plt.legend(loc = 'lower left')
                        plt.legend()
                        plt.show()
                
                        if debye_flag_6 == 1:
                            plt.plot(delay_stage, peak_expansion[6,:], label = 'SeventhOrder Expansion')
                            plt.legend()
                            plt.show() 
    return(debye, peak_expansion, static_q)

'''
# Debye-Waller intensity fits no offset
def exp_debye(t, tau, A, q):
    return (1-A)*np.exp( ( (-1*t+0) *q**2 ) /tau) +A 
    
def exp_debye_fit(x, data, tau, A, q):

    def dif(lltau, llA):
        dif = data - exp_debye(x, lltau, llA, q)
        return(dif)
    
    def sumsq(ltau, lA):
        err = dif(ltau, lA)
        return np.sum(err*err)

    res = opt.minimize(sumsq, x0 = (tau, A))
    return res
'''

# Debye-Waller intensity fits no offset
def exp_debye(t, tau, A, q):
    sol = []
    for i in t:
        if i<0:
            sol.append(1)
        else:
            sol.append( (1-A)*np.exp( ( (-1*i+0) *q**2) /tau) + A )
    return sol  
def exp_debye_fit(x, data, tau_seed_debye, A_seed_debye, q):
    exp_debye_2 = lambda lt, ltau, lA: exp_debye(lt, ltau, lA, q)
    tau = tau_seed_debye
    A = A_seed_debye      
        
    popt, pcov = opt.curve_fit(exp_debye_2, x, data, p0 = (tau, A))
    return popt, pcov

# Diffuse intensity Fits no offset
def exp_diffuse(t, tau, A):
    return A * (1-np.exp((-1*t+0)/tau)) * (0.5) * (np.sign(t) + 1)
def exp_diffuse_fit(x, data, tau_seed_diffuse, A_seed_diffuse):
    tau = tau_seed_diffuse
    A = A_seed_diffuse
    popt, pcov = opt.curve_fit(exp_diffuse, x, data, p0 = (tau, A))
    return popt, pcov
    
# Thermalization analysis of diffuse linescan
def thermalization(times, lineScan_short, times_timescan, recyp_lat, timeZeroShift, 
                   peak_intensity_norm, lower, upper, length, debye_x_lowLim, debye_x_upLim, 
                   debye_y_lowLim, debye_y_upLim, fitInten_y_lowLim, fitInten_y_upLim, 
                   diffuse_y_lowLim, diffuse_y_upLim, tau_y_lowLim, tau_y_upLim, tau_seed_debye,
                   A_seed_debye, tau_seed_diffuse, A_seed_diffuse, debye_region_bottom, 
                   diffuse_region_bottom, highlights_debye, highlights_diffuse, static_q, scale_diffuse = 0, 
                   ROI_title = 'ROI', plot = False, start = 'center', half = False, scaling_factors = False):
    # Shift times to align time zero
    times = times - timeZeroShift
    times = times[1:]
    
    if start == 'center':
        # Sum over regions determined in previous step making them symmetric
        middle = (upper+lower)/2
        bins = []
        bins.append(int(lower))
        temp = middle - length/2
        while lower < temp:
            bins.insert(1, int(temp))
            temp = temp - length
        temp = middle + length/2
        while temp < upper:
            bins.append(int(temp))
            temp = temp + length
        bins.append(int(upper))
        sums = []
        for i in range(np.shape(np.array(lineScan_short))[0]):
            sums.append([])
            for j in range(len(bins)-1):
                sums[i].append(np.sum(lineScan_short[i][bins[j]:bins[j+1]]))
        sums = np.array(sums)
    elif start == 'right':
        # sum over regions starting at right hand side and moving left
        bins = []
        bins.append(lower)
        temp = upper
        while lower < temp:
            bins.insert(1, int(temp))
            temp = temp - length
        sums = []
        for i in range(np.shape(np.array(lineScan_short))[0]):
            sums.append([])
            for j in range(len(bins)-1):
                sums[i].append(np.sum(lineScan_short[i][bins[j]:bins[j+1]]))
        sums = np.array(sums)
    elif start == 'left':
        # sum over regions starting at left hand side and moving right
        bins = []
        temp = lower
        while upper > temp:
            bins.append(int(temp))
            temp = temp + length
        bins.append(upper)
        sums = []
        for i in range(np.shape(np.array(lineScan_short))[0]):
            sums.append([])
            for j in range(len(bins)-1):
                sums[i].append(np.sum(lineScan_short[i][bins[j]:bins[j+1]]))
        sums = np.array(sums)
    else:
        i = 0 
        bins = []
        temp = lower
        while upper > temp:
            bins.append(int(temp))
            try:
                temp = temp + start[i]
            except:
                temp = temp + length
            i += 1
        bins.append(upper)
        sums = []
        for i in range(np.shape(np.array(lineScan_short))[0]):
            sums.append([])
            for j in range(len(bins)-1):
                sums[i].append(np.sum(lineScan_short[i][bins[j]:bins[j+1]]))
        sums = np.array(sums)
        if scaling_factors:
            sums = sums/np.array(scaling_factors)
        
    # Scale maximums of diffuse thermalizations to compare
    scale_const = []
    if scale_diffuse:
        for i in range(np.shape(sums)[1]):
            sums[:,i],scale_elem = Scale(sums[:,i], [-1,-2,-3,-4])
            scale_const.append(scale_elem)
                        
    # Fit the debye-waller data
    fits_debye = []
    error_debye = []
    #debye_region_top = -1
    debye_region_top = np.shape(peak_intensity_norm)[1]
    for i in range(np.shape(peak_intensity_norm)[0]):
        fits_debye.append(exp_debye_fit(times_timescan[debye_region_bottom:debye_region_top], 
                                        peak_intensity_norm[i,debye_region_bottom:debye_region_top],
                                        tau_seed_debye, A_seed_debye, static_q[i+1])[0])
        error_debye.append(np.sqrt(exp_debye_fit(times_timescan[debye_region_bottom:debye_region_top], 
                                                 peak_intensity_norm[i,debye_region_bottom:debye_region_top],
                                                 tau_seed_debye, A_seed_debye, static_q[i+1])[1][0,0]))

    # Fit the diffuse data
    fits_diffuse = []
    error_diffuse = []

    #diffuse_region_top = -1
    diffuse_region_top = np.shape(sums)[0]
    for i in range(np.shape(sums)[1]):
        fits_diffuse.append(exp_diffuse_fit(times[diffuse_region_bottom:diffuse_region_top], 
                                            sums[diffuse_region_bottom:diffuse_region_top,i],
                                            tau_seed_diffuse, A_seed_diffuse)[0])
        error_diffuse.append(np.sqrt(exp_diffuse_fit(times[diffuse_region_bottom:diffuse_region_top], 
                                                     sums[diffuse_region_bottom:diffuse_region_top,i],
                                                     tau_seed_diffuse, A_seed_diffuse)[1]))
    
    # Calculate pixel values for positions of time constant calculations                                                
    pixel_timeConst = list(bins)
    for i in range(len(bins)-1):
        pixel_timeConst[i] = round((pixel_timeConst[i] + pixel_timeConst[i+1])/2)
    del pixel_timeConst[-1]
    bins = np.array(bins)
    pixel_timeConst = np.array(pixel_timeConst)
    percentGamma_timeConst = (pixel_timeConst)/recyp_lat["gamma_minus_gamma"]
    
    if plot:
        # Plot figure similar to paper
        f, (ax1, ax2) = plt.subplots(2, 1, sharex = True)

        # Subplot with debye-waller
        colors_0 = ['c', 'b', 'r', 'k', 'g']
        label_debye = ["220", "400", "420", "440", "600"]
        #for i in range(np.shape(peak_intensity_norm)[0]):
        for i in highlights_debye:
                ax1.scatter(times_timescan, peak_intensity_norm[i], c= colors_0[i], label = label_debye[i])

        # Plot vertical lines showing fit region
        debye_region_bottom_ps = times_timescan[debye_region_bottom]
        debye_region_top_ps = times_timescan[debye_region_top-1]
        ax1.axvline(debye_region_bottom_ps, linestyle = ':', c = 'k')
        ax1.axvline(debye_region_top_ps, linestyle = ':', c = 'k')
        ax1.set_xlim(-23,101)
        ax1.set_ylabel("Intensity Bragg Peak (A.U.)")
        ax1.set_title("Timescale Comparison\n" + ROI_title)
        ax1.legend(loc="upper right")


        # debye-waller fits
        #for i in range(np.shape(peak_intensity_norm)[0]):
        for i in highlights_debye:
            ax1.plot(np.arange(debye_region_bottom_ps-1, debye_region_top_ps, 0.05), exp_debye(np.arange(debye_region_bottom_ps-1, debye_region_top_ps, 0.05), *fits_debye[i], static_q[i+1]), c = colors_0[i])   
        ax1.set_xlim(debye_x_lowLim, debye_x_upLim)
        ax1.set_ylim(debye_y_lowLim, debye_y_upLim)
        # Subplot with diffuse regions
        #colors = ['y', 'm', 'c', 'r', 'g', 'b']

        label = ["$\Gamma$", "$X$", '$\Gamma_2$']
        col = ["b", "g", "r"]
        #for i in range(np.shape(sums)[1]):
        for j,i in enumerate(highlights_diffuse):
            try:
                ax2.scatter(times, sums[:,i], label = label[j] + " / " + str( int(round(scale_const[i])) ) , c=col[j] )
            except:
                ax2.scatter(times, sums[:,i], label = label[j] , c=col[j] )
                
    
        # Plot vertical line showing fit region
        diffuse_region_bottom_ps = times[diffuse_region_bottom]
        diffuse_region_top_ps = times[diffuse_region_top-1]
        ax2.axvline(diffuse_region_bottom_ps, linestyle = ':', c = 'k')
        ax2.axvline(diffuse_region_top_ps, linestyle = ':', c = 'k')
        if scale_diffuse:
            ax2.set_ylim(-0.2, 1.5)
        else:
            ax2.set_ylim(fitInten_y_lowLim, fitInten_y_upLim)
        ax2.set_ylabel("Intensity Diffuse")
        ax2.set_xlabel("Time Delay (PS)")
        ax2.legend(loc = "lower right")


        # diffuse fits
        #for i in range(np.shape(sums)[1]):
        for i in highlights_diffuse:   
            ax2.plot(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.05), exp_diffuse(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.05), *fits_diffuse[i]))
        plt.show()   


        plt.figure()
        ax1 = plt.subplot(111)
        for i in range(np.shape(np.array(lineScan_short))[0]):
            ax1.plot(lineScan_short[i], label = str(round(times[i], 1)) + ' ps')    
        ax1.set_ylim(diffuse_y_lowLim, diffuse_y_upLim)
        #ax1.legend(loc = 'upper center', prop = {'size':12})
        ax1.set_xlabel('Pixel (' + str(length) +  ' Pixels Wide For Time Constants)')
        ax1.set_ylabel('Intensity (A.U)')
        ax1.set_title(ROI_title)
        ax1.legend(bbox_to_anchor = [1.17,0], loc = 'lower left', borderaxespad = 0, title = "Scattering Intensity")
        plt.subplots_adjust(right = 0.67)

        ax2 = ax1.twinx()
        ax2.errorbar(pixel_timeConst, np.array(fits_diffuse)[:,0], yerr = np.array(error_diffuse)[:,0,0], fmt = 'o-', label = "Time Const")
        ax2.legend(bbox_to_anchor = [1.13,1], loc = 'upper left', borderaxespad = 0)
        for i in bins:
            ax2.axvline(i, linestyle = ':', c='k')
        for i in highlights_diffuse:
            ax2.axvline(pixel_timeConst[i], linestyle = '-', c = 'r')

        ax2.set_ylim(tau_y_lowLim,tau_y_upLim)
        ax2.set_ylabel('Time Constant (PS)')
        ticks = [0, recyp_lat["k_minus_gamma"], recyp_lat["x_minus_gamma"], recyp_lat["x_minus_gamma"]+recyp_lat["x_minus_k"], recyp_lat["gamma_minus_gamma"]]
        labels = ['$\Gamma$', 'K', 'X', 'K', '$\Gamma$']
        plt.xticks(ticks, labels)
        if not half:
            ax2.set_xlim(0,recyp_lat["gamma_minus_gamma"])
        else:
            ax2.set_xlim(0,recyp_lat["x_minus_gamma"])
        plt.show()
        
    return(fits_debye, error_debye, fits_diffuse, error_diffuse, percentGamma_timeConst, bins, sums, lineScan_short, scale_const)
    
# Plot thermalization comparison 
    
def thermPlot(times, recyp_lat, diffuse_region_bottom,
              lineScan_short, tau_diffuse, error_diffuse, gamma_to_gamma_percent, sums, highlights_diffuse, bins, 
              length, diffuse_y_lowLim, diffuse_y_upLim, fitInten_y_lowLim, fitInten_y_upLim, tau_y_lowLim, tau_y_upLim, highlights_linescan, scale_const, half = False,
              flag_2 = 0, lineScan_short_2 = 0, tau_diffuse_2 = 0, error_diffuse_2 = 0, gamma_to_gamma_percent_2 = 0, sums_2 = 0, highlights_diffuse_2 = 0, bins_2 = 0, 
              length_2 = 0, diffuse_y_lowLim_2 = 0, diffuse_y_upLim_2 = 0, fitInten_y_lowLim_2 = 0, fitInten_y_upLim_2 = 0, tau_y_lowLim_2 = 0, tau_y_upLim_2 = 0, highlights_linescan_2 = 0, scale_const_2 = 0, half_2 = False,
              scale_diffuse=0, ROI_title = 'ROI', label = ["$\Gamma$", "$X$", "$\Gamma_2$"]):  
    # Plot figure similar to paper
    f, (ax1) = plt.subplots()

    # Necessary info 
    col = ["b", "g", "r"]
    diffuse_region_top = np.shape(sums)[0]
    times = times[1:]
    
    # Pixel values associated with time constants
    pixel_timeConst = list(bins)
    for i in range(len(bins)-1):
        pixel_timeConst[i] = round((pixel_timeConst[i] + pixel_timeConst[i+1])/2)
    del pixel_timeConst[-1]
    bins = np.array(bins)
    pixel_timeConst = np.array(pixel_timeConst)
    if flag_2:
        pixel_timeConst_2 = list(bins_2)
        for i in range(len(bins_2)-1):
            pixel_timeConst_2[i] = round((pixel_timeConst_2[i] + pixel_timeConst_2[i+1])/2)
        del pixel_timeConst_2[-1]
        bins_2 = np.array(bins_2)
        pixel_timeConst_2 = np.array(pixel_timeConst_2)

    
    #for i in range(np.shape(sums)[1]):
    j = 0
    for i in highlights_diffuse:
        try:
            ax1.plot(times, sums[:,i], 'o', label = '(' + label[j] + ')' + " / " + str( int(round(scale_const[i])) ), c = col[j])
        except:
            ax1.plot(times, sums[:,i], 'o', label = label[j] , c=col[j])
        j+=1                          

    if flag_2:
        for i in highlights_diffuse_2:
            try:
                ax1.plot(times, sums_2[:,i], 'o', label = '(' + label[j] + ')' + " / " + str( int(round(scale_const_2[i])) ), c = col[j])  
            except:
                ax1.plot(times, sums_2[:,i], 'o', label = label[j] , c=col[j])  
            j+=1                        
    
    # Plot vertical line showing fit region
    diffuse_region_bottom_ps = times[diffuse_region_bottom]
    diffuse_region_top_ps = times[diffuse_region_top-1]
    ax1.axvline(diffuse_region_bottom_ps, linestyle = ':', c = 'k')
    ax1.axvline(diffuse_region_top_ps, linestyle = ':', c = 'k')
    if scale_diffuse:
        ax1.set_ylim(-0.1, 1.5)
    else:
        ax1.set_ylim(fitInten_y_lowLim, fitInten_y_upLim)
    ax1.set_xlim(-1, 16.5)
    ax1.set_ylabel("Intensity Diffuse")
    ax1.set_xlabel("Time Delay (PS)")
    ax1.legend(loc = "lower right")
    ax1.set_title(ROI_title + ' Fits')


    # diffuse fits
    #for i in range(np.shape(sums)[1]):
    for i in highlights_diffuse:   
        ax1.plot(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.05), exp_diffuse(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.05), *tau_diffuse[i]))
    if flag_2:
        for i in highlights_diffuse_2:   
            ax1.plot(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.05), exp_diffuse(np.arange(diffuse_region_bottom_ps-1, diffuse_region_top_ps, 0.05), *tau_diffuse_2[i]))
    plt.show()   


    plt.figure()    
    ax1 = plt.subplot(111)
    for i in highlights_linescan:
        ax1.plot(lineScan_short[i], label = str(round(times[i], 1)) + ' ps')    
    ax1.set_ylim(diffuse_y_lowLim, diffuse_y_upLim)
    #ax1.legend(loc = 'upper center', prop = {'size':12})
    ax1.set_xlabel('Pixel (' + str(length) +  ' Pixels Wide For Time Constants)')
    ax1.set_ylabel('Intensity (A.U)')
    ax1.set_title(ROI_title)
    ax1.legend(bbox_to_anchor = [1.17,0], loc = 'lower left', borderaxespad = 0, title = "Scattering Intensity")
    plt.subplots_adjust(right = 0.67)

    ax2 = ax1.twinx()
    ax2.errorbar(pixel_timeConst, np.array(tau_diffuse)[:,0], yerr = np.array(error_diffuse)[:,0,0], fmt = 'o-', label = "Time Const")
    ax2.legend(bbox_to_anchor = [1.13,1], loc = 'upper left', borderaxespad = 0)
    for i in highlights_diffuse:
        ax2.axvline(pixel_timeConst[i], linestyle = '-', c = 'r')        

    ax2.set_ylim(tau_y_lowLim,tau_y_upLim)
    ax2.set_ylabel('Time Constant (PS)')
    ticks = [0, recyp_lat["k_minus_gamma"], recyp_lat["x_minus_gamma"], recyp_lat["x_minus_gamma"]+recyp_lat["x_minus_k"], recyp_lat["gamma_minus_gamma"]]
    labels = ['$\Gamma$', 'K', 'X', 'K', '$\Gamma$']
    plt.xticks(ticks, labels)
    if not half:
        ax2.set_xlim(0,recyp_lat["gamma_minus_gamma"])
    else:
        ax2.set_xlim(0,recyp_lat["x_minus_gamma"])
    plt.show()    
    
    if flag_2:
        plt.figure()    
        ax1 = plt.subplot(111)
        for i in highlights_linescan_2:
            ax1.plot(lineScan_short_2[i], label = str(round(times[i], 1)) + ' ps')    
        ax1.set_ylim(diffuse_y_lowLim_2, diffuse_y_upLim_2)
        #ax1.legend(loc = 'upper center', prop = {'size':12})
        ax1.set_xlabel('Pixel (' + str(length_2) +  ' Pixels Wide For Time Constants)')
        ax1.set_ylabel('Intensity (A.U)')
        ax1.set_title(ROI_title)
        ax1.legend(bbox_to_anchor = [1.17,0], loc = 'lower left', borderaxespad = 0, title = "Scattering Intensity")
        plt.subplots_adjust(right = 0.67)

        ax2 = ax1.twinx()
        ax2.errorbar(pixel_timeConst_2, np.array(tau_diffuse_2)[:,0], yerr = np.array(error_diffuse_2)[:,0,0], fmt = 'o-', label = "Time Const")
        ax2.legend(bbox_to_anchor = [1.13,1], loc = 'upper left', borderaxespad = 0)
        for i in highlights_diffuse_2:
            ax2.axvline(pixel_timeConst_2[i], linestyle = '-', c = 'r')

        ax2.set_ylim(tau_y_lowLim_2,tau_y_upLim_2)
        ax2.set_ylabel('Time Constant (PS)')
        ticks = [0, recyp_lat["k_minus_gamma"], recyp_lat["x_minus_gamma"], recyp_lat["x_minus_gamma"]+recyp_lat["x_minus_k"], recyp_lat["gamma_minus_gamma"]]
        labels = ['$\Gamma$', 'K', 'X', 'K', '$\Gamma$']
        plt.xticks(ticks, labels)
        if not half_2:
            ax2.set_xlim(0,recyp_lat["gamma_minus_gamma"])
        else:
            ax2.set_xlim(0,recyp_lat["x_minus_gamma"])
        plt.show()  
        
def debPlot(times_timescan, peak_intensity_norm, highlights_debye, 
            fits_debye, tau_error_debye, debye_region_bottom, 
            debye_x_lowLim, debye_x_upLim, debye_y_lowLim, debye_y_upLim, static_q):
    # Plot debye-waller and fitted exponentials
    f, (ax1) = plt.subplots()

    # Subplot with debye-waller
    colors_0 = ['c', 'b', 'r', 'k', 'g']
    label_debye = ["220", "400", "420", "440", "600"]
    #for i in range(np.shape(peak_intensity_norm)[0]):
    for i in highlights_debye:
        ax1.scatter(times_timescan, peak_intensity_norm[i], c= colors_0[i], label = label_debye[i]+'/200')

    # Plot vertical lines showing fit region
    debye_region_top = np.shape(peak_intensity_norm)[1]
    debye_region_bottom_ps = times_timescan[debye_region_bottom]
    debye_region_top_ps = times_timescan[debye_region_top-1]
    ax1.axvline(debye_region_bottom_ps, linestyle = ':', c = 'k')
    ax1.axvline(debye_region_top_ps, linestyle = ':', c = 'k')
    ax1.set_xlim(-23,101)
    ax1.set_ylabel("Intensity Bragg Peak (A.U.)")
    ax1.set_title("Timescale Comparison\nDebye-Waller")
    ax1.set_xlabel("Time Delay (ps)")
    ax1.legend(loc="upper right")

    # debye-waller fits
    #for i in range(np.shape(peak_intensity_norm)[0]):
    for i in highlights_debye:
        ax1.plot(np.arange(debye_region_bottom_ps-1, debye_region_top_ps, 0.05), 
                 exp_debye(np.arange(debye_region_bottom_ps-1, debye_region_top_ps, 0.05), *fits_debye[i], static_q[i+1]), c = colors_0[i])   
    ax1.set_xlim(debye_x_lowLim, debye_x_upLim)
    ax1.set_ylim(debye_y_lowLim, debye_y_upLim)
    
    # Plot time constants of debye-waller and their errors
    fits_debye = np.array(fits_debye)    
    fits = []
    error = []
    label = []
    for i in highlights_debye:
        fits.append(fits_debye[i,0])
        error.append(tau_error_debye[i])
        label.append(label_debye[i])
                
    f, ax1 = plt.subplots()
    ax1.errorbar(range(np.shape(highlights_debye)[0]), fits, yerr = error, fmt = 'o-')
    plt.xticks(range(np.shape(highlights_debye)[0]), label)
    ax1.set_xlim(-1, np.shape(highlights_debye)[0])
    ax1.set_title('Debye-Waller Time Constants')
    ax1.set_xlabel('Diffraction Peak')
    ax1.set_ylabel('Time Constant')
    
# Ni dispersion plots    
def dispPlot(times, recyp_lat, diffuse_region_bottom,
              lineScan_short, dispersion,
              diffuse_y_lowLim, diffuse_y_upLim, fitInten_y_lowLim, fitInten_y_upLim, disp_y_lowLim, disp_y_upLim, highlights_linescan, half = False,
              ROI_title = 'ROI', interpolate_disp = False):  
                  
    # Must manually make copy of the input dictionary otherwise the dictionary
    # will be altered in the main file
    dispersion = copy.deepcopy(dispersion)
    times = times[1:]
     
    if interpolate_disp:
        # implement cubic spline 
        f_T1 = interpolate.interp1d(dispersion['T1']['X'], dispersion['T1']['Y'], kind = 'linear')
        f_T2 = interpolate.interp1d(dispersion['T2']['X'], dispersion['T2']['Y'], kind = 'linear')
        f_L = interpolate.interp1d(dispersion['L']['X'], dispersion['L']['Y'], kind = 'linear')
        # Populate the dictionary with interpolated values
        x = np.arange(0, 120.5, 0.01)
        dispersion['T1']['X'] = list(x)
        dispersion['T1']['Y'] = list(f_T1(x)) 
        dispersion['T2']['X'] = list(x)
        dispersion['T2']['Y'] = list(f_T2(x)) 
        dispersion['L']['X'] = list(x)
        dispersion['L']['Y'] = list(f_L(x))
    
    # Create dispersion to plot    
    if half:
        T1_x = dispersion['T1']['X']
        T1_y = dispersion['T1']['Y']
        T2_x = dispersion['T2']['X']
        T2_y = dispersion['T2']['Y']
        L_x = dispersion['L']['X']
        L_y = dispersion['L']['Y']
    else:
        temp = dispersion['T1']['X']    
        temp = temp[::-1]
        temp = np.array(temp)
        temp = -1*(temp - temp[0])
        temp = temp + 121
        T1_x = dispersion['T1']['X']
        T1_x.extend(temp)
 
        temp = dispersion['T1']['Y']    
        temp = temp[::-1]
        T1_y = dispersion['T1']['Y']
        T1_y.extend(temp)        
        
        temp = dispersion['T2']['X']    
        temp = temp[::-1]
        temp = np.array(temp)
        temp = -1*(temp - temp[0])
        temp = temp + 121
        T2_x = dispersion['T2']['X']
        T2_x.extend(temp)
        
        temp = dispersion['T2']['Y']    
        temp = temp[::-1]        
        T2_y = dispersion['T2']['Y']
        T2_y.extend(temp)        
        
        temp = dispersion['L']['X']    
        temp = temp[::-1]
        temp = np.array(temp)
        temp = -1*(temp - temp[0])
        temp = temp + 121
        L_x = dispersion['L']['X']
        L_x.extend(temp)

        temp = dispersion['L']['Y']    
        temp = temp[::-1]        
        L_y = dispersion['L']['Y']
        L_y.extend(temp)
        
    f, (ax1, ax2) = plt.subplots(2,1, sharex = True)    

    plt.subplots_adjust(right = 0.7)

    ax1.plot(T1_x, T1_y, '--', label = "T1")
    ax1.plot(L_x, L_y, '--', label = "L")
    ax1.plot(T2_x, T2_y, '--', label = "T2")

    ax1.legend(bbox_to_anchor = [1.05,1], loc = 'upper left', borderaxespad = 0)
       

    ax1.set_ylim(disp_y_lowLim,disp_y_upLim)
    ax1.set_ylabel('$\\nu$ (THz)')
    ax1.set_title('Dispersion Branches')
    ticks = [0, recyp_lat["k_minus_gamma"], recyp_lat["x_minus_gamma"], recyp_lat["x_minus_gamma"]+recyp_lat["x_minus_k"], recyp_lat["gamma_minus_gamma"]]
    labels = ['$\Gamma$', 'K', 'X', 'K', '$\Gamma$']
    plt.xticks(ticks, labels)
    if not half:
        ax1.set_xlim(0,recyp_lat["gamma_minus_gamma"])
    else:
        ax1.set_xlim(0,recyp_lat["x_minus_gamma"])
    plt.show()    
    
    for i in highlights_linescan:
        ax2.plot(lineScan_short[i], label = str(round(times[i], 1)) + ' ps')    
    ax2.set_ylim(diffuse_y_lowLim, diffuse_y_upLim)
    #ax1.legend(loc = 'upper center', prop = {'size':12})
    ax2.set_ylabel('Intensity (A.U)')
    ax2.set_title(ROI_title)
    ax2.legend(bbox_to_anchor = [1.05,0], loc = 'lower left', borderaxespad = 0, title = "Scattering Intensity")
    

    
    
    #
  
    
    
    

    
    
    
    