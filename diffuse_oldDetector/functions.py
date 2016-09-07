'''
Tyler Chase
09/07/2016
- functions for symmetrizing diffraction pattern
- function for masking out edges and center of phosphor screen
- function for forming linescan and plotting ROI
- function for converting from delay stage to picoseconds
- function for expanding using the 2-D interpolation
- function for calculating the peak separation for a given order
'''





# Diffuse Tools
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
import matplotlib.patches as patches
from scipy import interpolate

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
    masked = image*mask
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
    
def ps_convert(start, stop, step, time_zero):
    # Delay Stage Settings
    delay_stage = np.arange(start, stop + (step)/2, step)
    speed_of_light = 299792458 # meters/second
    delay = (np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12 # Delay in ps
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