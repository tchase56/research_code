import numpy as np
import matplotlib.pyplot as plt 

''' Values to change for each run'''
load_address = r'E:\tchase56\TaS2\20160611\TimeScan\scan6\images-ANDOR1\\'
save_address = r'E:\tchase56\TaS2\20160611\TimeScan\scan6\images-ANDOR1\\'
# Delay stage settings
time_zero = 63.2

delayStage_start, delayStage_end, delayStage_step  = [62.90, 63.80, 0.0150] 

# Determine delay in ps from delay stage
delay_stage = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
speed_of_light = 299792458 # meters/second
delay = np.round((np.array(delay_stage) - time_zero)*2*10**(-3)*(1.0/speed_of_light)*10**12,2) # Delay in ps


# Load Mask and Images
boundary_mask = np.load(load_address + 'boundary_mask.npy')
peak_mask = np.load(load_address + 'peak_mask.npy')
peak_mask_inv = -1.0*(peak_mask-1.0)
throughbeam_mask = np.load(load_address + 'throughbeam_mask.npy')
background_mask = np.load(load_address + 'mask_detector_background.npy')
images = np.load(load_address + 'aligned_rotated.npy')

intensity = []
bragg_intensity = []
superstructure = []
background = []
bragg_peaks = []
throughbeam = []
throughbeam_intensity = []
background_intensity = []
for i in range(0,np.shape(images)[0]):
    superstructure.append(images[i][:,:] * boundary_mask*peak_mask)
    bragg_peaks.append(images[i][:,:] * boundary_mask*peak_mask_inv)
    throughbeam.append(images[i][:,:] * throughbeam_mask)
    background.append(images[i][:,:] * background_mask)
    intensity.append(np.sum(superstructure[i]))
    bragg_intensity.append(np.sum(bragg_peaks[i]))
    throughbeam_intensity.append(np.sum(throughbeam[i]))
    background_intensity.append(np.sum(background[i]))
    
# Define normalization
intensity = np.array(intensity)
bragg_intensity = np.array(bragg_intensity)
norm = intensity + bragg_intensity

  
plt.figure()
plt.subplot(121)
plt.title("superstructure intensity")
plt.plot(delay, intensity-background_intensity)
plt.xlabel("time(PS)")
plt.ylabel("Intensity(AU)")
plt.subplot(122)
plt.title("bragg intensity")
plt.plot(delay, bragg_intensity-background_intensity)
plt.xlabel("time(PS)")
plt.ylabel("Intensity(AU)")
plt.show()


plt.figure()
plt.subplot(222)
plt.title("superstructure intensity")
plt.plot(delay, intensity)
plt.xlabel("time(PS)")
plt.ylabel("Intensity(AU)")
plt.subplot(223)
plt.title("bragg intensity")
plt.plot(delay, bragg_intensity)
plt.xlabel("time(PS)")
plt.ylabel("Intensity(AU)")
plt.subplot(224)
plt.title("throughbeam intensity")
plt.plot(delay, throughbeam_intensity)
plt.xlabel("time(PS)")
plt.ylabel("Intensity(AU)")
plt.subplot(221)
plt.title("background intensity")
plt.plot(delay, background_intensity)
plt.xlabel("time(PS)")
plt.ylabel("Intensity(AU)")
plt.show()

plt.figure()
plt.subplot(221)
plt.title("background")
plt.imshow(images[i][:,:]*background_mask).set_clim(100,2000)
plt.subplot(222)
plt.title("superstructure intensity")
plt.imshow(images[i][:,:]*boundary_mask*peak_mask).set_clim(100,2000)
plt.subplot(223)
plt.title("bragg intensity")
plt.imshow(images[i][:,:]*boundary_mask*peak_mask_inv).set_clim(100,2000)
plt.subplot(224)
plt.title("throughbeam intensity")
plt.imshow(images[i][:,:]*throughbeam_mask).set_clim(100,2000)
plt.show()