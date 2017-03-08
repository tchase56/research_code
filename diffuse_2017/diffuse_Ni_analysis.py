# Diffuse Online Analysis 

# Libraries
import functions_20170302_Ni
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Plots
plot_diffuse = 0
clim = 5

# inputs
save_address = r'C:\Users\tchase56\Documents\UED\UED_2017\Ni\20170227\diffuse\\\\'

Ni_dispersion = np.load(r'C:\Users\tchase56\Documents\UED\UED_2017\Ni\20170227\Analysis\20170301\Ni_dispersion\Ni_dispersion.npy').item()
runsToAlign = []
#loadAlignTimescan = 0
runsToAnalyze = [1,2,3,4,5,6,7,8,9]
loadDebye = 0
loadBackground = 1
loadLinescan = 1
ROI = 'inner' # make sure loadLinescan has been run with the associated ROI 
                   #'inner', 'transverse', or 'longitudinal'


# Debye Constants
gaussian = 1
points = 5
numPeaks = 6
# Background Flags
subtract_background = 1
normalize = 1
expansion_interpolation = 1
symmetrize = 1
radius_mask = 600
# Delay stage information
#time_zero = 64.377
time_zero = 64.392
'''
delay_stage = np.array([64.305,  64.343,  64.417,  64.455,
        64.492,  64.530,  64.567,  64.642 ,
        64.717,  64.792,  64.942,  65.092,  66.778]) 
'''
delay_stage = np.array([64.3395, 64.3470, 64.3545, 64.3620, 64.3695, 64.3845, \
                        64.3920, 64.3995, 64.4070, 64.4145, 64.4519, 64.4894, \
                        64.5269, 64.5644, 64.6393, 64.7143, 64.7892, 64.9391, \
                        65.0890, 66.7753])
# Initial values for least square exponential fits
scale_diffuse = 0        
scaling_factors_tran = [1,1,2,1,1]      # scale fit intensities for non uniform ROI lengths     
scaling_factors_long = [1,1,1,1]
scaling_factors_inner = [1,1,1,1]     
tau_seed_debye = 0.5
A_seed_debye = 0.5
tau_seed_diffuse = 0.5
if scale_diffuse:
    A_seed_diffuse = 1
else:
    A_seed_diffuse = 1000
debye_region_bottom = 7
diffuse_region_bottom = 6

# Parameters for exponential fits for transverse ROI
timeZeroShift = 0
lower_tran = 25
upper_tran = 216     # (upper-lower) must be factor of 2
length_tran = 15   # must be factor of 2
debye_x_lowLim_tran, debye_x_upLim_tran, debye_y_lowLim_tran, debye_y_upLim_tran = [-1, 17, 0.91, 1.01]
fitInten_y_lowLim_tran, fitInten_y_upLim_tran = [-10, 2500]
diffuse_y_lowLim_tran, diffuse_y_upLim_tran = [-10, 150]
tau_y_lowLim_tran, tau_y_upLim_tran = [-0.75, 2.0]
highlights_debye_tran = [0, 1, 2, 3, 4] 
highlights_diffuse_tran = [0,2,4]
highlights_linescan_tran = [0,8,13,17]
start_tran = [30, 35, 60, 36, 30]     # 'center', 'right', or 'left'
half_tran = False
# Parameters for exponential fits for longitudinal ROI
lower_long = 19
upper_long = 120     # (upper-lower) must be factor of 2
length_long = 15     # must be factor of 2
debye_x_lowLim_long, debye_x_upLim_long, debye_y_lowLim_long, debye_y_upLim_long = [-1, 17, 0.90, 1.05]
fitInten_y_lowLim_long, fitInten_y_upLim_long = [-10, 1700]
diffuse_y_lowLim_long, diffuse_y_upLim_long = [-25, 225]
tau_y_lowLim_long, tau_y_upLim_long = [-0.75, 2.5]
highlights_debye_long = [0, 1, 2, 3, 4] 
highlights_diffuse_long = [1,3]
highlights_linescan_long = [0,8,13,17]
start_long = [5, 30, 36, 30]
half_long = True
# Parameters for exponential fits for inner ROI
lower_inner = 19
upper_inner = 120
length_inner = 15
debye_x_lowLim_inner, debye_x_upLim_inner, debye_y_lowLim_inner, debye_y_upLim_inner = [-1, 17, 0.90, 1.05]
fitInten_y_lowLim_inner, fitInten_y_upLim_inner = [0, 2800]
diffuse_y_lowLim_inner, diffuse_y_upLim_inner = [-25, 200]
tau_y_lowLim_inner, tau_y_upLim_inner = [-0.25, 1.8]
highlights_debye_inner = [0, 1, 2, 3, 4] 
highlights_diffuse_inner = [1,3]
highlights_linescan_inner = [0,8,13,17]
start_inner = [5, 30, 36, 30]
half_inner = True
    
# Convert Delay stage values to picoseconds
#delay_stage = np.arange(delayStage_start, delayStage_end + (delayStage_step)/2, delayStage_step)
times = functions_20170302_Ni.ps_convert(delay_stage, time_zero) 

if len(runsToAlign)>0:   
    ### Align unaligned runs
    load_addresses_align = []
    for i in range(len(runsToAlign)):
        temp = save_address + r'scan' + str(runsToAlign[i]) + '\images-ANDOR1\\'
        load_addresses_align.append(temp)      
    # initialize run arrays
    runs = []
    runs_aligned_centered = [] 
    # attain first averaged and averaged/aligned run to form ROIs
    runs.append(functions_20170302_Ni.load_timescan(load_addresses_align[0], delay_stage))
    temp, peak_region, background_region = functions_20170302_Ni.load_timescan(load_addresses_align[0], delay_stage, align = 1, center = 1)
    runs_aligned_centered.append(temp)
    print('scan ' + str(runsToAlign[0]))    
    # Load runs, average, center and average
    for i in range(1,len(runsToAlign)):
        print('scan ' + str(runsToAlign[i]))
        runs.append(functions_20170302_Ni.load_timescan(load_addresses_align[i], delay_stage))
        temp, _, _ = functions_20170302_Ni.load_timescan(load_addresses_align[i], delay_stage, align = 1, center = 1, peak_region = peak_region, background_region = background_region)
        runs_aligned_centered.append(temp)      
    for i in range(len(runsToAlign)):
        np.save(load_addresses_align[i] + 'noAlign_quick.npy', runs[i])
        np.save(load_addresses_align[i] + 'align_center_quick.npy', runs_aligned_centered[i])

    

    


### Analyze runs
# Create array of load addresses
load_addresses_analyze = []
runs_unAligned = []
runs_aligned = []
for i in range(len(runsToAnalyze)):
    temp = save_address + r'scan' + str(runsToAnalyze[i]) + '\images-ANDOR1\\'
    load_addresses_analyze.append(temp)
for i in range(len(runsToAnalyze)):
    runs_unAligned.append(np.load(load_addresses_analyze[i] + 'noAlign_quick.npy'))
    runs_aligned.append(np.load(load_addresses_analyze[i] + 'align_center_quick.npy'))
    
# Average over the runs to analyze 
runs_ave = np.average(runs_unAligned, axis = 0)
runs_aligned_centered_ave = np.average(runs_aligned, axis = 0)

# Do Debye-Waller Analysis and Peak-Expansion Analysis
if not loadDebye:
    peak_expansion = np.load(save_address + "peak_expansion.npy")
    static_q = np.load(save_address + "static_q.npy")
    if not gaussian:
        debye = np.load(save_address + "debye_intensity.npy")
    else:
        debye = np.load(save_address + "debye_intensity_gaussian.npy")
else:
    #[debye, peak_expansion] = functions_20170302_Ni.debyeWaller(save_address, delay_stage_timescan, runs_aligned_centered_ave_timescan, points, numPeaks = numPeaks, gaussian = gaussian, plot = True)
    [debye, peak_expansion, static_q] = functions_20170302_Ni.debyeWaller(save_address, delay_stage, runs_aligned_centered_ave, points, numPeaks = numPeaks, gaussian = gaussian, plot = True)

# Subtract detector background, normalize detector intensity, interpolate lattice expansion, symmetrize diffraction
if not loadBackground:
    # Load the data
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 1 and symmetrize == 1:
        rotated = np.load(save_address + 'expanded_norm_subBack_symmetrized.npy')
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 1 and symmetrize == 1:
        rotated = np.load(save_address + 'norm_subBack_symmetrized.npy')
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 1 and symmetrize == 1:
        rotated = np.load(save_address + 'expanded_subBack_symmetrized.npy')
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 0 and symmetrize == 1:
        rotated = np.load(save_address + 'expanded_norm_symmetrized.npy')
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 1 and symmetrize == 0:
        rotated = np.load(save_address + 'expanded_norm_subBack.npy')
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 1 and symmetrize == 1:
        rotated = np.load(save_address + 'subBack_symmetrized.npy')        
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 0 and symmetrize == 1:
        rotated = np.load(save_address + 'norm_symmetrized.npy')
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 1 and symmetrize == 0:
        rotated = np.load(save_address + 'norm_subBack.npy')
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 0 and symmetrize == 1:
        rotated = np.load(save_address + 'expanded_symmetrized.npy')
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 1 and symmetrize == 0:
        rotated = np.load(save_address + 'expanded_subBack.npy')
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 0 and symmetrize == 0:
        rotated = np.load(save_address + 'expanded_norm.npy')
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 0 and symmetrize == 1:
        rotated = np.load(save_address + 'symmetrized.npy')
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 1 and symmetrize == 0:
        rotated = np.load(save_address + 'subBack.npy')
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 0 and symmetrize == 0:
        rotated = np.load(save_address + 'norm.npy')
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 0 and symmetrize == 0:
        rotated = np.load(save_address + 'expanded.npy')
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 0 and symmetrize == 0:
        rotated = np.load(save_address + 'rotated.npy')        
else:
    rotated = functions_20170302_Ni.backgroundNormExpansionSym(save_address, runs_ave, runs_aligned_centered_ave, times, radius_mask, subtract_background, normalize, expansion_interpolation, symmetrize)
    
if plot_diffuse:
    for i in range(1,len(rotated)):
        plt.figure()
        plt.imshow(rotated[i]-rotated[0]).set_clim(0,clim)
        plt.colorbar()
        plt.title('{:.2f} ps'.format(times[i]))
        
# Select region of interest for high symmetry linescan analysis
if not loadLinescan:
    recyp_lat = np.load(save_address + 'recyprocal_lattice_points_gammaToGamma.npy').item()
    ROI_tran = np.load(save_address + 'ROI_tran_linescan.npy')
    ROI_long = np.load(save_address + 'ROI_long_linescan.npy')
    ROI_inner = np.load(save_address + 'ROI_inner_linescan.npy')
    ROI_x_min_tran = ROI_tran[0]
    ROI_y_min_tran = ROI_tran[1]
    width_tran = ROI_tran[2]
    height_tran = ROI_tran[3]
    ROI_x_min_long = ROI_long[0]
    ROI_y_min_long = ROI_long[1]
    width_long = ROI_long[2]
    height_long = ROI_long[3]
    ROI_x_min_inner = ROI_inner[0]
    ROI_y_min_inner = ROI_inner[1]
    width_inner = ROI_inner[2]
    height_inner = ROI_inner[3]
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 1 and symmetrize == 1:
        diffuse_tran = np.load(save_address + 'diffuseTran_expanded_norm_subBack_symmetrized.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_expanded_norm_subBack_symmetrized.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_expanded_norm_subBack_symmetrized.npy')
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 1 and symmetrize == 1:
        diffuse_tran = np.load(save_address + 'diffuseTran_norm_subBack_symmetrized.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_norm_subBack_symmetrized.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_norm_subBack_symmetrized.npy')
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 1 and symmetrize == 1:
        diffuse_tran = np.load(save_address + 'diffuseTran_expanded_subBack_symmetrized.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_expanded_subBack_symmetrized.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_expanded_subBack_symmetrized.npy')
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 0 and symmetrize == 1:
        diffuse_tran = np.load(save_address + 'diffuseTran_expanded_norm_symmetrized.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_expanded_norm_symmetrized.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_expanded_norm_symmetrized.npy')
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 1 and symmetrize == 0:
        diffuse_tran = np.load(save_address + 'diffuseTran_expanded_norm_subBack.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_expanded_norm_subBack.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_expanded_norm_subBack.npy')
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 1 and symmetrize == 1:
        diffuse_tran = np.load(save_address + 'diffuseTran_subBack_symmetrized.npy')        
        diffuse_long = np.load(save_address + 'diffuseLong_subBack_symmetrized.npy')        
        diffuse_inner = np.load(save_address + 'diffuseInner_subBack_symmetrized.npy')        
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 0 and symmetrize == 1:
        diffuse_tran = np.load(save_address + 'diffuseTran_norm_symmetrized.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_norm_symmetrized.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_norm_symmetrized.npy')
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 1 and symmetrize == 0:
        diffuse_tran = np.load(save_address + 'diffuseTran_norm_subBack.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_norm_subBack.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_norm_subBack.npy')
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 0 and symmetrize == 1:
        diffuse_tran = np.load(save_address + 'diffuseTran_expanded_symmetrized.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_expanded_symmetrized.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_expanded_symmetrized.npy')
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 1 and symmetrize == 0:
        diffuse_tran = np.load(save_address + 'diffuseTran_expanded_subBack.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_expanded_subBack.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_expanded_subBack.npy')
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 0 and symmetrize == 0:
        diffuse_tran = np.load(save_address + 'diffuseTran_expanded_norm.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_expanded_norm.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_expanded_norm.npy')
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 0 and symmetrize == 1:
        diffuse_tran = np.load(save_address + 'diffuseTran_symmetrized.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_symmetrized.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_symmetrized.npy')
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 1 and symmetrize == 0:
        diffuse_tran = np.load(save_address + 'diffuseTran_subBack.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_subBack.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_subBack.npy')
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 0 and symmetrize == 0:
        diffuse_tran = np.load(save_address + 'diffuseTran_norm.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_norm.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_norm.npy')
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 0 and symmetrize == 0:
        diffuse_tran = np.load(save_address + 'diffuseTran_expanded.npy')
        diffuse_long = np.load(save_address + 'diffuseLong_expanded.npy')
        diffuse_inner = np.load(save_address + 'diffuseInner_expanded.npy')
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 0 and symmetrize == 0:
        diffuse_tran = np.load(save_address + 'diffuseTran_rotated.npy')   
        diffuse_long = np.load(save_address + 'diffuseLong_rotated.npy')   
        diffuse_inner = np.load(save_address + 'diffuseInner_rotated.npy')   
   
    # Mask out the center and edges of diffraction pattern
    mask = functions_20170302_Ni.make_mask(rotated[0], 220, int(1024/2-66))
   
    # Plot the region of interest we are integrating to form the linescan
    # Region of interest for transverse ROI
    plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(mask*(rotated[-1]-rotated[0]), origin = 'upper').set_clim(0, 5)
    plt.title('Transverse ROI', fontsize = 20)
    ax1.add_patch(
        patches.Rectangle((ROI_x_min_tran, ROI_y_min_tran), width_tran, height_tran, fill = False, edgecolor = 'white', lw = 5)
    )
    ax1.set_xlim(220, 1024-220)
    ax1.set_ylim(1024-220, 220)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)  
    plt.show()   
    
    # Region of interest for longitudinal ROI
    plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(mask*(rotated[-1]-rotated[0]), origin = 'upper').set_clim(0, 5)
    plt.title('Longitudinal ROI', fontsize = 20)
    ax1.add_patch(
        patches.Rectangle((ROI_x_min_long, ROI_y_min_long), width_long, height_long, fill = False, edgecolor = 'white', lw = 5)
    )
    ax1.set_xlim(220, 1024-220)
    ax1.set_ylim(1024-220, 220)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)  
    plt.show() 
    
    # Region of interest for inner ROI
    plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(mask*(rotated[-1]-rotated[0]), origin = 'upper').set_clim(0, 5)
    plt.title('Inner ROI', fontsize = 20)
    ax1.add_patch(
        patches.Rectangle((ROI_x_min_inner, ROI_y_min_inner), width_inner, height_inner, fill = False, edgecolor = 'white', lw = 5)
    )
    ax1.set_xlim(220, 1024-220)
    ax1.set_ylim(1024-220, 220)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)  
    plt.show() 
else:
    # Calculate recyprocal lattice points and save them
    recyp_lat = functions_20170302_Ni.recyprocal_lattice_points(rotated[0], halfFlag = 0)
    np.save(save_address + 'recyprocal_lattice_points_gammaToGamma.npy', recyp_lat)
    # Subtract reference at negative time delay
    runs_sub = functions_20170302_Ni.difference_image(rotated)
    # Form diffuse linescans and the associated regions of interest, then save
    diffuse_tran, ROI_tran = functions_20170302_Ni.formDiffuseLinescan(times, rotated, runs_sub, plot = True, 
                                                                    message = "choose two peaks along virtical for transverse ROI")
    diffuse_long, ROI_long = functions_20170302_Ni.formDiffuseLinescan(times, rotated, runs_sub, plot = True, 
                                                                    message = "choose two peaks along virtical for longitudinal ROI")
    diffuse_inner, ROI_inner = functions_20170302_Ni.formDiffuseLinescan(times, rotated, runs_sub, plot = True, 
                                                                    message = "choose two peaks along virtical for inner ROI")
    np.save(save_address + 'ROI_tran_linescan.npy', ROI_tran)
    np.save(save_address + 'ROI_long_linescan.npy', ROI_long)
    np.save(save_address + 'ROI_inner_linescan.npy', ROI_inner)
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 1 and symmetrize == 1:
        np.save(save_address + 'diffuseTran_expanded_norm_subBack_symmetrized.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_expanded_norm_subBack_symmetrized.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_expanded_norm_subBack_symmetrized.npy', diffuse_inner)
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 1 and symmetrize == 1:
        np.save(save_address + 'diffuseTran_norm_subBack_symmetrized.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_norm_subBack_symmetrized.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_norm_subBack_symmetrized.npy', diffuse_inner)
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 1 and symmetrize == 1:
        np.save(save_address + 'diffuseTran_expanded_subBack_symmetrized.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_expanded_subBack_symmetrized.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_expanded_subBack_symmetrized.npy', diffuse_inner)
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 0 and symmetrize == 1:
        np.save(save_address + 'diffuseTran_expanded_norm_symmetrized.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_expanded_norm_symmetrized.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_expanded_norm_symmetrized.npy', diffuse_inner)
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 1 and symmetrize == 0:
        np.save(save_address + 'diffuseTran_expanded_norm_subBack.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_expanded_norm_subBack.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_expanded_norm_subBack.npy', diffuse_inner)
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 1 and symmetrize == 1:
        np.save(save_address + 'diffuseTran_subBack_symmetrized.npy', diffuse_tran)      
        np.save(save_address + 'diffuseLong_subBack_symmetrized.npy', diffuse_long)        
        np.save(save_address + 'diffuseInner_subBack_symmetrized.npy', diffuse_inner)        
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 0 and symmetrize == 1:
        np.save(save_address + 'diffuseTran_norm_symmetrized.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_norm_symmetrized.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_norm_symmetrized.npy', diffuse_inner)
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 1 and symmetrize == 0:
        np.save(save_address + 'diffuseTran_norm_subBack.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_norm_subBack.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_norm_subBack.npy', diffuse_inner)
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 0 and symmetrize == 1:
        np.save(save_address + 'diffuseTran_expanded_symmetrized.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_expanded_symmetrized.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_expanded_symmetrized.npy', diffuse_inner)
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 1 and symmetrize == 0:
        np.save(save_address + 'diffuseTran_expanded_subBack.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_expanded_subBack.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_expanded_subBack.npy', diffuse_inner)
    if expansion_interpolation == 1 and normalize == 1 and subtract_background == 0 and symmetrize == 0:
        np.save(save_address + 'diffuseTran_expanded_norm.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_expanded_norm.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_expanded_norm.npy', diffuse_inner)
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 0 and symmetrize == 1:
        np.save(save_address + 'diffuseTran_symmetrized.npy', diffuse_tran)
        np.save(save_address + 'diffuseTran_symmetrized.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_symmetrized.npy', diffuse_inner)
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 1 and symmetrize == 0:
        np.save(save_address + 'diffuseTran_subBack.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_subBack.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_subBack.npy', diffuse_inner)
    if expansion_interpolation == 0 and normalize == 1 and subtract_background == 0 and symmetrize == 0:
        np.save(save_address + 'diffuseTran_norm.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_norm.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_norm.npy', diffuse_inner)
    if expansion_interpolation == 1 and normalize == 0 and subtract_background == 0 and symmetrize == 0:
        np.save(save_address + 'diffuseTran_expanded.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_expanded.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_expanded.npy', diffuse_inner)
    if expansion_interpolation == 0 and normalize == 0 and subtract_background == 0 and symmetrize == 0:
        np.save(save_address + 'diffuseTran_rotated.npy', diffuse_tran)
        np.save(save_address + 'diffuseLong_rotated.npy', diffuse_long)
        np.save(save_address + 'diffuseInner_rotated.npy', diffuse_inner)


static_q = static_q/static_q    
# Fit timescales of thermalization for linescans
times_timescan = times
tau_debye_tran, tau_error_debye_tran, tau_diffuse_tran, tau_error_diffuse_tran,\
 gamma_to_gamma_percent_tran, bins_tran, sums_tran, lineScan_short_tran, scale_const_tran = functions_20170302_Ni.thermalization(times, diffuse_tran, times_timescan, recyp_lat, 
                                                                            timeZeroShift, debye, lower_tran, upper_tran, length_tran, 
                                                                            debye_x_lowLim_tran, debye_x_upLim_tran, debye_y_lowLim_tran, debye_y_upLim_tran, 
                                                                            fitInten_y_lowLim_tran, fitInten_y_upLim_tran, diffuse_y_lowLim_tran, 
                                                                            diffuse_y_upLim_tran, tau_y_lowLim_tran, tau_y_upLim_tran, tau_seed_debye, 
                                                                            A_seed_debye, tau_seed_diffuse, A_seed_diffuse, debye_region_bottom,
                                                                            diffuse_region_bottom, highlights_debye_tran, highlights_diffuse_tran, static_q, scale_diffuse,
                                                                            ROI_title = 'Transverse ROI', plot = False, start = start_tran, half = half_tran, scaling_factors = scaling_factors_tran)                                                                                                                    
tau_debye_long, tau_error_debye_long, tau_diffuse_long, tau_error_diffuse_long,\
 gamma_to_gamma_percent_long, bins_long, sums_long, lineScan_short_long, scale_const_long = functions_20170302_Ni.thermalization(times, diffuse_long, times_timescan, recyp_lat, 
                                                                            timeZeroShift, debye, lower_long, upper_long, length_long, 
                                                                            debye_x_lowLim_long, debye_x_upLim_long, debye_y_lowLim_long, debye_y_upLim_long, 
                                                                            fitInten_y_lowLim_long, fitInten_y_upLim_long, diffuse_y_lowLim_long, 
                                                                            diffuse_y_upLim_long, tau_y_lowLim_long, tau_y_upLim_long, tau_seed_debye, 
                                                                            A_seed_debye, tau_seed_diffuse, A_seed_diffuse, debye_region_bottom,
                                                                            diffuse_region_bottom, highlights_debye_long, highlights_diffuse_long, static_q, scale_diffuse,
                                                                            ROI_title = 'Longitudinal ROI', plot = False, start = start_long, half = half_long, scaling_factors = scaling_factors_long)                                                                                                                  
tau_debye_inner, tau_error_debye_inner, tau_diffuse_inner, tau_error_diffuse_inner, \
gamma_to_gamma_percent_inner, bins_inner, sums_inner, lineScan_short_inner, scale_const_inner = functions_20170302_Ni.thermalization(times, diffuse_inner, times_timescan, recyp_lat, 
                                                                             timeZeroShift, debye, lower_inner, upper_inner, length_inner, 
                                                                             debye_x_lowLim_inner, debye_x_upLim_inner, debye_y_lowLim_inner, debye_y_upLim_inner, 
                                                                             fitInten_y_lowLim_inner, fitInten_y_upLim_inner, diffuse_y_lowLim_inner, 
                                                                             diffuse_y_upLim_inner, tau_y_lowLim_inner, tau_y_upLim_inner, tau_seed_debye, 
                                                                             A_seed_debye, tau_seed_diffuse, A_seed_diffuse, debye_region_bottom,
                                                                             diffuse_region_bottom, highlights_debye_inner, highlights_diffuse_inner, static_q, scale_diffuse, 
                                                                             ROI_title = 'Inner ROI', plot = False, start = start_inner, half = half_inner, scaling_factors = scaling_factors_inner)
                                                                             
''' Print Time Constants'''
# Print the time constants for the debye-waller
peak_names = ['220', '400', '420', '440', '600']
# Print time constant fits for Debye-Waller
print('\nDebye-Waller')
print("\nPeak\t" + "time_constant")
for i in range(len(peak_names)):
    print( peak_names[i] + '\t' + str(round(np.array(tau_debye_tran)[i,0],2)) + '+-' + str(round(tau_error_debye_tran[i],2)) )
print('\n')
# Print the time constants for the diffuse scattering transverse ROI
print('ROI Transverse')
print("\ngamma_to_gamma_%\t" + 'bin_width\t' + "time_constant")
for i in range(np.shape(tau_error_diffuse_tran)[0]):
    print(str(round(gamma_to_gamma_percent_tran[i],10)) + '\t\t' + str(bins_tran[i+1]-bins_tran[i]) + '\t\t' + str(round(np.array(tau_diffuse_tran)[i,0],2)) + ' +- ' + str(round(np.array(tau_error_diffuse_tran)[i,0,0], 2))) 
print('\n')
# Print the time constants for the diffuse scattering transverse ROI
print('ROI Longitudinal')
print("\ngamma_to_gamma_%\t" + 'bin_width\t' + "time_constant")
for i in range(np.shape(tau_error_diffuse_long)[0]):
    print(str(round(gamma_to_gamma_percent_long[i],10)) + '\t\t' + str(bins_long[i+1]-bins_long[i]) + '\t\t' + str(round(np.array(tau_diffuse_long)[i,0],2)) + ' +- ' + str(round(np.array(tau_error_diffuse_long)[i,0,0], 2))) 
print('\n')
# Print the time constants for the diffuse scattering transverse ROI
print('ROI Inner')
print("\ngamma_to_gamma_%\t" + 'bin_width\t' + "time_constant")
for i in range(np.shape(tau_error_diffuse_inner)[0]):
    print(str(round(gamma_to_gamma_percent_inner[i],10)) + '\t\t' + str(bins_inner[i+1]-bins_inner[i]) + '\t\t' + str(round(np.array(tau_diffuse_inner)[i,0],2)) + ' +- ' + str(round(np.array(tau_error_diffuse_inner)[i,0,0], 2))) 
print('\n')

# Plot time constant comparison
# Compare gamma and x for transverse ROI
functions_20170302_Ni.thermPlot(times, recyp_lat,  diffuse_region_bottom,
          lineScan_short_tran, tau_diffuse_tran, tau_error_diffuse_tran, gamma_to_gamma_percent_tran, sums_tran, highlights_diffuse_tran, bins_tran,
          length_tran, diffuse_y_lowLim_tran, diffuse_y_upLim_tran, fitInten_y_lowLim_tran, fitInten_y_upLim_tran, tau_y_lowLim_tran, tau_y_upLim_tran, highlights_linescan_tran, scale_const_tran, half = half_tran,
          scale_diffuse = scale_diffuse, ROI_title = 'Transverse_ROI', label = ['16% $\Gamma$ to $\Gamma$', '49% $\Gamma$ to $\Gamma$', '83% $\Gamma$ to $\Gamma$'])          
# Compare gamma and x for longitudinal ROI
functions_20170302_Ni.thermPlot(times, recyp_lat,  diffuse_region_bottom,
          lineScan_short_long, tau_diffuse_long, tau_error_diffuse_long, gamma_to_gamma_percent_long, sums_long, highlights_diffuse_long, bins_long,
          length_long, diffuse_y_lowLim_long, diffuse_y_upLim_long, fitInten_y_lowLim_long, fitInten_y_upLim_long, tau_y_lowLim_long, tau_y_upLim_long, highlights_linescan_long, scale_const_long, half = half_long,
          scale_diffuse = scale_diffuse, ROI_title = 'Longitudinal_ROI', label = ['16% $\Gamma$ to $\Gamma$', '49% $\Gamma$ to $\Gamma$'])
# Compare gamma and x for Inner ROI
functions_20170302_Ni.thermPlot(times, recyp_lat,  diffuse_region_bottom,
          lineScan_short_inner, tau_diffuse_inner, tau_error_diffuse_inner, gamma_to_gamma_percent_inner, sums_inner, highlights_diffuse_inner, bins_inner,
          length_inner, diffuse_y_lowLim_inner, diffuse_y_upLim_inner, fitInten_y_lowLim_inner, fitInten_y_upLim_inner, tau_y_lowLim_inner, tau_y_upLim_inner, highlights_linescan_inner, scale_const_inner, half = half_inner,
          scale_diffuse = scale_diffuse, ROI_title = 'Inner_ROI', label = ['16% $\Gamma$ to $\Gamma$', '49% $\Gamma$ to $\Gamma$'])
          
# Compare x from transverse ROI and x from longitudinal ROI         
functions_20170302_Ni.thermPlot(times, recyp_lat,  diffuse_region_bottom,
          lineScan_short_tran, tau_diffuse_tran, tau_error_diffuse_tran, gamma_to_gamma_percent_tran, sums_tran, [2], bins_tran,
          length_tran, diffuse_y_lowLim_tran, diffuse_y_upLim_tran, fitInten_y_lowLim_tran, fitInten_y_upLim_tran, tau_y_lowLim_tran, tau_y_upLim_tran, highlights_linescan_tran, scale_const_tran, half = half_tran,
          flag_2 = True, lineScan_short_2 = lineScan_short_inner, tau_diffuse_2 = tau_diffuse_inner, error_diffuse_2 = tau_error_diffuse_inner, gamma_to_gamma_percent_2 = gamma_to_gamma_percent_inner, sums_2 = sums_inner, highlights_diffuse_2 = [3], bins_2 = bins_inner,
          length_2 = length_inner, diffuse_y_lowLim_2 = diffuse_y_lowLim_inner, diffuse_y_upLim_2 = diffuse_y_upLim_inner, fitInten_y_lowLim_2 = fitInten_y_lowLim_inner, fitInten_y_upLim_2 = fitInten_y_upLim_inner, tau_y_lowLim_2 = tau_y_lowLim_inner, tau_y_upLim_2 = tau_y_upLim_inner, highlights_linescan_2 = highlights_linescan_inner, scale_const_2 = scale_const_inner, half_2 = half_inner,
          scale_diffuse = scale_diffuse, ROI_title = 'Transverse vs Inner X ROI', label = ['49% $\Gamma$ to $\Gamma$ Transverse', '49% $\Gamma$ to $\Gamma$ Inner'])
          
# Compare x from transverse ROI and x from longitudinal ROI         
functions_20170302_Ni.thermPlot(times, recyp_lat,  diffuse_region_bottom,
          lineScan_short_tran, tau_diffuse_tran, tau_error_diffuse_tran, gamma_to_gamma_percent_tran, sums_tran, [0,4], bins_tran,
          length_tran, diffuse_y_lowLim_tran, diffuse_y_upLim_tran, fitInten_y_lowLim_tran, fitInten_y_upLim_tran, tau_y_lowLim_tran, tau_y_upLim_tran, highlights_linescan_tran, scale_const_tran, half = half_tran,
          flag_2 = True, lineScan_short_2 = lineScan_short_inner, tau_diffuse_2 = tau_diffuse_inner, error_diffuse_2 = tau_error_diffuse_inner, gamma_to_gamma_percent_2 = gamma_to_gamma_percent_inner, sums_2 = sums_inner, highlights_diffuse_2 = [1], bins_2 = bins_inner,
          length_2 = length_inner, diffuse_y_lowLim_2 = diffuse_y_lowLim_inner, diffuse_y_upLim_2 = diffuse_y_upLim_inner, fitInten_y_lowLim_2 = fitInten_y_lowLim_inner, fitInten_y_upLim_2 = fitInten_y_upLim_inner, tau_y_lowLim_2 = tau_y_lowLim_inner, tau_y_upLim_2 = tau_y_upLim_inner, highlights_linescan_2 = highlights_linescan_inner, scale_const_2 = scale_const_inner, half_2 = half_inner,
          scale_diffuse = scale_diffuse, ROI_title = 'Transverse vs Inner $\Gamma$ ROI', label = ['16% $\Gamma$ to $\Gamma$ Transverse', '83% $\Gamma$ to $\Gamma$ Transverse', '16% $\Gamma$ to $\Gamma$ Inner'])
          
# Plot debye-waller graphs
functions_20170302_Ni.debPlot(times, debye, highlights_debye_tran, 
                              tau_debye_tran, tau_error_debye_tran, debye_region_bottom,
                              debye_x_lowLim_tran, debye_x_upLim_tran, debye_y_lowLim_tran, debye_y_upLim_tran, static_q)
                              
functions_20170302_Ni.debPlot(times, debye, [1,2,4], 
                              tau_debye_tran, tau_error_debye_tran, debye_region_bottom,
                              debye_x_lowLim_tran, debye_x_upLim_tran, debye_y_lowLim_tran, debye_y_upLim_tran, static_q)
  
                            
# Plot diffuse with dispersion
# Compare gamma and x for longitudinal ROI
disp_y_lowLim = 0
disp_y_upLim = 10                          

[disp_y_lowLim, disp_y_upLim] = [0,10]
functions_20170302_Ni.dispPlot(times, recyp_lat,  diffuse_region_bottom,
          lineScan_short_tran, Ni_dispersion,
          -10, 100, fitInten_y_lowLim_tran, fitInten_y_upLim_tran, disp_y_lowLim, disp_y_upLim, highlights_linescan_tran, half = half_tran,
          ROI_title = 'Transverse_ROI', interpolate_disp = True)          

functions_20170302_Ni.dispPlot(times, recyp_lat,  diffuse_region_bottom,
          lineScan_short_inner, Ni_dispersion,
          diffuse_y_lowLim_inner, 150, fitInten_y_lowLim_inner, fitInten_y_upLim_inner, disp_y_lowLim, disp_y_upLim, highlights_linescan_inner, half = half_inner,
          ROI_title = 'Inner_ROI', interpolate_disp = True)
          
functions_20170302_Ni.dispPlot(times, recyp_lat,  diffuse_region_bottom,
          lineScan_short_long, Ni_dispersion,
          diffuse_y_lowLim_long, 201, fitInten_y_lowLim_long, fitInten_y_upLim_long, disp_y_lowLim, disp_y_upLim, highlights_linescan_long, half = half_long,
          ROI_title = 'Longitudinal_ROI', interpolate_disp = True) 


        

    

