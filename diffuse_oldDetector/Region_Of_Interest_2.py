"""
RegionOfInterest
Return a set of regions for selected areas using ginput.
Original by Xiazhe Shen
Adapted by Tyler Chase 03/29/2016
     -superfluous code removed
     -regions defined by one click instead of four clicks
"""

import matplotlib.pyplot as plt
from pylab import ginput

def GetRegionOfInterest(image, nRegion, contrastFactor=1.0, halfLength=40, 
                        flagSaveFig=False, saveName='./peak_number.jpg', message = ''):
    vmax = max(image.max(axis=0))*contrastFactor
    plt.figure(figsize=(15,15))
    plt.imshow(image, vmax=vmax)
    plt.title(message)
    plt.axis('off')
    plt.bar(0, 1, 1, 0, color="w", alpha=0.3)
    #above line added to temporarily fix graph reloading problem
    peakRegion = []
    backgroundRegion = []
    
    # Read in ROIs for diffraction peaks using ginput
    for i in range(nRegion):
        pts = ginput(1)
        colMin = max(int(pts[0][0])-halfLength, 0)
        colMax = min(int(pts[0][0])+halfLength, image.shape[0]-1)
        rowMin = max(int(pts[0][1])-halfLength, 0)
        rowMax = min(int(pts[0][1])+halfLength, image.shape[1]-1)
        height = rowMax - rowMin
        width  = colMax - colMin
        # Plot chosen ROI on image
        plt.bar(colMin, height, width, rowMin, color="w", alpha=0.3)
        plt.annotate(str(i+1), xy=(colMax, rowMax), color="w")
        peakRegion.append([rowMin, rowMax, colMin, colMax])
        
        # Read in ROIs for background near a given peak using ginput
    for i in range(nRegion):
        pts = ginput(1)
        colMin = max(int(pts[0][0])-halfLength, 0)
        colMax = min(int(pts[0][0])+halfLength, image.shape[0]-1)
        rowMin = max(int(pts[0][1])-halfLength, 0)
        rowMax = min(int(pts[0][1])+halfLength, image.shape[1]-1)
        height = rowMax - rowMin
        width  = colMax - colMin
        plt.bar(colMin, height, width, rowMin, color="r", alpha=0.3)
        plt.annotate(str(i+1), xy=(colMax, rowMax), color="r")
        backgroundRegion.append([rowMin, rowMax, colMin, colMax])

    plt.close("all")
    
    # Return ROIs for peaks and background near peaks
    return [peakRegion, backgroundRegion]