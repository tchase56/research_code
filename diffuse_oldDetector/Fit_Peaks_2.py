'''
Fit_Peaks
Fit peaks using a 1d gaussian
Original by Xiazhe Shen
Adapted by Tyler Chase 03/29/2016
     -some superfluous code removed
     -adapted for latest version of python
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
#from pylab import *
import pylab as py
from scipy.optimize import curve_fit
from statistics import median 


# Gaussian function for determining peak positions
def OneDimensionalGaussian(x, amp, mu, sigma, offset):
    return amp*np.exp(-0.5*(x-mu)**2/(sigma**2)) + offset



def FitOnePeak(image, region, background, flagPlot=False):
    print("1D Gaussian fitting")
    
    # Determine the average background per pixel
    rowMin = background[0]
    rowMax = background[1]
    colMin = background[2]
    colMax = background[3]    
    
    background_subImage = image[rowMin:rowMax+1, colMin:colMax+1]
    backgroundPerPixel = median(background_subImage.ravel())
    
    
    # Read in tentative row and column minimums and maximums
    rowMin = region[0]
    rowMax = region[1]
    colMin = region[2]
    colMax = region[3]


    # Determine the gaussian center for both the rows and columns
    subImage = image[rowMin:rowMax+1, colMin:colMax+1]-backgroundPerPixel
    colPro = subImage.sum(axis=0)
    muCol = py.argmax(colPro)
    rowPro = subImage.sum(axis=1)
    muRow = py.argmax(rowPro)
    muRow = rowMin + muRow
    muCol = colMin + muCol


    # Determine actual row and column minimums and maximums
    halfLength = 50
    colMin = muCol - halfLength
    colMax = muCol + halfLength
    rowMin = muRow - halfLength
    rowMax = muRow + halfLength
    if rowMin < 0:
        rowMin = 0
    elif rowMax > image.shape[0]:
        rowMax = image.shape[0]
    if colMin < 0:
        colMin = 0
    elif colMax > image.shape[1]:
        colMax = image.shape[1]
    
            
    # Determine gaussian parameters for fitting a peak
    subImage = image[rowMin:rowMax+1, colMin:colMax+1]
    rowPro = subImage.sum(axis=1)
    colPro = subImage.sum(axis=0)
    colProFiltered = ndimage.filters.gaussian_filter(colPro, 2)
    rowProFiltered = ndimage.filters.gaussian_filter(rowPro, 2)
    muRow = py.argmax(rowProFiltered)
    muCol = py.argmax(colProFiltered)   
    if rowMin == 0:
        offsetRow = py.mean(subImage[-1,:])
    elif rowMax == image.shape[0]:
        offsetRow = py.mean(subImage[0,:])
    else:
        offsetRow = 0.5 * ( py.mean(subImage[0,:]) + py.mean(subImage[-1,:]) )  
    if colMin ==0:
        offsetCol = py.mean(subImage[:,-1])
    elif colMax == image.shape[1]:
        offsetCol = py.mean(subImage[:,0])
    else:
        offsetCol = 0.5 * ( py.mean(subImage[:,0]) + py.mean(subImage[:,-1]) )
    ampRow = rowProFiltered[muRow] - offsetRow
    sigRow = 10
    ampCol = colProFiltered[muCol] - offsetCol
    sigCol = 10


    # Fit a 1D gaussian to the row and column projections (if the regression fails use the gaussian parameters used above)
    fitRangeRow = np.arange(len(rowProFiltered))
    fitValueRow = np.array(rowProFiltered)
    # If fit a runtime error occurs (most likely regression times out) print out error and use guesses determined about
    try:
        paraRow, covRow = curve_fit(OneDimensionalGaussian, fitRangeRow, fitValueRow, p0=(ampRow, muRow, sigRow, offsetRow))
    except Exception as e:
        print("%s\n"%e)
        paraRow = [ampRow, muRow, sigRow, offsetRow]
        
    fitRangeCol = np.arange(len(colProFiltered))
    fitValueCol = np.array(colProFiltered)
    # If fit a runtime error occurs (most likely regression times out) print out error and use guesses determined about
    try:
        paraCol, covCol = curve_fit(OneDimensionalGaussian, fitRangeCol, fitValueCol, p0=(ampCol, muCol, sigCol, offsetCol))
    except Exception as e:
        print("%s\n"%e)
        paraCol= [ampCol, muCol, sigCol, offsetCol]


    # Determine fits for if plotting is flagged
    rowAxis = np.arange(len(rowPro))
    colAxis = np.arange(len(colPro))
    rowProFitted = OneDimensionalGaussian(fitRangeRow, *paraRow)
    colProFitted = OneDimensionalGaussian(fitRangeCol, *paraCol)
    colPro = colPro/colPro.max()*len(colPro)*0.3
    rowPro = rowPro/rowPro.max()*len(rowPro)*0.3
    rowProFitted = rowProFitted/rowProFitted.max()*len(rowProFitted)*0.3
    colProFitted = colProFitted/colProFitted.max()*len(colProFitted)*0.3


    # If plotting is flagged plot the fits
    if flagPlot:
        fig=plt.figure()
        plt.imshow(subImage, origin='lower')
        plt.colorbar()
        plt.plot(colAxis, colPro, 'ro')
        plt.plot(rowPro, rowAxis, 'go')
        plt.plot(colAxis, colProFitted, 'w-')
        plt.plot(rowProFitted, rowAxis, 'w-')
        plt.plot(paraCol[1], paraRow[1], 'wx', markersize=20)
        plt.xlim([0, subImage.shape[0]])
        plt.ylim([0, subImage.shape[1]])
        plt.xlabel('X [pixel]')
        plt.ylabel('Y [pixel]')
        fig.tight_layout()
        plt.show()

    # Correct the peak positions by adding the start of the region of interest
    paraCol[1] = paraCol[1] + colMin
    paraRow[1] = paraRow[1] + rowMin


    # Return the gaussian fit parameters
    return [paraCol, paraRow]
    
# Takes in image, peak ROIs, background ROIs, and returns fit parameters, peak intensity, and background intensity
def FitPeaks(image, regionSet, backgroundRegion, halfLength=40, flagPlot=False):
    peakSet = []
    peakSquareIntensitySum = [] # The sum of intensity around the peak locations
    backgroundSquareIntensitySum = []
    for i in range(len(regionSet)):
        peakSet.append(FitOnePeak(image, regionSet[i], backgroundRegion[i], flagPlot) )
        colMin = max(peakSet[i][0][1] - halfLength, 0)
        colMax = min(peakSet[i][0][1] + halfLength, image.shape[0]-1)
        rowMin = max(peakSet[i][1][1] - halfLength, 0)
        rowMax = min(peakSet[i][1][1] + halfLength, image.shape[1]-1)
        peakSquare = image[rowMin:rowMax, colMin:colMax]
        peakSquareIntensitySum.append(sum(peakSquare.sum(axis=0)))
        
        rowMin = backgroundRegion[i][0]
        rowMax = backgroundRegion[i][1]
        colMin = backgroundRegion[i][2]
        colMax = backgroundRegion[i][3]
        backgroundSquare = image[rowMin:rowMax, colMin:colMax]
        backgroundSquareIntensitySum.append(sum(backgroundSquare.sum(axis=0)))
#        background_average = sum(backgroundSquare.sum(axis=0))/(np.shape(backgroundSquare[0])*np.shape(backgroundSquare[1]))
        


        # If plotting is flagged plot the fits
        if flagPlot:
            fig = plt.figure(figsize=(18,6))
            plt.subplot(131)
            plt.title('Peak')
            plt.imshow(peakSquare, origin="lower")
            plt.colorbar()
            plt.xlabel('X [pixel]')
            plt.ylabel('Y [pixel]')
            plt.subplot(132)
            plt.title('Background')
            plt.imshow(backgroundSquare, origin="lower")
            plt.colorbar()
            plt.xlabel('X [pixel]')
            plt.ylabel('Y [pixel]')
            plt.subplot(133)
            plt.title('Difference')
            plt.imshow(peakSquare-backgroundSquare, origin="lower")
            plt.colorbar()
            plt.xlabel('X [pixel]')
            plt.ylabel('Y [pixel]')
            fig.tight_layout()
            plt.show()

    return [peakSet, peakSquareIntensitySum, backgroundSquareIntensitySum]
    


