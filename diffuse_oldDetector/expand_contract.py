'''
Tyler Chase
07/20/2016
take in a centered diffraction pattern as well as a factor and expand or contract 
the image by the factor using expand() to account for lattice expansion 
'''

from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

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

    
def main():
    # Form a test image for expansion
    x = np.array(range(100))
    y = np.array(range(100))
    z = np.zeros((100,100))
    z[30, 30] = 1
    z[69, 69] = 1
    z[30, 69] = 1
    z[69, 30] = 1

    # Expand the image by 2, then contract the image by 2, then expand and 
    # contract the image by 2
    factor = 2
    z_2 = expand(z, factor)
    z_3 = expand(z, 1.0/factor)
    z_4 = expand(expand(z, factor),1.0/factor)

    # Plot the initial image, the expanded image, the contracted image,
    # Then the expanded and contracted image
    plt.figure()
    plt.subplot(221)
    plt.imshow(z, interpolation = 'none')
    plt.title('Original')
    plt.subplot(222)
    plt.imshow(z_2, interpolation = 'none')
    plt.title('Contracted')
    plt.subplot(223)
    plt.imshow(z_3, interpolation = 'none')
    plt.title('Expanded')
    plt.subplot(224)
    plt.imshow(z_4, interpolation = 'none')
    plt.title('Contracted then Expanded Back')
    plt.show()
if __name__ == "__main__" :
    main()