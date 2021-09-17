from itertools import combinations
import numpy as np
L = np.array([[1, 2, 3, 4], [1, 2, 3, 4]])

# https://stackoverflow.com/questions/32208359/is-there-a-multi-dimensional-version-of-arange-linspace-in-numpy
#np.arange equivalent
xy = np.mgrid[-5:5.1:0.5, -5:5.1:0.5].reshape(2,-1).T
print(xy)

#np.linspace equivalent
xy = np.mgrid[-5:5:21j, -5:5:21j].reshape(2,-1).T
print(xy)

#newshapeint or tuple of ints
#The new shape should be compatible with the original shape.
# If an integer, then the result will be a 1-D array of that length.
# One shape dimension can be -1.
# In this case, the value is inferred from the length of the array and remaining dimensions.
L = np.array([[1, 2, 3, 4], [1, 2, 3, 4]])
test = np.mgrid[0:L.shape[0]+.1:1, 0:L.shape[1]+.1:1].reshape(len(L.shape),-1).T
print(test)
