import numpy as np
from matplotlib.mlab import PCA

dataMatrix = np.array([[1.0, 0.01, 0.01], [2.0, -0.02, -0.01], [3.0, 0.01, 0.02], [4.0, -0.03, -0.02], [5.0, 0.02, 0], [6.0, 0.04, 0.01], [7.0, -0.02, 0.03]])
myPCA = PCA(dataMatrix)

print('PCA fracs:\n', myPCA.fracs)