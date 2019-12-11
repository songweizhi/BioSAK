#  http://stackoverflow.com/questions/14391959/heatmap-in-matplotlib-with-pcolor

import matplotlib.pyplot as plt
import numpy as np

density = [[70.7, 68.95, 68.99, 69.49, 69.46, 69.45, 69.71, 68.64, 70.27, 69.26, 68.3, 69.39], [68.97, 71.09, 69.23, 69.01, 67.31, 69.58, 69.17, 71.04, 69.98, 69.24, 69.13, 69.29], [68.99, 69.24, 70.48, 69.3, 68.9, 69.07, 69.45, 69.57, 69.7, 69.55, 69.18, 69.37], [69.5, 69.02, 69.31, 70.23, 67.53, 69.33, 69.53, 69.3, 70.04, 69.4, 69.34, 69.5], [69.46, 67.31, 68.71, 67.53, 0, 68.9, 69.14, 69.2, 69.59, 68.87, 0, 70.41], [69.44, 69.6, 69.06, 69.31, 68.9, 71.45, 69.74, 69.78, 69.77, 69.93, 69.5, 69.43], [69.69, 69.16, 69.46, 69.52, 69.14, 69.74, 72.38, 69.38, 70.02, 72.77, 70.38, 69.52], [68.82, 71.04, 69.57, 69.3, 69.2, 69.77, 69.38, 0, 70.28, 69.46, 71.14, 69.41], [70.36, 69.96, 69.7, 70.04, 69.59, 69.77, 70.02, 70.28, 0, 69.86, 69.1, 69.89], [69.25, 69.24, 69.52, 69.41, 68.87, 69.92, 72.77, 69.47, 69.82, 72.96, 70.0, 69.39], [68.3, 69.13, 69.18, 69.37, 0, 69.64, 70.38, 71.14, 69.1, 69.93, 0, 69.7], [69.39, 69.26, 69.37, 69.49, 70.38, 69.44, 69.53, 69.42, 69.87, 69.38, 69.68, 69.88]]
data = np.array(density)

column_labels = list('ABCDEFGHIJKL')
row_labels = list('ABCDEFGHIJKL')

fig, ax = plt.subplots()
heatmap = ax.pcolor(data)
plt.colorbar(heatmap)

ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False) # put the major ticks at the middle of each cell
ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False) # put the major ticks at the middle of each cell
ax.invert_yaxis() #setup the direction of Y-axis
ax.xaxis.tick_top() #put X-axis on the top of the figure
ax.set_xticklabels(row_labels, minor=False)
ax.set_yticklabels(column_labels, minor=False)

plt.show()


