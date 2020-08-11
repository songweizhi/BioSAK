import argparse
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec


pwd_figure = '/Users/songweizhi/Desktop/multi_subplots.png'

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 9], width_ratios=[9, 1])
ax0 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[2])
ax3 = plt.subplot(gs[3])


# ax0
ax0.bar([1,2,3,4,5,6,7,8,9], [4,2,3,4,5,2,6,3,5], color='skyblue', width=1, linewidth=0)
plt.xlim(0, 10)
plt.ylim(0, 10)
# hide x-axis
#x_axis = ax0.axes.get_xaxis()
#x_axis.set_visible(False)

# ax2
ax2.scatter([1,2,3,4,5,6,7,8,9], [1,2,3,4,5,6,7,8,9], s=32, linewidths=0, marker='o', c='orange')
ax2.scatter([9,8,7,6,5,4,3,2,1], [1,2,3,4,5,6,7,8,9], s=32, linewidths=0, marker='o', c='lightblue')
plt.xlim(0, 10)
plt.ylim(0, 10)

# ax3
ax3.barh([1,2,3,4,5,6,7,8,9], [4,2,3,4,5,2,6,3,5], color='skyblue', linewidth=0)
plt.ylim(0, 10)
plt.xlim(0, 10)
# hide x-axis
#y_axis = ax3.axes.get_yaxis()
#y_axis.set_visible(False)


plt.tight_layout()
plt.savefig(pwd_figure)
plt.close()

