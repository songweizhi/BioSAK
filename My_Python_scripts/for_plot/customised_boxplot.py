import numpy as np
import matplotlib.pyplot as plt


a = np.random.normal(0, 2, 20)
b = np.random.normal(-2, 7, 20)
data = [a, b]

box_plot = plt.boxplot(data,
                       patch_artist=True,
                       whiskerprops=dict(color='lightblue', linewidth=2),
                       capprops=dict(color='lightblue'))


# set the color pf box
for box in box_plot['boxes']:
    box.set(linewidth=0)
    box.set_facecolor('lightblue')


# add dots, https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.plot.html
plt.plot(np.random.normal(1, 0.02, len(a)), a, '.', alpha=0.8, marker='o', color='orange', markersize=6, markeredgewidth=0)
plt.plot(np.random.normal(2, 0.02, len(b)), b, '.', alpha=0.8, marker='o', color='orange', markersize=6, markeredgewidth=0)


# for i in [1,2]:
#     y = data[i-1]
#     x = np.random.normal(i, 0.02, len(y))
#     plt.plot(x, y, 'r.', alpha=0.2)

plt.show()
