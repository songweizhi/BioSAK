import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def boxplot_with_dots(num_lol, name_list, output_plot):

    data = [np.array(i) for i in num_lol]

    box_plot = plt.boxplot(data, labels=name_list, patch_artist=True,
                           whiskerprops=dict(color='lightblue', linewidth=2), capprops=dict(color='lightblue'))

    # set the color pf box
    for box in box_plot['boxes']:
        box.set(linewidth=0)
        box.set_facecolor('lightblue')

    # add dots, https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    col_index = 1
    for num_arrary in data:
        plt.plot(np.random.normal(col_index, 0.02, len(num_arrary)), num_arrary, '.', alpha=0.8, color='orange',
                 markersize=6, markeredgewidth=0)
        col_index += 1

    # write out
    plt.tight_layout()
    plt.savefig(output_plot, bbox_inches='tight', dpi=300)
    plt.close()


num_lol     = [[1, 2, 3, 4], [4, 5, 6, 7], [6, 7, 8, 9]]
name_list   = ['apple', 'banana', 'orange']
output_plot = '/Users/songweizhi/Desktop/boxplot.png'

boxplot_with_dots(num_lol, name_list, output_plot)
