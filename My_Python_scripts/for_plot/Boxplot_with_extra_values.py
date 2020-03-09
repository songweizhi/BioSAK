import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def get_boxplot(num_lol, label_list, label_rotation, output_plot):

    num_lol_arrary = [np.array(i) for i in num_lol]
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(num_lol_arrary)
    ax.set_xticklabels(label_list, rotation=label_rotation, fontsize=8)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='+', color='black', alpha=0.7, markersize=3)

    # add dots, https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    plt.plot(1, 7, alpha=1, marker='^', markersize=10, markeredgewidth=0, color='coral')
    plt.plot(2, 2, alpha=1, marker='v', markersize=10, markeredgewidth=0, color='deepskyblue')
    plt.plot(3, 4, alpha=1, marker='v', markersize=10, markeredgewidth=0, color='deepskyblue')

    plt.tight_layout()
    fig.savefig(output_plot, bbox_inches='tight', dpi=300)
    plt.close()


num_lol = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
label_list = ['Apl', 'Car', 'Cli']
label_rotation = 0
output_plot = '/Users/songweizhi/Desktop/Boxplot_with_extra_values.png'
get_boxplot(num_lol, label_list, label_rotation, output_plot)


