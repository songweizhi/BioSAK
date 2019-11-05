import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def box_plotter(num_list_in, label_list, output_plot):

    # turn num list into arrary
    MAG_HGT_num_lol_arrary = [np.array(i) for i in num_list_in]

    # get plot
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(MAG_HGT_num_lol_arrary)

    # set x tick labels
    if len(label_list) <= 10:
        rotation_value = 0
    else:
        rotation_value = 270
    ax.set_xticklabels(label_list, rotation=rotation_value, fontsize=8)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='+', color='black', alpha=0.7, markersize=3)

    plt.tight_layout()
    fig.savefig(output_plot, bbox_inches='tight', dpi=300)
    plt.close()


num_list_in = [[1,2,3,4,5], [4,5,6,7,8], [6,7,8,9,10]]
label_list = ['small', 'medium', 'high']
output_plot = '/Users/songweizhi/Desktop/test.png'



box_plotter(num_list_in, label_list, output_plot)



