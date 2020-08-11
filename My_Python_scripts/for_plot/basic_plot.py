import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


########################################################################################################################

def get_bar_plot(num_list_1, num_list_2, label_list, label_rotation, output_plot):

    # Set position of bar on X axis
    barWidth = 0.4
    r1 = np.arange(len(num_list_1))
    r2 = [x + barWidth for x in r1]

    plt.bar(r1, num_list_1, color='orange', width=barWidth, edgecolor='white', label='Number')
    plt.bar(r2, num_list_2, color='lightgreen', width=barWidth, edgecolor='white', label='Size')
    plt.xticks([r + barWidth for r in range(len(num_list_1))], label_list, rotation=label_rotation, fontsize=10)
    plt.title('The number and total size of bins')
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()


num_list_1 = [65, 29, 5, 99, 126, 111, 8, 60]
num_list_2 = [247, 86, 10, 370, 454, 439, 17, 172]
label_list = ['Apl', 'Car', 'Cli', 'Cos', 'Irc', 'Rho', 'Sty', 'swt']
label_rotation = 0
output_plot = '/Users/songweizhi/Desktop/test.png'
get_bar_plot(num_list_1, num_list_2, label_list, label_rotation, output_plot)


########################################################################################################################

def get_boxplot(num_lol, label_list, label_rotation, output_plot):

    num_lol_arrary = [np.array(i) for i in num_lol]
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(num_lol_arrary)
    ax.set_xticklabels(label_list, rotation=label_rotation, fontsize=8)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='+', color='black', alpha=0.7, markersize=3)

    plt.tight_layout()
    fig.savefig(output_plot, bbox_inches='tight', dpi=300)
    plt.close()


num_lol = [[1,2,3], [4,5,6], [7,8,9]]
label_list = ['Apl', 'Car', 'Cli']
label_rotation = 0
output_plot = '/Users/songweizhi/Desktop/boxplot.png'
get_boxplot(num_lol, label_list, label_rotation, output_plot)


########################################################################################################################

def get_categorical_scatter_plot(num_list, label_list, label_rotation, output_plot):

    num_index = list(range(1, len(num_list) + 1))

    plt.scatter(num_index, num_list)
    plt.xticks(num_index, label_list, rotation=label_rotation)
    plt.margins(0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()


num_list = [1, 2, 4, 3]
label_list = ['Frogs', 'Hogs', 'Bogs', 'Slogs']
get_categorical_scatter_plot(num_list, label_list, 0, '/Users/songweizhi/Desktop/categorical_scatter_plot.png')


# # add vertical lines
# for x_16s in ribosomal_RNA_16S_pos_list_dict['2.10_chromosome']:
#     plt.axvline(x=x_16s)


# # add horizontal lines
# for y_16s in ribosomal_RNA_16S_pos_list_dict['D2_c']:
#     plt.hlines(y=y_16s, xmin=0, xmax=ref_length_dict['2.10_chromosome'])
