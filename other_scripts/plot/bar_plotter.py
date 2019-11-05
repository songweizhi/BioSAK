import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def bar_plotter(label_list, num_list, output_plot, y_axis_label):

    x_range = range(len(label_list))
    plt.bar(x_range, num_list, tick_label=label_list, align='center', alpha=0.2, linewidth=0)
    plt.xticks(x_range, label_list, rotation=270, fontsize=10, horizontalalignment='center')
    plt.title('')
    plt.ylabel(y_axis_label)
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()

