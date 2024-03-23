import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


boxplot_usage = '''
===================== boxplot example commands =====================

python3 boxplot_with_dots.py -i demo_boxplot.txt -o boxplot.pdf
python3 boxplot_with_dots.py -i demo_boxplot.txt -o boxplot.pdf -d

# format of input file (tab separated, no row name)
Apple	Banana	Orange
1	4	6
2	5	7
3	6	8
4	7	9

====================================================================
'''


def boxplot_with_dots(args):

    data_file   = args['i']
    add_dots    = args['d']
    output_plot = args['o']

    # read in file
    data_dict = dict()
    col_index = dict()
    line_num_index = 0
    for each_line in open(data_file):
        line_num_index += 1
        line_split = each_line.strip().split('\t')
        if line_num_index == 1:
            col_index = {key: i for key, i in enumerate(line_split)}
            data_dict = {i: [] for key, i in enumerate(line_split)}
        else:
            element_pos = 0
            for each_element in line_split:
                grp = col_index[element_pos]
                print('%s\t%s\t%s' % (grp, element_pos, each_element))
                data_dict[grp].append(float(each_element))
                element_pos += 1

    # get label
    name_list = sorted(list(data_dict.keys()))
    num_lol = [data_dict[i] for i in name_list]

    # plot
    data = [np.array(i) for i in num_lol]
    box_plot = plt.boxplot(data, labels=name_list, patch_artist=True, whiskerprops=dict(color='black', linewidth=1), capprops=dict(color='black'))

    # set the color pf box
    for box in box_plot['boxes']:
        box.set(linewidth=1)
        box.set_facecolor('white')

    # add dots, https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    if add_dots is True:
        col_index = 1
        for num_arrary in data:
            plt.plot(np.random.normal(col_index, 0.02, len(num_arrary)), num_arrary, '.', alpha=0.8, color='orange', markersize=6, markeredgewidth=0)
            col_index += 1

    # write out
    plt.tight_layout()
    plt.savefig(output_plot, bbox_inches='tight', dpi=300)
    plt.close()


if __name__ == '__main__':

    boxplot_parser = argparse.ArgumentParser()
    boxplot_parser.add_argument('-i',   required=True,                          help='input files')
    boxplot_parser.add_argument('-d',   required=False, action="store_true",    help='show dots')
    boxplot_parser.add_argument('-o',   required=True,                          help='output plot')
    args = vars(boxplot_parser.parse_args())
    boxplot_with_dots(args)
