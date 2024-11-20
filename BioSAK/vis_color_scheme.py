import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


vis_color_scheme_usage = '''
=================== vis_color_scheme example commands ===================

BioSAK vis_color_scheme -i lineage.txt -c tax_color.txt -o output.pdf
BioSAK vis_color_scheme -i Sponge_full_lineage_GTDB_format.txt -c color_code_sponge.txt -o color_code_sponge.pdf

=========================================================================
'''


def vis_color_scheme(args):

    tax_full_lineage_txt = args['i']
    tax_to_color_txt     = args['c']
    plot_width           = args['x']
    plot_height          = args['y']
    output_plot          = args['o']

    # read in taxon to color info
    taxon_to_color_dict = dict()
    for each_taxon in open(tax_to_color_txt):
        each_taxon_split = each_taxon.strip().split('\t')
        taxon_to_color_dict[each_taxon_split[0]] = each_taxon_split[1]

    max_lineage_len = 0
    full_lineage_list = []
    max_label_len = 0
    for each_lineage in open(tax_full_lineage_txt):
        each_lineage_split = each_lineage.strip().split(';')
        if len(each_lineage_split) > max_lineage_len:
                max_lineage_len = len(each_lineage_split)
        full_lineage_list.append(each_lineage_split)
        for each_rank in each_lineage_split:
                if len(each_rank) > max_label_len:
                    max_label_len = len(each_rank)

    fig = plt.figure(1, figsize=(plot_width, plot_height))
    ax = fig.add_subplot(111)

    y_pos = 0.97
    for each_lineage in full_lineage_list:
        x_pos = 0.005
        for each_rank in each_lineage:
            current_color = taxon_to_color_dict.get(each_rank, 'white')

            if current_color == 'white':
                label_txt = each_rank
                label_txt += ' '*((max_label_len - len(each_rank)) + 10)
            else:
                extra_space_str = ' '*(max_label_len - len(each_rank))
                label_txt = '%s %s%s' % (each_rank, extra_space_str, current_color)

            if current_color == 'white':
                ax.text(x_pos, y_pos, label_txt, color='black', family='monospace', bbox=dict(facecolor=current_color, edgecolor='black'))
            else:
                ax.text(x_pos, y_pos, label_txt, color='black', family='monospace', bbox=dict(facecolor=current_color, edgecolor=current_color))

            x_pos += (1/max_lineage_len)
        y_pos -= (1/(len(full_lineage_list)+2))

    # export plot
    plt.margins(tight=True)
    plt.savefig(output_plot, dpi=1000)
    plt.close()


if __name__ == '__main__':

    vis_color_scheme_parser = argparse.ArgumentParser(usage=vis_color_scheme_usage)
    vis_color_scheme_parser.add_argument('-i',  required=True,                          help='lineage file')
    vis_color_scheme_parser.add_argument('-c',  required=True,                          help='taxon color file')
    vis_color_scheme_parser.add_argument('-x',  required=False,type=int, default=42,    help='plot width, default is 42')
    vis_color_scheme_parser.add_argument('-y',  required=False,type=int, default=22,    help='plot height, default is 22')
    vis_color_scheme_parser.add_argument('-o',  required=True,                          help='output pdf')
    args = vars(vis_color_scheme_parser.parse_args())
    vis_color_scheme(args)
