import math
import random
import argparse
import seaborn as sns


iTOL_usage = '''
==================================== iTOL example commands ====================================

# Example commands
BioSAK iTOL -ColorStrip -lg MagTaxon.txt -lt Phylum -out ColorStrip_taxon.txt
BioSAK iTOL -ColorRange -lg MagTaxon.txt -lt Phylum -out ColorRange_taxon.txt
BioSAK iTOL -SimpleBar -lv MagSize.txt -scale 0-3-6-9 -lt Size -out SimpleBar_size.txt
BioSAK iTOL -Heatmap -lm MagAbundance.txt -lt Abundance -out Heatmap_abundance.txt

# Leaf to Group file format (-lg, tab separated, no header)
NorthSea_bin001	Alphaproteobacteria
NorthSea_bin003	Verrucomicrobiae

# Group to Color file format (-gc, tab separated, no header)
Alphaproteobacteria	#CCCC00
Verrucomicrobiae	#9999FF

# Leaf to Value file format (-lv, tab separated, no header)
NorthSea_bin001	6.15
NorthSea_bin002	6.63

# Leaf to Matrix file format (-lm, tab separated, header required!!!)
MAG_id Sample_A   Sample_B   Sample_C
NorthSea_bin001	6.15    2.23    1.56
NorthSea_bin002	6.63    1.72    2.55

# Note!!!
The purpose for developing this module is to generate iTOL-recognizable file for your dataset, 
parameters (e.g. colour, font size and strip width) provided in the output file might need 
further adjustment. The best way to optimize your plot is to visualize your tree with the 
default output and optimize the parameters by looking at the tree.

===============================================================================================
'''


def get_color_list(color_num):

    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']

    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a', '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']

    else:
        color_num_each = math.ceil(color_num/8) + 2

        color_list_1 = sns.color_palette('Blues',  n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy',   n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green',  n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds',   n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive',  n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()

        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6, color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    return color_list_to_return_sorted


def iTOL(args):

    # read in arguments
    ColorStrip      = args['ColorStrip']
    ColorRange      = args['ColorRange']
    SimpleBar       = args['SimpleBar']
    Heatmap         = args['Heatmap']
    LeafGroup       = args['lg']
    GroupColor      = args['gc']
    LeafValue       = args['lv']
    LeafMatrix      = args['lm']
    scale_str       = args['scale']
    LegendTitle     = args['lt']
    FileOut         = args['out']

    # General
    STRIP_WIDTH                 = 100
    MARGIN                      = 20

    # SimpleBar
    SimpleBar_COLOR             = 'grey'
    SimpleBar_WIDTH             = 300
    SimpleBar_HEIGHT_FACTOR     = 0.8
    SimpleBar_BORDER_WIDTH      = 0
    SimpleBar_SCALE_COLOR       = '#696969'
    SimpleBar_SCALE_WIDTH       = 1
    SimpleBar_SCALE_DASHED      = 1
    SimpleBar_SCALE_FontSize    = 2

    # Heatmap
    Heatmap_STRIP_WIDTH         = 60

    # check the number of specified file type
    True_num = 0
    for file_type in [ColorStrip, ColorRange, SimpleBar, Heatmap]:
        if file_type is True:
            True_num += 1

    if True_num == 0:
        print('Please specify one file type, choose from -ColorStrip, -ColorRange, -SimpleBar or -Heatmap')
        exit()
    if True_num > 1:
        print('Please specify one file type ONLY, choose from -ColorStrip, -ColorRange, -SimpleBar or -Heatmap')
        exit()

    # Prepare ColorStrip and ColorRange file
    if (ColorStrip is True) or (ColorRange is True):

        Leaf_to_Group_dict = {}
        Group_list = []
        for each_leaf in open(LeafGroup):
            each_leaf_split = each_leaf.strip().split('\t')
            Leaf_to_Group_dict[each_leaf_split[0]] = each_leaf_split[1]
            # get Group_list
            if each_leaf_split[1] not in Group_list:
                Group_list.append(each_leaf_split[1])

        Group_to_Color_dict = {}
        if GroupColor is None:
            color_list = get_color_list(len(Group_list))
            Group_to_Color_dict = dict(zip(Group_list, color_list))
        else:
            for each_group in open(GroupColor):
                each_group_split = each_group.strip().split('\t')
                Group_to_Color_dict[each_group_split[0]] = each_group_split[1]

        group_list = [i for i in Group_to_Color_dict]
        color_list = [Group_to_Color_dict[i] for i in group_list]

        FileOut_handle = open(FileOut, 'w')

        # write out header
        if ColorStrip is True:
            FileOut_handle.write('DATASET_COLORSTRIP\n')
            FileOut_handle.write('SEPARATOR TAB\n')
            FileOut_handle.write('DATASET_LABEL\t%s_ColorStrip\n' % LegendTitle)
        if ColorRange is True:
            FileOut_handle.write('TREE_COLORS\n')
            FileOut_handle.write('SEPARATOR TAB\n')
            FileOut_handle.write('DATASET_LABEL\t%s_ColorRange\n' % LegendTitle)

        # write out strip attributes
        if ColorStrip is True:
            FileOut_handle.write('\n# customize strip attributes here\n')
            FileOut_handle.write('STRIP_WIDTH\t%s\n' % STRIP_WIDTH)
            FileOut_handle.write('MARGIN\t%s\n'      % MARGIN)

        # write out legend info
        FileOut_handle.write('\n# customize legend here\n')
        FileOut_handle.write('LEGEND_TITLE\t%s\n' % LegendTitle)
        FileOut_handle.write('LEGEND_SHAPES\t%s\n' % '\t'.join(['1' for i in Group_to_Color_dict]))
        FileOut_handle.write('LEGEND_COLORS\t%s\n' % '\t'.join(color_list))
        FileOut_handle.write('LEGEND_LABELS\t%s\n' % '\t'.join(group_list))

        # write out data info
        FileOut_handle.write('\n# provide data here\nDATA\n')
        for leaf in Leaf_to_Group_dict:
            leaf_group = Leaf_to_Group_dict[leaf]
            leaf_color = Group_to_Color_dict[leaf_group]

            if ColorStrip is True:
                FileOut_handle.write('%s\t%s\t%s\n' % (leaf, leaf_color, leaf_group))
            if ColorRange is True:
                FileOut_handle.write('%s\trange\t%s\t%s\n' % (leaf, leaf_color, leaf_group))

        FileOut_handle.close()


    # Prepare SimpleBar file
    if SimpleBar is True:

        if scale_str is None:
            print('Please provide scale values for barchart, in format 0-3-6-9')
            exit()

        # read in leaf value into dict
        leaf_value_dict = {}
        max_value = None
        for leaf_value in open(LeafValue):
            leaf_value_split = leaf_value.strip().split('\t')
            leaf_value_dict[leaf_value_split[0]] = float(leaf_value_split[1])

            # get max value
            if max_value == None:
                max_value = float(leaf_value_split[1])
            else:
                if float(leaf_value_split[1]) > max_value:
                    max_value = float(leaf_value_split[1])

        SimpleBar_FileOut_handle = open(FileOut, 'w')

        # write out header
        SimpleBar_FileOut_handle.write('DATASET_SIMPLEBAR\n')
        SimpleBar_FileOut_handle.write('# Reference: https://itol.embl.de/help/dataset_simplebar_template.txt\n')
        SimpleBar_FileOut_handle.write('\nSEPARATOR TAB\n')

        # write out SimpleBar attributes
        SimpleBar_FileOut_handle.write('\n# customize barchart attributes here\n')
        SimpleBar_FileOut_handle.write('DATASET_LABEL\t%s\n'    % LegendTitle)
        SimpleBar_FileOut_handle.write('COLOR\t%s\n'            % SimpleBar_COLOR)
        SimpleBar_FileOut_handle.write('WIDTH\t%s\n'            % SimpleBar_WIDTH)
        SimpleBar_FileOut_handle.write('MARGIN\t%s\n'           % MARGIN)
        SimpleBar_FileOut_handle.write('HEIGHT_FACTOR\t%s\n'    % SimpleBar_HEIGHT_FACTOR)
        SimpleBar_FileOut_handle.write('BORDER_WIDTH\t%s\n'     % SimpleBar_BORDER_WIDTH)

        # write out scale attributes
        scale_attributes_list = []
        for scale_value in scale_str.split('-'):
            scale_attributes = '%s-%s-%s-%s-%s-%s' % (scale_value, scale_value, SimpleBar_SCALE_COLOR, SimpleBar_SCALE_WIDTH, SimpleBar_SCALE_DASHED, SimpleBar_SCALE_FontSize)
            scale_attributes_list.append(scale_attributes)

        SimpleBar_FileOut_handle.write('\n# customize scale attributes here\n')
        SimpleBar_FileOut_handle.write('# format: VALUE-LABEL-COLOR-WIDTH-DASHED-LABEL_SCALE_FACTOR, LABEL_SCALE_FACTOR controls font size\n')
        SimpleBar_FileOut_handle.write('DATASET_SCALE\t%s\n' % '\t'.join(scale_attributes_list))

        # write out data info
        SimpleBar_FileOut_handle.write('\n# provide data here\n')
        SimpleBar_FileOut_handle.write('DATA\n')

        for leaf in leaf_value_dict:
            SimpleBar_FileOut_handle.write('%s\t%s\n' % (leaf, leaf_value_dict[leaf]))

        SimpleBar_FileOut_handle.close()

    # Prepare Heatmap file
    if Heatmap is True:

        # read in leaf matrix into dict
        n = 0
        col_name_list = []
        leaf_matrix_dict = {}
        for leaf_matrix in open(LeafMatrix):
            leaf_matrix_split = leaf_matrix.strip().split('\t')
            if n == 0:
                col_name_list = leaf_matrix_split[1:]
            else:
                leaf_matrix_dict[leaf_matrix_split[0]] = leaf_matrix_split[1:]
            n += 1

        Heatmap_FileOut_handle = open(FileOut, 'w')

        # write out header
        Heatmap_FileOut_handle.write('DATASET_HEATMAP\n')
        Heatmap_FileOut_handle.write('# Reference https://itol.embl.de/help/dataset_heatmap_template.txt\n')
        Heatmap_FileOut_handle.write('\nSEPARATOR TAB\n')

        # write out heatmap attributes
        Heatmap_FileOut_handle.write('\n# customize heatmap attributes here\n')
        Heatmap_FileOut_handle.write('MARGIN\t%s\n'          % MARGIN)
        Heatmap_FileOut_handle.write('STRIP_WIDTH\t%s\n'     % Heatmap_STRIP_WIDTH)

        # write out legend info
        Heatmap_FileOut_handle.write('\n# customize legend here\n')
        Heatmap_FileOut_handle.write('AUTO_LEGEND\t1\n')
        Heatmap_FileOut_handle.write('DATASET_LABEL\t%s\n' % LegendTitle)
        Heatmap_FileOut_handle.write('USE_MID_COLOR\t1\n')
        Heatmap_FileOut_handle.write('COLOR_MIN\t#2980B9\n')
        Heatmap_FileOut_handle.write('COLOR_MID\t#ECF0F1\n')
        Heatmap_FileOut_handle.write('COLOR_MAX\t#E74C3C\n')
        Heatmap_FileOut_handle.write('\n# customize value range here. By default, color gradients will be calculated based on dataset values\n')
        Heatmap_FileOut_handle.write('# USER_MIN_VALUE	0\n')
        Heatmap_FileOut_handle.write('# USER_MID_VALUE	5\n')
        Heatmap_FileOut_handle.write('# USER_MAX_VALUE	10\n')
        Heatmap_FileOut_handle.write('\n# customize column name here\n')
        Heatmap_FileOut_handle.write('FIELD_LABELS\t%s\n' % '\t'.join(col_name_list))

        # write out data info
        Heatmap_FileOut_handle.write('\n# Provide data here\n')
        Heatmap_FileOut_handle.write('DATA\n')
        for leaf in leaf_matrix_dict:
            Heatmap_FileOut_handle.write('%s\t%s\n' % (leaf, '\t'.join(leaf_matrix_dict[leaf])))

        Heatmap_FileOut_handle.close()


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-ColorStrip',  required=False, action='store_true',   help='ColorStrip')
    parser.add_argument('-ColorRange',  required=False, action='store_true',   help='ColorRange')
    parser.add_argument('-SimpleBar',   required=False, action='store_true',   help='SimpleBar')
    parser.add_argument('-Heatmap',     required=False, action='store_true',   help='Heatmap')
    parser.add_argument('-lg',          required=False, default=None,          help='Leaf Group')
    parser.add_argument('-gc',          required=False, default=None,          help='Specify Group Color (optional)')
    parser.add_argument('-lv',          required=False, default=None,          help='Leaf Value')
    parser.add_argument('-lm',          required=False, default=None,          help='Leaf Matrix')
    parser.add_argument('-scale',       required=False, default=None,          help='Scale Values, in format 0-3-6-9')
    parser.add_argument('-lt',          required=False, default=None,          help='Legend Title')
    parser.add_argument('-out',         required=True,                         help='Output filename')

    args = vars(parser.parse_args())

    iTOL(args)

