import os
import sys
import argparse


iTOL_usage = '''
=========================================== iTOL example commands ===========================================

BioSAK iTOL -ColorStrip -LeafGroup LeafGroup.txt -GroupColor GroupColor.txt -prefix NorthSea -LegendTitle Phylum 
BioSAK iTOL -ColorRange -LeafGroup LeafGroup.txt -GroupColor GroupColor.txt -prefix NorthSea -LegendTitle Phylum 

# LeafGroup file format (tab separated)
NorthSea_bin070	Alphaproteobacteria
NorthSea_bin054	Actinobacteria
NorthSea_bin103	Verrucomicrobiae

# GroupColor file format (tab separated)
Alphaproteobacteria	#CCCC00
Verrucomicrobiae	#9999FF
Actinobacteria	#CC6600

=============================================================================================================
'''

'''

cd /Users/songweizhi/Desktop/iTOL
python ~/PycharmProjects/BioSAK/BioSAK/iTOL.py -ColorStrip -LeafGroup NorthSea_LeafGroup.txt -GroupColor NorthSea_GroupColor.txt -prefix NorthSea -LegendTitle MAG_Class 
python ~/PycharmProjects/BioSAK/BioSAK/iTOL.py -ColorRange -LeafGroup NorthSea_LeafGroup.txt -GroupColor NorthSea_GroupColor.txt -prefix NorthSea -LegendTitle MAG_Class 

'''


def iTOL(args):

    # read in arguments
    FilePrefix      = args['prefix']
    ColorStrip      = args['ColorStrip']
    ColorRange      = args['ColorRange']
    LeafGroup       = args['LeafGroup']
    GroupColor      = args['GroupColor']
    LegendTitle     = args['LegendTitle']
    STRIP_WIDTH     = 100
    MARGIN          = 20

    if (ColorStrip is True) and (ColorRange is True):
        print('Please specify one at a time')
        exit()

    elif (ColorStrip is True) or (ColorRange is True):

        file_out = ''
        if ColorStrip is True:
            file_out = '%s_ColorStrip.txt' % FilePrefix
        if ColorRange is True:
            file_out = '%s_ColorRange.txt' % FilePrefix

        Leaf_to_Group_dict = {}
        for each_leaf in open(LeafGroup):
            each_leaf_split = each_leaf.strip().split('\t')
            Leaf_to_Group_dict[each_leaf_split[0]] = each_leaf_split[1]

        Group_to_Color_dict = {}
        for each_group in open(GroupColor):
            each_group_split = each_group.strip().split('\t')
            Group_to_Color_dict[each_group_split[0]] = each_group_split[1]

        group_list = [i for i in Group_to_Color_dict]
        color_list = [Group_to_Color_dict[i] for i in group_list]

        file_out_handle = open(file_out, 'w')

        # write out header
        if ColorStrip is True:
            file_out_handle.write('DATASET_COLORSTRIP\n')
            file_out_handle.write('SEPARATOR TAB\n')
            file_out_handle.write('DATASET_LABEL\t%s_ColorStrip\n' % LegendTitle)
        if ColorRange is True:
            file_out_handle.write('TREE_COLORS\n')
            file_out_handle.write('SEPARATOR TAB\n')
            file_out_handle.write('DATASET_LABEL\t%s_ColorRange\n' % LegendTitle)

        # write out strip attributes
        if ColorStrip is True:
            file_out_handle.write('\n# customize strip attributes here\n')
            file_out_handle.write('STRIP_WIDTH\t%s\n' % STRIP_WIDTH)
            file_out_handle.write('MARGIN\t%s\n'      % MARGIN)

        # write out legend info
        file_out_handle.write('\n# customize legend here\n')
        file_out_handle.write('LEGEND_TITLE\t%s\n' % LegendTitle)
        file_out_handle.write('LEGEND_SHAPES\t%s\n' % '\t'.join(['1' for i in Group_to_Color_dict]))
        file_out_handle.write('LEGEND_COLORS\t%s\n' % '\t'.join(color_list))
        file_out_handle.write('LEGEND_LABELS\t%s\n' % '\t'.join(group_list))

        # write out data info
        file_out_handle.write('\n# provide data here\nDATA\n')
        for leaf in Leaf_to_Group_dict:
            leaf_group = Leaf_to_Group_dict[leaf]
            leaf_color = Group_to_Color_dict[leaf_group]

            if ColorStrip is True:
                file_out_handle.write('%s\t%s\t%s\n' % (leaf, leaf_color, leaf_group))
            if ColorRange is True:
                file_out_handle.write('%s\trange\t%s\t%s\n' % (leaf, leaf_color, leaf_group))

        file_out_handle.close()


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-prefix',      required=True,                         help='Output prefix')
    parser.add_argument('-ColorStrip',  required=False, action='store_true',   help='ColorStrip')
    parser.add_argument('-ColorRange',  required=False, action='store_true',   help='ColorRange')
    parser.add_argument('-LeafGroup',   required=False, default=None,          help='Leaf Group')
    parser.add_argument('-GroupColor',  required=False, default=None,          help='Group Color')
    parser.add_argument('-LegendTitle', required=False, default=None,          help='Legend Title')

    args = vars(parser.parse_args())

    iTOL(args)
