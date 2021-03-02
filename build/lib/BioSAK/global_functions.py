import os
import shutil
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

time_format = '[%Y-%m-%d %H:%M:%S] '


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)

    os.mkdir(folder_to_create)


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)
    return folder_list


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def ctg_depth_and_gbk_to_gene_depth(ctg_depth_file, gbk_file, skip_depth_file_header, gene_depth_file_folder):

    gbk_file_path, gbk_file_basename, gbk_file_extension = sep_path_basename_ext(gbk_file)
    pwd_depth_file = '%s/%s.depth' % (gene_depth_file_folder, gbk_file_basename)

    # read in depth
    ctg_depth_dict = {}
    line = 0
    for ctg in open(ctg_depth_file):

        ctg_split = ctg.strip().split('\t')

        if skip_depth_file_header is True:
            if line > 0:
                ctg_depth_dict[ctg_split[0]] = float(ctg_split[1])
        else:
            ctg_depth_dict[ctg_split[0]] = float(ctg_split[1])

        line += 1

    # get gene depth
    gene_depth_file_handle = open(pwd_depth_file, 'w')
    gene_depth_file_handle.write('Gene\tDepth\n')
    for seq_record in SeqIO.parse(gbk_file, 'genbank'):

        seq_id = seq_record.id
        seq_depth = ctg_depth_dict[seq_id]

        for feature in seq_record.features:
            if feature.type == 'CDS':
                gene_id = feature.qualifiers['locus_tag'][0]
                for_out = '%s\t%s\n' % (gene_id, seq_depth)
                gene_depth_file_handle.write(for_out)

    gene_depth_file_handle.close()


def barh_plotter(num_list, label_list, query_seq_num, query_ko_NA, fig_width, fig_height, plot_file):

    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)

    y_pos = range(len(num_list))
    ax.barh(y_pos, num_list, height=0.8, align='center', alpha=0.2, linewidth=0)
    ax.set_yticks([])  # not show yticks
    ax.invert_xaxis()  # line up bar on right
    ax.invert_yaxis()  # put first number on top
    ax.axis('tight')   # remove extra spaces at the top and bottom, equal to: ax.margins(0, 0)
    # ax.margins(0, 0.01) # customize space percentage

    ax.set_xlabel('Number of gene')
    ax.set_title('Query genes number: %s, genes without KO: %s' % (query_seq_num, query_ko_NA))

    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(label_list)

    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()
    plt.clf()


def AnnotateNorm(file_in, skip_header, value_column, Divisor_value, file_out, file_out_header):

    file_out_handle = open(file_out, 'w')
    file_out_handle.write(file_out_header)
    line_num = 0
    for each_line in open(file_in):

        each_line_split = each_line.strip().split('\t')
        value_str = each_line_split[value_column - 1]

        if (skip_header is True and line_num > 0) or (skip_header is False):
            value_pct = float(value_str) * 100 / Divisor_value
            each_line_split[value_column - 1] = str(float("{0:.2f}".format(value_pct)))
            file_out_handle.write('%s\n' % '\t'.join(each_line_split))

        line_num += 1

    file_out_handle.close()


def get_gene_list_TotalDepth(gene_list, gene_to_depth_dict):

    total_depth = 0
    for gene in gene_list:
        gene_depth = gene_to_depth_dict[gene]
        total_depth += gene_depth

    return total_depth

