import os
import glob
import argparse
from Bio import SeqIO

# Make ribbon diagrams of conserved linkages between genomes
# https://github.com/conchoecia/odp


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


ribbon_usage = '''
==================== ribbon example commands ====================

BioSAK ribbon -h

python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/ribbon.py -o demo -fa example_files -f -odp /scratch/PI/ocessongwz/Software/odp/scripts/odp -odp_rbh_to_ribbon /scratch/PI/ocessongwz/Software/odp/scripts/odp_rbh_to_ribbon

python3 ribbon.py -co chrom_order2.txt -o demo2 -fa example_files -pep example_files -chrom example_files -f -odp /scratch/PI/ocessongwz/Software/odp/scripts/odp -odp_rbh_to_ribbon /scratch/PI/ocessongwz/Software/odp/scripts/odp_rbh_to_ribbon
python3 ribbon.py -co chrom_order2.txt -o demo4 -fa example_files -pep example_files -chrom example_files -f -odp /scratch/PI/ocessongwz/Software/odp/scripts/odp -odp_rbh_to_ribbon /scratch/PI/ocessongwz/Software/odp/scripts/odp_rbh_to_ribbon

Note:
1. Species names can't have '_' char

=================================================================
'''


def ribbon(args):

    fa_in                   = args['fa']
    fa_x                    = args['fax']
    gff_in                  = args['gff']
    gff_x                   = args['gffx']
    gbk_in                  = args['gbk']
    gbk_x                   = args['gbkx']
    pep_in                  = args['pep']
    pep_x                   = args['pepx']
    chrom_in                = args['chrom']
    chrom_x                 = args['chromx']
    op_dir                  = args['o']
    species_order_txt       = args['so']
    chrom_order_txt         = args['co']
    thread_num              = args['t']
    force_overwrite         = args['f']
    minscafsize             = args['m']
    odp_exe                 = args['odp']
    odp_rbh_to_ribbon_exe   = args['odp_rbh_to_ribbon']
    plot_all                = args['plot_all']

    ################################################# define file name #################################################

    get_macrosynteny_plot_wd    = '%s/get_macrosynteny_plot'                % op_dir
    get_ribbon_diagram_wd       = '%s/get_ribbon_diagram'                   % op_dir
    config_yaml_macrosynteny    = '%s/config.yaml'                          % get_macrosynteny_plot_wd
    config_yaml_ribbon          = '%s/config.yaml'                          % get_ribbon_diagram_wd

    ####################################################################################################################

    fa_f_re   = '%s/*.%s' % (fa_in, fa_x)
    fa_f_list = glob.glob(fa_f_re)

    species_to_chrom_dict = dict()
    for each_fa in fa_f_list:
        _, _, fa_base, _ = sep_path_basename_ext(each_fa)
        qualified_ctg_set = []
        for each_seq in SeqIO.parse(each_fa, 'fasta'):
            seq_len = len(each_seq.seq)
            if seq_len >= minscafsize:
                qualified_ctg_set.append(each_seq.id)
        species_to_chrom_dict[fa_base] = qualified_ctg_set

    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()

    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % get_macrosynteny_plot_wd)
    os.system('mkdir %s' % get_ribbon_diagram_wd)

    ####################################################################################################################

    species_order_list = []
    if species_order_txt is not None:
        for each_s in open(species_order_txt):
            species_order_list.append(each_s.strip())
    else:
        species_order_list = sorted(list(species_to_chrom_dict.keys()))

    chrom_order_dict = dict()
    if chrom_order_txt is not None:
        for each_c in open(chrom_order_txt):
            each_c_split = each_c.strip().split('\t')
            s_id = each_c_split[0]
            c_id = each_c_split[1]
            if s_id not in chrom_order_dict:
                chrom_order_dict[s_id] = []
            chrom_order_dict[s_id].append(c_id)

    fa_in_abspath                    = os.path.abspath(fa_in)
    op_dir_abspath                   = os.path.abspath(op_dir)
    get_macrosynteny_plot_wd_abspath = os.path.abspath(get_macrosynteny_plot_wd)
    get_ribbon_diagram_wd_abspath    = os.path.abspath(get_ribbon_diagram_wd)

    rbh_file_dir                = '%s/odp/step2-figures/synteny_nocolor'    % get_macrosynteny_plot_wd_abspath

    ################################### get config.yaml for getting macrosynteny plot ##################################

    # get config.yaml for getting macrosynteny plot
    config_yaml_macrosynteny_handle = open(config_yaml_macrosynteny, 'w')
    config_yaml_macrosynteny_handle.write('ignore_autobreaks: True\n')
    config_yaml_macrosynteny_handle.write('diamond_or_blastp: "diamond"\n')
    config_yaml_macrosynteny_handle.write('duplicate_proteins: "pass"\n')
    config_yaml_macrosynteny_handle.write('plot_LGs: False\n')
    config_yaml_macrosynteny_handle.write('plot_sp_sp: True\n')
    config_yaml_macrosynteny_handle.write('\n')

    # write out species section
    config_yaml_macrosynteny_handle.write('species:\n')
    for each_species in species_order_list:
        current_fa    = '%s/%s.fa'    % (fa_in_abspath, each_species)
        current_pep   = '%s/%s.pep'   % (fa_in_abspath, each_species)
        current_chrom = '%s/%s.chrom' % (fa_in_abspath, each_species)
        config_yaml_macrosynteny_handle.write('  %s:\n' % each_species)
        config_yaml_macrosynteny_handle.write('    proteins:    %s\n' % current_pep)
        config_yaml_macrosynteny_handle.write('    chrom:       %s\n' % current_chrom)
        config_yaml_macrosynteny_handle.write('    genome:      %s\n' % current_fa)
        config_yaml_macrosynteny_handle.write('    minscafsize: %s\n' % minscafsize)
    config_yaml_macrosynteny_handle.close()

    #################################### get config.yaml for getting ribbon diagram ####################################

    config_yaml_ribbon_handle = open(config_yaml_ribbon, 'w')
    config_yaml_ribbon_handle.write('chr_sort_order: optimal-chr-or\n')
    if plot_all is True:
        config_yaml_ribbon_handle.write('plot_all: True\n')
    else:
        config_yaml_ribbon_handle.write('plot_all: False\n')
    config_yaml_ribbon_handle.write('\n')

    # write out species_order section
    config_yaml_ribbon_handle.write('species_order:\n')
    for species in species_order_list:
        config_yaml_ribbon_handle.write('  - %s\n' % species)
    config_yaml_ribbon_handle.write('\n')

    # write out chromorder section
    config_yaml_ribbon_handle.write('chromorder:\n')

    if chrom_order_txt is not None:
        for each_s in chrom_order_dict:
            config_yaml_ribbon_handle.write('  %s:\n' % each_s)
            for each_chrom in chrom_order_dict[each_s]:
                config_yaml_ribbon_handle.write('    - %s\n' % each_chrom)
        config_yaml_ribbon_handle.write('\n')

    else:
        for each_species in species_order_list:
            current_chrom_list = species_to_chrom_dict.get(each_species, [])
            if len(current_chrom_list) > 0:
                config_yaml_ribbon_handle.write('  %s:\n' % each_species)
                for each_chrom in current_chrom_list:
                    config_yaml_ribbon_handle.write('    - %s\n' % each_chrom)
            config_yaml_ribbon_handle.write('\n')

    # write out rbh_directory section
    config_yaml_ribbon_handle.write('rbh_directory: %s\n' % rbh_file_dir)
    config_yaml_ribbon_handle.write('\n')

    # write out species section
    config_yaml_ribbon_handle.write('species:\n')
    for each_species in species_order_list:
        current_fa    = '%s/%s.fa'    % (fa_in_abspath, each_species)
        current_pep   = '%s/%s.pep'   % (fa_in_abspath, each_species)
        current_chrom = '%s/%s.chrom' % (fa_in_abspath, each_species)
        config_yaml_ribbon_handle.write('  %s:\n' % each_species)
        config_yaml_ribbon_handle.write('    proteins:    %s\n' % current_pep)
        config_yaml_ribbon_handle.write('    chrom:       %s\n' % current_chrom)
        config_yaml_ribbon_handle.write('    genome:      %s\n' % current_fa)
        config_yaml_ribbon_handle.write('    minscafsize: %s\n' % minscafsize)
    config_yaml_ribbon_handle.close()

    ################################################## run snakemake ###################################################

    snakemake_cmd_macrosynteny  = 'snakemake --cores %s --snakefile %s' % (thread_num, odp_exe)
    snakemake_cmd_ribbon        = 'snakemake --cores %s --snakefile %s' % (thread_num, odp_rbh_to_ribbon_exe)

    # get macrosynteny plot
    print(snakemake_cmd_macrosynteny)
    os.chdir(get_macrosynteny_plot_wd_abspath)
    os.system(snakemake_cmd_macrosynteny)

    # get ribbon diagram
    print(snakemake_cmd_ribbon)
    os.chdir(get_ribbon_diagram_wd_abspath)
    os.system(snakemake_cmd_ribbon)

    print('Done!')


if __name__ == '__main__':

    ribbon_parser = argparse.ArgumentParser(usage=ribbon_usage)
    ribbon_parser.add_argument('-o',                    required=True,                                  help='output directory')
    ribbon_parser.add_argument('-fa',                   required=True,                                  help='fa file directory')
    ribbon_parser.add_argument('-fax',                  required=False, default='fa',                   help='fa file extension, default: fa')
    ribbon_parser.add_argument('-gff',                  required=False,                                 help='gff file directory')
    ribbon_parser.add_argument('-gffx',                 required=False, default='gff',                  help='gff file extension, default: gff')
    ribbon_parser.add_argument('-gbk',                  required=False,                                 help='gbk file directory')
    ribbon_parser.add_argument('-gbkx',                 required=False, default='gbk',                  help='gbk file extension, default: gbk')
    ribbon_parser.add_argument('-pep',                  required=False,                                 help='pep file directory')
    ribbon_parser.add_argument('-pepx',                 required=False, default='pep',                  help='pep file extension, default: pep')
    ribbon_parser.add_argument('-chrom',                required=False,                                 help='chrom file directory')
    ribbon_parser.add_argument('-chromx',               required=False, default='chrom',                help='chrom file extension, default: chrom')
    ribbon_parser.add_argument('-so',                   required=False, default=None,                   help='species order in the ribbon diagram')
    ribbon_parser.add_argument('-co',                   required=False, default=None,                   help='chromosome order in the ribbon diagram')
    ribbon_parser.add_argument('-m',                    required=False, type=int, default=10,           help='minscafsize, default: 10')
    ribbon_parser.add_argument('-t',                    required=False, type=int, default=1,            help='number of core, default: 1')
    ribbon_parser.add_argument('-f',                    required=False, action="store_true",            help='force overwrite')
    ribbon_parser.add_argument('-plot_all',             required=False, action="store_true",            help='plot_all')
    ribbon_parser.add_argument('-odp',                  required=False, default='odp',                  help='executable file odp, default: odp')
    ribbon_parser.add_argument('-odp_rbh_to_ribbon',    required=False, default='odp_rbh_to_ribbon',    help='executable file odp_rbh_to_ribbon, default: odp_rbh_to_ribbon')
    args = vars(ribbon_parser.parse_args())
    ribbon(args)
